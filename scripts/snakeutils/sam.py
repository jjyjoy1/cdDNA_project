import os
import subprocess
import operator
import importlib.resources

from pss.seqtools import convert_sense
import pss.samtools_win

# Constants for indexing fields in a SAM record
_IDX_QNAME = 0
_IDX_FLAG = 1
_IDX_RNAME = 2
_IDX_POS = 3
_IDX_MAPQ = 4
_IDX_CIGAR = 5
_IDX_RNEXT = 6
_IDX_PNEXT = 7
_IDX_TLEN = 8
_IDX_SEQ = 9
_IDX_QUAL = 10

# On Windows, we use the bundled samtools.exe
if os.name == 'nt':
    with importlib.resources.path(pss.samtools_win, 'samtools.exe') as path:
        _SAM_TOOLS_CALL = str(path)
else:
    _SAM_TOOLS_CALL = 'samtools'


class SamError(Exception):
    """
    Custom exception for errors raised within the Sam class.
    """
    pass


class Sam:
    """
    A convenience class for building and writing SAM/BAM alignments.
    It collects alignments in memory (via add() calls), then writes them out
    to SAM or BAM using samtools.

    Attributes
    ----------
    _alignments : list of list
        Each sub-list represents one SAM record with
        [QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL, RG_TAG].
    _refs_len : dict of {str : int}
        Maps reference (e.g. chromosome) names to their lengths.
    _sample_ids : dict of {str : str}
        Maps read group IDs (RG) to sample IDs (SM).
    _sorted : bool
        Whether the in-memory alignments are currently sorted by
        (reference, position, QNAME). Marked False after each new add().
    """

    def __init__(self):
        """
        Initialize an empty Sam object with no references, read groups, or alignments.
        """
        self._alignments = []
        self._refs_len = {}
        self._sample_ids = {}
        self._sorted = False

    # -------------------------------------------------------------------------
    # Methods for adding references and read groups
    # -------------------------------------------------------------------------

    def add_read_group(self, read_group, sample_id):
        """
        Register a read group (RG) along with its sample ID (SM).

        Parameters
        ----------
        read_group : str
            The read group identifier (RG:Z:???).
        sample_id : str
            The sample identifier (SM:???).

        Raises
        ------
        SamError
            If the read_group was previously registered with a different sample ID.
        """
        existing_sample = self._sample_ids.setdefault(read_group, sample_id)
        if existing_sample != sample_id:
            raise SamError(
                f"Read group '{read_group}' already exists with sample '{existing_sample}', "
                f"but tried to add it again with '{sample_id}'"
            )

    def add_reference(self, ref_name, length):
        """
        Register a reference (e.g. chromosome name) and its length.

        Parameters
        ----------
        ref_name : str
            The reference identifier (e.g. 'chr1').
        length : int
            The length of this reference.

        Raises
        ------
        SamError
            If the reference was previously registered with a different length.
        """
        existing_len = self._refs_len.setdefault(ref_name, length)
        if existing_len != length:
            raise SamError(
                f"Reference '{ref_name}' is already set to length {existing_len}, "
                f"cannot redefine with length {length}"
            )

    # -------------------------------------------------------------------------
    # Internal checks
    # -------------------------------------------------------------------------

    def _assert_in(self, dict_label, dictionary, key):
        """
        Ensure that 'key' is present in 'dictionary'. Otherwise, raise SamError.

        Parameters
        ----------
        dict_label : str
            Label for the dictionary (e.g. 'references' or 'groups').
        dictionary : dict
        key : str
            The key to check for.

        Raises
        ------
        SamError
        """
        if key not in dictionary:
            raise SamError(f"'{key}' not found in Sam {dict_label}")

    def _assert_reference(self, ref_name):
        self._assert_in('references', self._refs_len, ref_name)

    def _assert_group(self, group_name):
        self._assert_in('groups', self._sample_ids, group_name)

    # -------------------------------------------------------------------------
    # Internal alignment storage
    # -------------------------------------------------------------------------

    def _add_align(self,
                   qname,
                   flag,
                   rname,
                   pos,
                   mapq,
                   cigar,
                   rnext,
                   pnext,
                   tlen,
                   seq,
                   qual,
                   read_group):
        """
        Low-level method to store a single alignment record.
        All references and read groups must already be registered.

        Parameters correspond directly to SAM fields:
          - qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
          - read_group is appended as a "RG:Z:..." tag
        """
        self._assert_reference(rname)
        self._assert_group(read_group)
        # Build a sub-list for the SAM record
        record = [
            qname,            # QNAME
            flag,             # FLAG
            rname,            # RNAME
            pos,              # POS
            mapq,             # MAPQ
            cigar,            # CIGAR
            rnext,            # RNEXT
            pnext,            # PNEXT
            tlen,             # TLEN
            seq,              # SEQ
            qual,             # QUAL
            f'RG:Z:{read_group}'  # read-group tag
        ]
        self._alignments.append(record)
        self._sorted = False

    def add(self,
            ref_name,
            query_name,
            seq_query,
            qual,
            cigar_list,
            reverse,
            pos,
            read_group='rg001'):
        """
        Public method to create and store a single alignment record.

        This is the primary user-facing method: specify the reference name,
        query name, read sequence/quality, CIGAR, orientation, position,
        and optionally the read group ID.

        Parameters
        ----------
        ref_name : str
            Reference name (e.g. chromosome).
        query_name : str
            Query/Read identifier (QNAME).
        seq_query : str
            The read's nucleotide sequence.
        qual : str
            ASCII-encoded Phred quality string for seq_query.
        cigar_list : list of str
            e.g. ['10=', '2I', '3=']. Will be concatenated into '10=2I3='.
        reverse : bool
            If True, the read is aligned on the reverse strand, so we reverse-complement
            the sequence and reverse the CIGAR segments.
        pos : int
            1-based alignment start (POS) on the reference.
        read_group : str, default 'rg001'
            Which read group this read belongs to.

        Raises
        ------
        SamError
            If the reference or read group has not been registered beforehand.
        """
        # Compute the bitwise FLAG
        flag = 0
        # For a reversed read, set 0x10 (16 decimal)
        if reverse:
            flag |= 16

            # Reverse-complement the sequence
            seq_query = convert_sense(seq_query, replace_non_dna=False)
            # Reverse the quality string
            qual = qual[::-1]
            # Reverse the CIGAR segments
            cigar_list = cigar_list[::-1]

        # Concatenate the CIGAR
        cigar_str = ''.join(cigar_list)

        # Add the alignment record
        self._add_align(
            qname=query_name,
            flag=flag,
            rname=ref_name,
            pos=pos,
            mapq=255,         # Hardcode a default MAPQ
            cigar=cigar_str,
            rnext='*',        # No multi-segment references
            pnext=0,
            tlen=0,
            seq=seq_query,
            qual=qual,
            read_group=read_group
        )

    # -------------------------------------------------------------------------
    # Sorting & Header Generation
    # -------------------------------------------------------------------------

    def _sort_alignments(self):
        """
        Sort the in-memory alignment records by (RNAME, POS, QNAME).
        """
        if not self._sorted:
            self._alignments.sort(
                key=operator.itemgetter(_IDX_RNAME, _IDX_POS, _IDX_QNAME)
            )
            self._sorted = True

    def make_header_str(self):
        """
        Build the SAM header lines, including @RG for each read group
        and @SQ for each reference sequence.

        Returns
        -------
        str
            Multi-line string representing the SAM header.
        """
        lines = []
        # Build read-group lines
        for rg, sm in sorted(self._sample_ids.items()):
            lines.append(f'@RG\tID:{rg}\tPL:Illumina\tSM:{sm}')
        # Build reference lines
        for ref_name, ref_len in sorted(self._refs_len.items()):
            lines.append(f'@SQ\tSN:{ref_name}\tLN:{ref_len}')
        return '\n'.join(lines) + '\n'

    # -------------------------------------------------------------------------
    # Writing to disk
    # -------------------------------------------------------------------------

    def _write_file(self, fname, bam):
        """
        Write all stored alignments to a SAM or BAM file using samtools.

        Parameters
        ----------
        fname : str
            Output filename (e.g. 'output.sam' or 'output.bam').
        bam : bool
            If True, write BAM. Otherwise, write SAM.

        Raises
        ------
        OSError, CalledProcessError
            If samtools invocation fails or file I/O fails.
        """
        self._sort_alignments()

        # Build the samtools call
        call = [_SAM_TOOLS_CALL, 'view', '--with-header']
        if bam:
            call.append('-b')

        mode = 'wb' if bam else 'w'
        with open(fname, mode) as f_out:
            # Create a subprocess for samtools
            with subprocess.Popen(
                call,
                stdin=subprocess.PIPE,
                stdout=f_out,
                stderr=f_out
            ) as samtools_proc:
                # Write header
                header_str = self.make_header_str()
                samtools_proc.stdin.write(header_str.encode())

                # Write each alignment line
                for align in self._alignments:
                    line_str = '\t'.join(str(x) for x in align) + '\n'
                    samtools_proc.stdin.write(line_str.encode())

                samtools_proc.communicate()

    def write_sam(self, fname):
        """
        Write the current alignments as a SAM file.

        Parameters
        ----------
        fname : str
            Path for the SAM output file.
        """
        self._write_file(fname, bam=False)

    def write_bam(self, fname):
        """
        Write the current alignments as a BAM file, and index it via samtools.

        Parameters
        ----------
        fname : str
            Path for the BAM output file.
        """
        self._write_file(fname, bam=True)
        # After writing the BAM, index it
        subprocess.run([_SAM_TOOLS_CALL, 'index', fname], check=True)



