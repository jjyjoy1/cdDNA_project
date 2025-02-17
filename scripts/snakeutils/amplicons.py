"""Classes for amplicon sequences and information.

Use the class method AmpliconGroup.from_files to create an
AmpliconGroup object. You can get the Amplicon objects from
the group using the __iter__ function or access them via the
Amplicon-IDs.

"""
import numpy as np
import os
import functools
from pss.seqtools import convert_sense
from pss.utils import iter_file_nocomment


def _split_genomic_coordinates(coord_str):
    """
    Parse a genome coordinate string into separate parts.

    The coordinate string typically looks like "GRCh37:5:102-419(:-1)?".
    Example: "GRCh37:5:102-419" or "GRCh38:chr5:100-400:-1"

    Parameters
    ----------
    coord_str : str
        The genomic coordinates in the format: genome:chrom:start-end[:orientation]

    Returns
    -------
    (genome, chrom, start, end, reverse_flag) : (str, str, int, int, bool)
        genome : e.g. "GRCh37"
        chrom : e.g. "5"
        start : the integer start position (1-based in data, but we often shift to 0-based).
        end : the integer end position (1-based).
        reverse_flag : True if orientation is "-1", otherwise False.
    """
    parts = coord_str.split(':')
    genome = parts[0]
    chrom = parts[1]
    start_str, end_str = parts[2].split('-')
    start, end = int(start_str), int(end_str)

    # If orientation is given, parse it. Otherwise default to forward (False).
    reverse_flag = False
    if len(parts) > 3:
        reverse_flag = (parts[3] == '-1')

    return genome, chrom, start, end, reverse_flag


def _stringify_orientation(forward):
    """
    Return '1' if forward is True, otherwise '-1'.

    Parameters
    ----------
    forward : bool

    Returns
    -------
    str
        '1' for forward, '-1' for reverse.
    """
    return '1' if forward else '-1'


def _get_forward_reverse(seq, forward):
    """
    If forward is True, return seq as is. 
    If forward is False, return its reverse complement via convert_sense().

    Parameters
    ----------
    seq : str
        The DNA sequence.
    forward : bool
        Whether to keep the sequence orientation or reverse-complement it.

    Returns
    -------
    str
    """
    return seq if forward else convert_sense(seq)


@functools.total_ordering
class Amplicon:
    """
    Represents a single amplicon with primer sequences, an amplicon body, 
    and various metadata (genome coordinates, gene symbol, transcript, etc.).

    Parameters
    ----------
    amplicon_id : str
        The unique identifier for this amplicon.
    forward_primer : str
        Nucleotide sequence of the forward primer.
    reverse_primer : str
        Nucleotide sequence of the reverse primer (provided in forward orientation).
    ampseq : str
        The amplicon body (i.e., the in-seq) located between the two primers.
        May contain intronic regions in lowercase, exonic in uppercase.
    genome_coordinates : str
        A string describing the genomic position, e.g. "GRCh37:5:102-419:-1".
    gene_symbol : str
        The gene symbol, e.g. "BRAF".
    transcript : str
        Transcript identifier or name.
    cds_start : int
        The start position of the CDS region within the amplicon (1-based).
    cds_end : int
        The end position of the CDS region (1-based).
    gene_orientation : str
        '1' or '-1', describing the gene's orientation. 
        If '-1', gene is on the reverse strand.
    uid_position : int or None
        If present, the start position of a UID region in the amplicon.
    error : int
        The default alignment error threshold to allow for read alignment.

    Notes
    -----
    - `casesensitive_ampseq` and `casesensitive_construct` maintain 
      the original uppercase/lowercase exonic/intronic markings.
    - The final stored `_ampseq` and `_construct` are uppercase-only 
      for uniform alignment usage.
    - The genome coordinates are adjusted to 0-based internally for 
      more convenient indexing. 
    """

    def __init__(
        self,
        amplicon_id,
        forward_primer,
        reverse_primer,
        ampseq,
        genome_coordinates,
        gene_symbol,
        transcript,
        cds_start,
        cds_end,
        gene_orientation,
        uid_position,
        error
    ):
        # Amplicon ID
        self._amplicon_id = amplicon_id.upper()

        # Store primers and convert them to uppercase
        self._fp = forward_primer.upper()
        self._rp = reverse_primer.upper()

        # Build the "reverse primer" in forward orientation 
        # (which is the reverse complement of the original _rp)
        rp_converted = convert_sense(self._rp)

        # Save the amplicon sequence (with original case) 
        self._casesensitive_ampseq = ampseq
        # Build a "construct" = forward_primer + ampseq + (reverse primer reversed)
        self._casesensitive_construct = self._fp + ampseq + rp_converted

        # Also store uppercase versions for alignment usage
        self._ampseq = self._casesensitive_ampseq.upper()
        self._construct = self._casesensitive_construct.upper()

        # Parse & store genomic coordinates
        self._genome, self._chromosome, self._genomic_start, self._genomic_end, self._reverse = \
            _split_genomic_coordinates(genome_coordinates)

        # Convert internal start from 1-based to 0-based
        self._genomic_start -= 1

        # Gene orientation (e.g. '1' or '-1')
        self._gene_reverse = (gene_orientation == '-1')
        self._gene_symbol = gene_symbol
        self._transcript = transcript
        self._cds_start = cds_start
        self._cds_end = cds_end

        # UID region (if any)
        self._uid_position = uid_position if uid_position else None

        # Default alignment error threshold
        self._err = error

        # For optional usage (protein seq, cDNA, etc.)
        self._prot_seq = None
        self._cdna = None
        self._ensembl_transcript = None

        # Template dictionary: {template_id: (sequence, error_threshold)}
        self._templates = {}

        # Build arrays describing exonic/intronic offsets
        self._transcript_positions, self._exon_shift = self._compute_transcript_positions()

    ############################################################################
    # Internal Calculation: Exon Mapping
    ############################################################################

    def _compute_transcript_positions(self):
        """
        Construct arrays describing each position's transcript coordinate 
        (if exonic) and intron shift.

        Returns
        -------
        (exon_pos, exon_shift) : (np.ndarray, np.ndarray)
            exon_pos[i] = the cDNA coordinate for self._casesensitive_ampseq[i]
                           or -1 if intronic
            exon_shift[i] = how many bases offset from the exon position 
                            (0 if exonic, +/- if intronic).
        """
        seq_len = len(self._casesensitive_ampseq)
        exon_pos = -np.ones(seq_len, dtype=int)
        exon_shift = np.zeros(seq_len, dtype=float)

        # Mark uppercase positions as exonic
        is_exon = np.array([nt.isupper() for nt in self._casesensitive_ampseq])
        # For exonic positions, define a cDNA coordinate
        # cDNA start is (cds_start - 1) because cds_start is 1-based
        exon_pos[is_exon] = np.arange(is_exon.sum()) + (self._cds_start - 1)

        # Fill intron offsets by scanning from right to left
        for i in reversed(range(seq_len - 1)):
            if exon_pos[i] < 0 and exon_pos[i + 1] >= 0:
                exon_pos[i] = exon_pos[i + 1]
                exon_shift[i] = exon_shift[i + 1] - 1

        # Fill intron offsets by scanning from left to right
        for i in range(seq_len - 1):
            if exon_pos[i + 1] < 0 and exon_pos[i] >= 0:
                exon_pos[i + 1] = exon_pos[i]
                exon_shift[i + 1] = exon_shift[i] + 1

        return exon_pos, exon_shift

    ############################################################################
    # Comparisons & String Representations
    ############################################################################

    def __eq__(self, other):
        """
        Two Amplicons are equal if they share the same ID.
        """
        return self._amplicon_id == other._amplicon_id

    def __lt__(self, other):
        """
        Order Amplicons lexicographically by their ID.
        """
        return self._amplicon_id < other._amplicon_id

    def __repr__(self):
        """
        A tab-separated representation including:
        ID, forward primer, in-seq, reverse primer (converted), 
        genome coordinates, orientation, start, end, error threshold.
        """
        return "\t".join([
            self._amplicon_id,
            self._fp,
            self._casesensitive_ampseq,
            convert_sense(self._rp),
            self.get_genome_coord_str(),
            self.geneorient_str,
            str(self._genomic_start),
            str(self._genomic_end),
            str(self._err)
        ])

    ############################################################################
    # Public Properties
    ############################################################################

    @property
    def id(self):
        """str: Amplicon ID (all uppercase)."""
        return self._amplicon_id

    @property
    def forward_primer(self):
        """str: Forward primer sequence (uppercase)."""
        return self._fp

    @property
    def reverse_primer(self):
        """str: Reverse primer sequence (uppercase). 
        (Note this is not the "reverse complement", just the 
         original reverse primer in forward orientation.)"""
        return self._rp

    @property
    def both_primers(self):
        """
        tuple(str, str): (forward_primer, reverse_primer).
        """
        return (self._fp, self._rp)

    @property
    def chromosome(self):
        """str: Chromosome identifier, e.g. '5', 'X', etc."""
        return self._chromosome

    @property
    def is_forward(self):
        """
        bool: True if the amplicon is on the forward strand 
        (per `_split_genomic_coordinates`).
        """
        return not self._reverse

    @property
    def gene_is_forward(self):
        """
        bool: True if the gene itself is on the forward strand 
        (based on `gene_orientation`).
        """
        return not self._gene_reverse

    @property
    def same_orientation_as_gene(self):
        """
        bool: True if the amplicon and the gene share orientation 
        (both forward or both reverse).
        """
        return (self._reverse == self._gene_reverse)

    @property
    def geneorient_str(self):
        """
        str: '1' if the gene is forward, '-1' if the gene is reverse.
        """
        return _stringify_orientation(not self._gene_reverse)

    @property
    def casesensitive_ampseq(self):
        """
        str: The raw in-seq (ampseq) exactly as provided, preserving 
        exons in uppercase and introns in lowercase.
        """
        return self._casesensitive_ampseq

    @property
    def casesensitive_construct(self):
        """
        str: Forward primer + raw in-seq + reverse primer (converted),
        preserving the original uppercase/lowercase for the in-seq.
        """
        return self._casesensitive_construct

    @property
    def error(self):
        """int: Default alignment error threshold."""
        return self._err

    @property
    def gene_symbol(self):
        """str: Gene symbol, e.g. 'BRAF'."""
        return self._gene_symbol

    @property
    def transcript(self):
        """str: The amplicon's transcript name or identifier."""
        return self._transcript

    @property
    def ensembl_transcript(self):
        """str: The Ensembl transcript ID, if set (defaults to None)."""
        return self._ensembl_transcript

    @property
    def gene_cds(self):
        """
        str or None: The entire coding DNA sequence of the gene (if available).
        """
        return self._cdna

    @property
    def cds_start(self):
        """int: 1-based start position of the CDS region inside the in-seq."""
        return self._cds_start

    @property
    def cds_end(self):
        """int: 1-based end position of the CDS region inside the in-seq."""
        return self._cds_end

    ############################################################################
    # Genome Coordinates
    ############################################################################

    @property
    def genomic_start_db(self):
        """
        int: 1-based genomic start position (converted from the stored 0-based start).
        """
        return self._genomic_start + 1

    @property
    def genomic_start_ampseq_db(self):
        """
        int: 1-based genomic position of where the in-seq (ampseq) starts 
        (skipping forward_primer).
        """
        return self.genomic_start_db + len(self.get_primer(forward=self.is_forward))

    @property
    def genomic_end_db(self):
        """int: 1-based genomic end position (as originally read)."""
        return self._genomic_end

    @property
    def genomic_range_db(self):
        """
        (int, int): The (start, end) positions in 1-based genomic coords.
        """
        return (self.genomic_start_db, self.genomic_end_db)

    @property
    def genomic_start_pss(self):
        """int: 0-based genomic start position."""
        return self._genomic_start

    @property
    def genomic_end_pss(self):
        """int: 0-based genomic end position."""
        return self._genomic_end

    ############################################################################
    # Sequence Access
    ############################################################################

    def get_primer(self, forward=True):
        """
        Retrieve either the forward primer or the reverse primer.

        Parameters
        ----------
        forward : bool
            True => forward primer, 
            False => reverse primer.

        Returns
        -------
        str
        """
        return self._fp if forward else self._rp

    def get_ampseq(self, forward=True):
        """
        The amplicon body (in-seq) in the specified orientation.

        If forward is False, returns its reverse complement.

        Parameters
        ----------
        forward : bool
            Whether to get the forward orientation or the reverse complement.

        Returns
        -------
        str
        """
        return _get_forward_reverse(self._ampseq, forward)

    def get_ampseq_with_margin(self, forward, margin):
        """
        Return the in-seq plus `margin` bases of each primer on both ends.

        The forward primer suffix & reverse primer prefix are included 
        around the amplicon body.

        Parameters
        ----------
        forward : bool
            Orientation for the returned sequence.
        margin : int
            Number of bases from each primer to include.

        Returns
        -------
        str
        """
        seq_with_margin = (
            self._fp[-margin:] + 
            self._ampseq + 
            convert_sense(self._rp)[:margin]
        )
        return _get_forward_reverse(seq_with_margin, forward)

    def get_construct(self, forward=True):
        """
        The full sequence (forward_primer + in-seq + reverse_primer) 
        in the specified orientation.

        Parameters
        ----------
        forward : bool

        Returns
        -------
        str
        """
        return _get_forward_reverse(self._construct, forward)

    def seqlen(self, include_primers=True):
        """
        Length of the amplicon sequence.

        Parameters
        ----------
        include_primers : bool
            Whether to include the primer lengths in the result.

        Returns
        -------
        int
        """
        if include_primers:
            return len(self._construct)
        return len(self._ampseq)

    ############################################################################
    # Template Management
    ############################################################################

    def set_template(self, template_id, seq, error, overwrite=False):
        """
        Add or update a template entry for this amplicon. 
        Also store its reverse complement for potential usage.

        Parameters
        ----------
        template_id : str
            A name/ID for this template.
        seq : str
            The template DNA sequence.
        error : int
            Error threshold for alignment.
        overwrite : bool, default=False
            If True, overwrite any existing template with the same ID.
        """
        add_func = self._templates.__setitem__ if overwrite else self._templates.setdefault
        add_func(template_id, (seq, error))

        # Also store the reverse complement
        rc_key = template_id + "_rc"
        add_func(rc_key, (convert_sense(seq), error))

    def in_templates(self, template_id):
        """
        Check if this amplicon has a template with the given ID.

        Parameters
        ----------
        template_id : str

        Returns
        -------
        bool
        """
        return template_id in self._templates

    def iter_templates(self):
        """
        Iterate over all templates stored in this amplicon.

        Yields
        ------
        (template_id, seq, err) : (str, str, int)
        """
        # Note: original code had an error using `for template_id, (seq, err) in self._templates:`
        # We'll fix it to `.items()`.
        for template_id, (seq, err) in self._templates.items():
            yield template_id, seq, err

    ############################################################################
    # Transcript (Gene) Info
    ############################################################################

    def add_transcript_info(self, transcript_info):
        """
        Assign gene/transcript/protein data from a dictionary. 
        Typically, the dictionary might map gene_symbol -> (ensembl_id, cDNA, prot_seq).

        Parameters
        ----------
        transcript_info : dict of {str : (str, str, str)}
            Maps gene_symbol => (ensembl_transcript, cdna, protein_seq)
        """
        if self._gene_symbol in transcript_info:
            enst, cds, prot = transcript_info[self._gene_symbol]
            self._ensembl_transcript = enst
            self._cdna = cds
            self._prot_seq = prot

    def cds_index(self, amplicon_pos):
        """
        Convert an amplicon-based position to a cDNA coordinate 
        and a shift if it's intronic.

        Parameters
        ----------
        amplicon_pos : int
            0-based offset from the forward primer.

        Returns
        -------
        (cdna_position, shift) : (int or None, float or None)
            cdna_position is the cDNA coordinate (0-based) or None 
            if outside the in-seq. shift is how many bases offset 
            from the exon boundary if intronic.
        """
        # Offset by length of forward primer
        i = amplicon_pos - len(self._fp)
        if i < 0 or i >= len(self._transcript_positions):
            return None, None

        return self._transcript_positions[i], self._exon_shift[i]

    def genomic_index(self, amplicon_pos):
        """
        Translate an amplicon position (0-based from forward primer) 
        into a 1-based genomic coordinate.

        Parameters
        ----------
        amplicon_pos : int

        Returns
        -------
        int
            The 1-based genomic position.
        """
        if self.is_forward:
            return self.genomic_start_db + amplicon_pos
        else:
            # If reverse, the position is measured from the other end
            total_len = len(self._casesensitive_construct)
            return self.genomic_start_db + (total_len - 1 - amplicon_pos)

    def genomic_to_amplicon_index(self, genomic_index):
        """
        Convert a 1-based genomic coordinate into an amplicon index 
        (0-based from forward primer).

        Parameters
        ----------
        genomic_index : int
            A 1-based genomic position.

        Returns
        -------
        int
        """
        if self.is_forward:
            return genomic_index - self.genomic_start_db
        else:
            total_len = len(self._casesensitive_construct)
            return (self.genomic_start_db + total_len - 1) - genomic_index

    def is_exonic_region(self, start, end):
        """
        Check if the range [start, end) in amplicon coordinates 
        is entirely exonic.

        Parameters
        ----------
        start : int
            Start index in the amplicon (0-based from forward primer).
        end : int
            End index (non-inclusive).

        Returns
        -------
        bool
            False if any part of that range is intronic or outside the in-seq.
        """
        lp = len(self._fp)
        start_inseq = start - lp
        end_inseq = end - lp

        if start_inseq < 0 or end_inseq > len(self._exon_shift):
            return False

        # If any intronic shift != 0, it's not purely exonic
        region_shift = self._exon_shift[start_inseq:end_inseq]
        return not np.any(region_shift)

    ############################################################################
    # Coordinate String Generators
    ############################################################################

    def get_genome_coord_str(self, start_shift=0, seq_len=None):
        """
        Generate a string describing the genomic region for part of 
        this amplicon.

        Parameters
        ----------
        start_shift : int, default=0
            Shift from the amplicon's start.
        seq_len : int or None, default=None
            The length of the region. If None, use the entire amplicon length 
            from `_genomic_end`.

        Returns
        -------
        str
            Something like "GRCh37:5:103-420" for a partial region.
        """
        start = self.genomic_start_db + start_shift
        if seq_len is None:
            end = self._genomic_end
        else:
            end = start + seq_len
        return f"{self._genome}:{self._chromosome}:{start}-{end}"

    def get_gene_coord_str(self):
        """
        Describe the gene region covered by this amplicon 
        in cDNA coordinates.

        Returns
        -------
        str
            Example: "BRAF:c.100_500"
        """
        return f"{self._gene_symbol}:c.{self._cds_start}_{self._cds_end}"



def _read_template_line(template_line):
    """Read the TEMPLATE line in an amplicon file

    Parameters
    ----------
    template_line : str
        Line read from amplicon which contains the keyword 'TEMPLATE'.

    Return
    ------
    str
        Unique name of the template.
    str
        DNA sequence.
    int
        Alignment error threshold.

    """
    identifier, seq, error = template_line.upper().split(",")[:3]
    mother = identifier.split("_TEMPLATE")[0]
    template_id = identifier.split(mother)[1][1:]
    return template_id, mother, seq[:97], int(error)


def _read_amplicons_csv(fname_amplicons,
                        err_thresh,
                        sep,
                        assay_prefix):
    """Import a given CSVamplicon file into a dictionary holding the
    amplicons as a list of chromosome IDs as used by the safeSEQ analysis.

    Parameters
    ----------
    fname_amplicons : str
        File path to the amplicon file.
    err_thresh : int
        Alignment erro threshold.
    sep : str
        Separator character used in the CSV file.
    assay_prefix : str
        Prefix used to identify the SafeSEQ assay panel.

    Return
    ------
    dict
        Collection of amplicons indexed by their id.
    list
        Unique list of chromosomes.

    """
    assay_prefix = assay_prefix.upper() if assay_prefix else ""
    amplicons = {}
    with open(fname_amplicons) as f_in:
        for line in iter_file_nocomment(f_in):
            if "TEMPLATE" in line:
                template_id, mother, seq, error = _read_template_line(line)
                amplicons[mother].set_template(template_id, seq, error)
                continue
            line_split = line.split(sep)
            (
                amplicon_id,
                forward_primer,
                reverse_primer,
                ampseq,
                coordinates,
                gene_symbol,
                transcript,
                cds_start,
                cds_end,
                gene_orientation
            ) = line_split[:10]
            try:
                uid_position = line_split[10]
                error = line_split[11]
            except IndexError:
                uid_position = None
                error = err_thresh
            amplicon_id = amplicon_id.upper()
            if assay_prefix in amplicon_id:
                amp = Amplicon(
                    amplicon_id,
                    forward_primer,
                    reverse_primer,
                    ampseq,
                    coordinates,
                    gene_symbol,
                    transcript,
                    int(cds_start),
                    int(cds_end),
                    gene_orientation,
                    uid_position,
                    int(error)
                )
                amplicons.setdefault(amplicon_id, amp)
    return amplicons


def _iter_strings(str_or_list):
    """Iterate a list of strings. Also works if only a string is passed.

    Parameters
    ----------
    str_or_list : str or list of str
         One or more arbitrary strings.

    Yields
    ------
    str
        The stings passed.

    """
    if str_or_list is None:
        return
    if isinstance(str_or_list, str):
        yield str_or_list
    else:
        for string in str_or_list:
            yield string


def _read_transcript_information(transcript_dir):
    """Read gene and protein sequences.

    Parameters
    ----------
    transcript_dir : str
        Path of directory with sequence files.

    Returns
    -------
    dict of str: (str, str, str)
        The transcript-ID is the key. Values are a database name, the
        coding DNA sequence and the amino acid sequence.

    """
    transcript_info = {}
    for fname in os.listdir(transcript_dir):
        if not fname.endswith(".fa"):
            continue
        with open(os.path.join(transcript_dir, fname)) as f:
            firstline = f.readline().strip().split(':')
            tid = firstline[0][1:]
            enst = firstline[1]
            cds = f.readline().strip()
            next(f)  # skipt the second header line
            prot_seq = f.readline().strip()
            transcript_info[tid] = (enst, cds, prot_seq)
    return transcript_info



class AmpliconGroup:
    """
    A container for managing multiple Amplicon objects. 
    Provides methods for reading Amplicons from CSV files, 
    loading template sequences, attaching transcript data, 
    filtering by ID, and more.

    Attributes
    ----------
    _amplicons : dict of {str : Amplicon}
        Maps amplicon IDs to Amplicon instances.
    """

    def __init__(self):
        """
        Create an empty AmpliconGroup. You can populate it by calling 
        `read_amplicons` or `from_files`.
        """
        self._amplicons = {}

    def __getitem__(self, amplicon_id):
        """
        Access a specific Amplicon by its ID.

        Parameters
        ----------
        amplicon_id : str
            The ID of the desired Amplicon.

        Returns
        -------
        Amplicon
            The Amplicon object with the given ID.

        Raises
        ------
        KeyError
            If no Amplicon exists for the specified ID.
        """
        return self._amplicons[amplicon_id]

    def __iter__(self):
        """
        Iterate over all Amplicon objects in this group.

        Yields
        ------
        Amplicon
        """
        return iter(self._amplicons.values())

    def __len__(self):
        """
        Number of Amplicon objects in this group.
        """
        return len(self._amplicons)

    def __bool__(self):
        """
        True if this AmpliconGroup contains any Amplicons, 
        otherwise False.
        """
        return bool(self._amplicons)

    def read_amplicons(self, fname_amplicons, err_thresh, sep, assay_prefix):
        """
        Read Amplicon definitions from a CSV file, adding them 
        to this AmpliconGroup.

        Parameters
        ----------
        fname_amplicons : str
            The path to the amplicon CSV file.
        err_thresh : int
            Default alignment error threshold for any amplicon that 
            doesn't specify one in its CSV row.
        sep : str
            The delimiter used in the CSV file (e.g., ',' or ';').
        assay_prefix : str
            A prefix used to filter Amplicon IDs. Only lines whose 
            Amplicon-ID contains (or matches) this prefix will be loaded.

        Notes
        -----
        - Amplicons from the CSV may also define 'template' entries inline.
        - For each row or template definition, an Amplicon is constructed 
          or updated in `_amplicons`.
        """
        # 
        # This helper is assumed to be defined elsewhere:
        # _read_amplicons_csv(fname_amplicons, err_thresh, sep, assay_prefix)
        #
        new_amplicons = _read_amplicons_csv(
            fname_amplicons, err_thresh, sep, assay_prefix
        )
        # Incorporate new amplicons into this AmpliconGroup
        self._amplicons.update(new_amplicons)

    def read_templates(self, fname_template):
        """
        Read additional template sequences from a file and associate 
        them with existing Amplicons in this group.

        Parameters
        ----------
        fname_template : str
            A CSV (or similarly structured) file containing lines 
            with 'TEMPLATE' definitions.

        Notes
        -----
        - Template lines typically look like 
          'MYAMP_TEMPLATE_01,ACTGACG...,3', where 'MYAMP' references 
          an existing Amplicon ID.
        """
        with open(fname_template, 'r') as f_in:
            for line in iter_file_nocomment(f_in):
                if "TEMPLATE" in line:
                    template_id, mother_id, seq, error = _read_template_line(line)
                    # Attach to the corresponding Amplicon
                    if mother_id in self._amplicons:
                        self._amplicons[mother_id].set_template(template_id, seq, error)

    def add_transcript_info(self, transcript_dir):
        """
        Load transcript/protein data from a directory and attach 
        it to each Amplicon in this group if possible.

        Parameters
        ----------
        transcript_dir : str
            A directory containing *.fa files with transcript 
            and protein sequences.

        Notes
        -----
        - The helper `_read_transcript_information` typically reads 
          each *.fa file, extracting the transcript info, which is 
          then added to each Amplicon if the gene symbol matches.
        """
        transcript_data = _read_transcript_information(transcript_dir)
        for amp in self._amplicons.values():
            amp.add_transcript_info(transcript_data)

    def get_amplicons_filtered(self, amp_id_substring):
        """
        Return a new AmpliconGroup containing only Amplicons 
        whose ID contains a specified substring.

        Parameters
        ----------
        amp_id_substring : str
            The substring to search for in Amplicon IDs.

        Returns
        -------
        AmpliconGroup
            A new AmpliconGroup with a subset of matching Amplicons.
        """
        filtered = AmpliconGroup()
        filtered._amplicons = {
            amp_id: amp
            for amp_id, amp in self._amplicons.items()
            if amp_id_substring in amp_id
        }
        return filtered

    def get_chromosomes(self):
        """
        Collect the unique set of chromosome names for all 
        Amplicons in this group.

        Returns
        -------
        list of str
            Chromosome identifiers, e.g. ["1", "2", "X", ...].
        """
        return list({amp.chromosome for amp in self._amplicons.values()})

    def get_max_primer_length(self):
        """
        Determine the maximum length among all forward and 
        reverse primers in this group.

        Returns
        -------
        int
            The length of the longest primer sequence, or 0 
            if no Amplicons are present.
        """
        if not self._amplicons:
            return 0
        return max(
            max(len(amp.forward_primer), len(amp.reverse_primer)) 
            for amp in self._amplicons.values()
        )

    def get_ids(self):
        """
        List all Amplicon IDs in this group.

        Returns
        -------
        list of str
        """
        return list(self._amplicons.keys())

    def get_all_sequences(self, add_reverse=True):
        """
        Gather the complete sequences (constructs) from all Amplicons. 
        Optionally include their reverse complements.

        Parameters
        ----------
        add_reverse : bool, default=True
            If True, each Amplicon's reverse-complement sequence 
            is also appended to the list.

        Returns
        -------
        list of str
            The sequences (in forward orientation, followed by 
            reverse orientation if add_reverse is True).
        """
        seqs = [amp.get_construct(forward=True) for amp in self._amplicons.values()]
        if add_reverse:
            seqs += [amp.get_construct(forward=False) for amp in self._amplicons.values()]
        return seqs

    def contains_id(self, amp_id):
        """
        Check if a given amplicon ID exists in this group.

        Parameters
        ----------
        amp_id : str
            The Amplicon ID to look for.

        Returns
        -------
        bool
        """
        return amp_id in self._amplicons

    @classmethod
    def from_files(cls,
                   fnames_amplicons,
                   fnames_templates=None,
                   transcript_dir=None,
                   err_thresh=20,
                   sep=',',
                   assay_prefix=None):
        """
        Create an AmpliconGroup by reading amplicon definitions from 
        one or more CSV files, plus optional template and transcript data.

        Parameters
        ----------
        fnames_amplicons : str or list of str
            One or more file paths containing amplicon definitions.
        fnames_templates : str or list of str, optional
            One or more file paths with additional 'TEMPLATE' definitions.
        transcript_dir : str, optional
            If provided, the directory from which to load transcript 
            and protein sequences.
        err_thresh : int, default=20
            Default alignment error threshold used when an amplicon line 
            doesn't specify its own threshold.
        sep : str, default=','
            The CSV delimiter used in the amplicon files.
        assay_prefix : str, optional
            A string used to filter Amplicon IDs (only lines with 
            an ID matching this prefix are loaded).

        Returns
        -------
        AmpliconGroup
            A newly constructed group containing the parsed Amplicons.

        Notes
        -----
        - After reading amplicons, if `fnames_templates` is provided, 
          each template file is processed to attach templates to 
          existing Amplicons.
        - If `transcript_dir` is given, the group tries to load 
          transcript data from that directory for each amplicon.
        """
        # Prepare an empty AmpliconGroup
        group = cls()

        # Convert single strings to lists if needed
        for fname_amp in _iter_strings(fnames_amplicons):
            group.read_amplicons(fname_amp, err_thresh, sep, assay_prefix)

        # If template files were provided, read them
        for fname_temp in _iter_strings(fnames_templates):
            group.read_templates(fname_temp)

        # If a transcript directory is provided, attach transcript data
        if transcript_dir:
            group.add_transcript_info(transcript_dir)

        return group


