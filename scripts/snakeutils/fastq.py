"""Handle fastq files

Attributes
----------
PHRED_SHIFT : int
    Shift of ascii encoding to phred score.
    Subtract this from the ascii integer of a phred
    character from the quality string to get the score.

"""
import numpy as np
import gzip
from collections import Counter


PHRED_SHIFT = 33


def phred_array(qual_str):
    """Get an integer array for a Phred score quality string.

    Parameters
    ----------
    phred_str

    Returns
    -------
    array of int

    """
    return np.array(list(qual_str.encode('ascii'))) - PHRED_SHIFT



##############################################################################
# Helper Classes and Functions (previously imported from .fastq)
##############################################################################

class FastqEntry:
    """
    A read entry from a FASTQ file.
    Parameters
    ----------
    header : str
        Header line from the FASTQ file (includes leading '@').
    read : str
        DNA sequence of the read.
    qual : str
        Phred quality string for the read.
    comment : str
        The third line of the FASTQ entry (often '+', but can hold comments).

    Attributes
    ----------
    header_str : str
        Raw header line (includes leading '@').
    header : str
        By default, the same as header_str (subclasses may parse it differently).
    read : str
        DNA sequence.
    qual : str
        Phred quality string.
    comment : str
        The '+' comment line (or other possible line content).
    """
    def __init__(self, header, read, qual, comment):
        self.header_str = header
        self.header = self.header_str
        self.read = read
        self.qual = qual
        self.comment = comment

    def __str__(self):
        return '\n'.join([
            str(self.header),
            self.read,
            self.comment,
            self.qual,  # blank line after
            ''
        ])


class ConsensusHeader:
    """
    A specialized header structure for consensus reads, parsed from lines
    such as "@sequencer:FCID:flowcell:UID:umi:INDEX:index:TYPE:umi_log_str:N:n_reads"

    Parameters
    ----------
    sequencer : str
        Sequencer ID.
    flowcell : str
        Flowcell ID.
    umi : str
        UMI sequence.
    index : str
        Index sequence information.
    umi_log_str : str
        String describing how the consensus was formed (e.g., "unique", "frequency", etc.).
    n_reads : int
        The number of raw reads that contributed to this consensus.
    """
    _header_entries = {
        'FCID': 'flowcell',
        'UID': 'umi',
        'INDEX': 'index',
        'TYPE': 'umi_log_str',
        'N': 'n_reads'
    }
    sep = ':'

    def __init__(self, sequencer, flowcell, umi, index, umi_log_str, n_reads):
        self.sequencer = sequencer
        self.flowcell = flowcell
        self.umi = umi
        self.index = index
        self.umi_log_str = umi_log_str
        self.n_reads = int(n_reads)

    def __str__(self):
        """
        Reconstruct the header as a string (with leading '@').
        """
        entries = ['@' + self.sequencer]
        for label, attr in self._header_entries.items():
            entries.append(label)
            entries.append(str(getattr(self, attr)))
        return self.sep.join(entries)

    @classmethod
    def from_string(cls, header_str):
        """
        Parse a FASTQ consensus header string into a ConsensusHeader object.
        """
        parts = header_str.split(cls.sep)
        # The first part is '@SEQUENCER'
        sequencer = parts[0][1:]  # remove leading '@'
        # The rest are label/value pairs
        kwargs = {
            cls._header_entries[label]: val
            for label, val in zip(parts[1:-1:2], parts[2::2])
        }
        return cls(sequencer, **kwargs)


class ConsensusEntry(FastqEntry):
    """
    A specialized FastqEntry for consensus reads, whose header is parsed
    into a ConsensusHeader instance.
    """
    def __init__(self, header, read, qual, comment):
        super().__init__(header, read, qual, comment)
        self.header = ConsensusHeader.from_string(self.header_str)


def _iter_fastq_simple(line_iterable):
    """
    Minimal FASTQ parser that yields 4-tuple of lines (header, read, comment, qual).
    Expects blocks of 4 lines each. Stops when exhausted.

    Parameters
    ----------
    line_iterable : iterable
        Yields lines of a FASTQ file.

    Yields
    ------
    (header, read, comment, qual) : tuple of str
    """
    lines = iter(line_iterable)
    while True:
        try:
            header = next(lines).strip()
            read = next(lines).strip()
            comment = next(lines).strip()
            qual = next(lines).strip()
            yield header, read, comment, qual
        except StopIteration:
            break  # no more lines


def _iter_fastq(line_iterable, entry_class):
    """
    General FASTQ parser that uses _iter_fastq_simple internally,
    returning objects of the specified entry_class.

    Parameters
    ----------
    line_iterable : iterable
        Yields lines of a FASTQ file.
    entry_class : class
        A class (e.g., FastqEntry or ConsensusEntry) for wrapping each entry.

    Yields
    ------
    entry_class instance
    """
    for head, read, comment, qual in _iter_fastq_simple(line_iterable):
        yield entry_class(head, read, qual, comment)


##############################################################################
# Required Functions (self-contained, no external .fastq imports)
##############################################################################

def experiment_info_from_fastq(fname_fastq):
    """
    Extract high-level experiment information (sequencer, flowcell, and index)
    from the very first read of a gzipped FASTQ file.

    Parameters
    ----------
    fname_fastq : str
        Path to a gzipped FASTQ file.

    Returns
    -------
    sequencer : str
        Name or ID of the sequencer.
    flowcell : str
        Flowcell identifier.
    index : str
        Index (the last field in the header).
    """
    with gzip.open(fname_fastq, 'rt') as fq:
        # We only need the first entry's header
        for header, _, _, _ in _iter_fastq_simple(fq):
            parts = header.strip().split(":")
            # The first part is "@SEQUENCER"
            sequencer = parts[0][1:]  # remove leading '@'
            flowcell = parts[2]
            index = parts[-1]
            return sequencer, flowcell, index
    # If file is empty, returns None implicitly.


def iter_fastq_consensus(line_iterable):
    """
    Iterate over lines from a gzipped FASTQ file containing consensus reads,
    yielding specialized ConsensusEntry objects.

    Parameters
    ----------
    line_iterable : iterable of str
        Lines of a FASTQ file. Each group of 4 lines is one entry.

    Yields
    ------
    ConsensusEntry
        A specialized FASTQ entry that stores consensus metadata in its header.
    """
    return _iter_fastq(line_iterable, entry_class=ConsensusEntry)


def fastq_consensus_count_distinct_reads(fname):
    """
    Count distinct consensus reads in a gzipped FASTQ file. The count stored is
    the sum of the number of raw reads that contributed to each consensus
    (parsed from each entry's header).

    Parameters
    ----------
    fname : str
        Path to a gzipped FASTQ file containing consensus reads.

    Returns
    -------
    Counter
        A mapping where the key is the consensus read (string),
        and the value is the total number of raw reads that contributed
        to that consensus read.
    """
    read_counter = Counter()
    with gzip.open(fname, 'rt') as f:
        # Use our local iter_fastq_consensus
        for cons_entry in iter_fastq_consensus(f):
            # Each ConsensusEntry has .read (sequence) and .header.n_reads (# raw reads)
            read_counter.update({cons_entry.read: cons_entry.header.n_reads})
    return read_counter


