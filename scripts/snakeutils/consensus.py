import os
import gzip
import math
from collections import Counter
from multiprocessing import Pool

import numpy as np

from pss.fastq import (
    experiment_info_from_fastq,
    iter_fastq_simple,
    ConsensusHeader,
    PHRED_SHIFT
)

# Constants
_ADAPTER_SEQ = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
_ADAPTER_DETLEN = 10
_DNA_ALPHABET_ARR = np.array(list(b'ACGT'))

_CONS_UNIQUE = 0
_CONS_FREQUENCY = 1
_CONS_WT = 2
_CONS_CONS = 3
_CONS_WT_NOISY = 4
_CONS_UNCLEAR = 5

_CONS_LOG_MSGS = {
    _CONS_UNIQUE: 'unique',
    _CONS_FREQUENCY: 'frequency',
    _CONS_WT: 'wt',
    _CONS_CONS: 'consensus',
    _CONS_WT_NOISY: 'wt_noisy',
    _CONS_UNCLEAR: 'unclear'
}


def _count_umis(fnames, umilen):
    """
    Count all UMIs in a list of gzipped FASTQ files.

    Parameters
    ----------
    fnames : list of str
        FASTQ filenames (gzipped).
    umilen : int
        Length of the UMI at the start of each read.

    Returns
    -------
    Counter
        Key: UMI string
        Value: number of occurrences across all files
    """
    umi_counts = Counter()
    for fname in fnames:
        with gzip.open(fname, 'rt') as fq:
            for _, read, _, _ in iter_fastq_simple(fq):
                umi_counts.update([read[:umilen]])
    return umi_counts


def _chunk_list(lst, n):
    """
    Split a list into chunks of size <= n.

    Parameters
    ----------
    lst : list
        A list of elements.
    n : int
        Max chunk size.

    Returns
    -------
    list of sets
        Each set up to length n, covering all elements in lst.
    """
    return [set(lst[i:i+n]) for i in range(0, len(lst), n)]


def _adapter_detect(read, detlen=_ADAPTER_DETLEN, adapter=_ADAPTER_SEQ):
    """
    Search for a partial adapter sequence near the end of 'read'.

    Parameters
    ----------
    read : str
    detlen : int
        The length of the adapter prefix to look for.
    adapter : str
        Adapter sequence.

    Returns
    -------
    int or None
        Position of the detected adapter prefix, or None if not found.
    """
    pos = read.rfind(adapter[:detlen])
    return pos if pos != -1 else None


def _ascii_list(string):
    """
    Convert a string to a list of ASCII (integer) values.

    Parameters
    ----------
    string : str

    Returns
    -------
    list of int
    """
    return list(string.encode('ascii'))


def _decode_ascii(ascii_arr):
    """
    Convert an array of ASCII (integer) values to a string.

    Parameters
    ----------
    ascii_arr : array-like of int

    Returns
    -------
    str
    """
    return ''.join(map(chr, ascii_arr))


def _avg_qual(qual_strings):
    """
    Compute the average Phred quality from multiple quality strings.

    Parameters
    ----------
    qual_strings : list of str
        Each string is a line of ASCII-encoded Phred scores.

    Returns
    -------
    str
        A single quality string representing the average of all input strings.
    """
    # Convert each quality string to an array of integers (subtract PHRED_SHIFT)
    q_arr = np.array([_ascii_list(q) for q in qual_strings]) - PHRED_SHIFT
    # Compute mean across all strings, round to int, re-add PHRED_SHIFT
    avg_q_arr = q_arr.mean(axis=0).astype(int) + PHRED_SHIFT
    return _decode_ascii(avg_q_arr)


def _avg_qual_filtered(qual_strings, reads, target_read):
    """
    Average Phred scores for only those reads that match `target_read`.

    Parameters
    ----------
    qual_strings : list of str
        Quality strings, in 1-1 correspondence with 'reads'.
    reads : list of str
        The actual DNA reads in the same order as 'qual_strings'.
    target_read : str
        We'll compute the average quality only for entries in 'reads' 
        that match this read.

    Returns
    -------
    str
        Average Phred quality string over all reads matching target_read.
    """
    filtered_quals = [
        q for (q, r) in zip(qual_strings, reads) if r == target_read
    ]
    return _avg_qual(filtered_quals)


def _div_zero(a, b):
    """
    Safe division, returning 0 if b==0.

    Returns
    -------
    float
    """
    return a / b if b else 0


class ConsensusStats:
    """
    Helper class to store, write, and read statistics from consensus computations.

    The "counter" attributes (reads_raw, reads_removed, reads_short, etc.)
    are summed across multiple runs if desired. The "identity" attributes
    (flowcell, fastq) identify the run or file.

    Attributes
    ----------
    (See the docstring in the original code for a description
     of each attribute.)
    """
    header_info = '\t'.join([
        'FlowCellID', 'fastqFile', 'rawReads', 'excludedReads',
        'shortReads', 'excludedUidReads', 'usableReads', 'shortReadRate',
        'usableReadRate', 'rawUIDs', 'smallUIDs', 'consensusUIDReads',
        'removedUIDReads', 'smallUIDReads', 'rawConsensus',
        'clearConsensus', 'UnclearConsensus'
    ])
    _counter_attr = [
        'reads_raw',
        'reads_removed',
        'reads_short',
        'umis_Q',
        'umis_removed',
        'umis_raw',
        'umis_small',
        'consensus_Q',
        'consensus_removed',
        'consensus_raw',
        'consensus_small',
        'consensus_clear',
        'consensus_unclear'
    ]
    _identity_attr = ['flowcell', 'fastq']

    def __init__(self):
        self.flowcell = None
        self.fastq = None
        for attr in self._counter_attr:
            setattr(self, attr, 0)

    def __add__(self, other):
        """
        Summation of two ConsensusStats objects. Sums counter attributes,
        preserves the identity attributes from self.
        """
        ret = type(self)()
        for attr in self._identity_attr:
            setattr(ret, attr, getattr(self, attr))
        for attr in self._counter_attr:
            val = getattr(self, attr) + getattr(other, attr)
            setattr(ret, attr, val)
        return ret

    @property
    def short_read_rate(self):
        return _div_zero(self.reads_short, self.reads_raw)

    @property
    def usable_read_rate(self):
        return _div_zero(self.umis_Q, self.reads_raw)

    @property
    def sufficient_uid_size_rate(self):
        return _div_zero(self.consensus_Q, self.umis_Q)

    @property
    def usable_consensus_rate(self):
        return _div_zero(self.consensus_clear, self.consensus_raw)

    def tab_str(self):
        """
        Format the object's data into a tab-separated string.
        Follows the order of 'header_info'.
        """
        return '\t'.join(str(val) for val in (
            self.flowcell,
            self.fastq,
            # reads
            self.reads_raw,
            self.reads_removed,
            self.reads_short,
            # umis
            self.umis_removed,
            self.umis_Q,
            self.short_read_rate,
            self.usable_read_rate,
            self.umis_raw,
            self.umis_small,
            # consensus
            self.consensus_Q,
            self.consensus_removed,
            self.consensus_small,
            self.consensus_raw,
            self.consensus_clear,
            self.consensus_unclear
        ))

    @classmethod
    def from_string(cls, tab_str):
        """
        Reconstruct a ConsensusStats from a tab-delimited line.
        """
        stats = cls()
        entries = tab_str.strip().split('\t')
        stats.flowcell = entries[0]
        stats.fastq = entries[1]
        # reads
        stats.reads_raw = int(entries[2])
        stats.reads_removed = int(entries[3])
        stats.reads_short = int(entries[4])
        # umis
        stats.umis_removed = int(entries[5])
        stats.umis_Q = int(entries[6])
        stats.umis_raw = int(entries[9])
        stats.umis_small = int(entries[10])
        # consensus
        stats.consensus_Q = int(entries[11])
        stats.consensus_removed = int(entries[12])
        stats.consensus_small = int(entries[13])
        stats.consensus_raw = int(entries[14])
        stats.consensus_clear = int(entries[15])
        stats.consensus_unclear = int(entries[16])
        return stats

    @classmethod
    def write_stats(cls, fname, stats_list):
        """
        Write multiple ConsensusStats objects to a file, one line each.
        """
        with open(fname, 'w') as fstats:
            fstats.write(cls.header_info + '\n')
            for stats in stats_list:
                fstats.write(stats.tab_str() + '\n')

    @classmethod
    def read_stats(cls, fname):
        """
        Read a file containing multiple lines of consensus statistics.
        Returns a list of ConsensusStats.
        """
        with open(fname) as fstats:
            next(fstats)  # skip header
            return [cls.from_string(line) for line in fstats]


class ConsensusCalculator:
    """
    Computes consensus sequences from FASTQ reads grouped by UMI (unique molecular identifier).

    Workflow:
      1. For each set of FASTQ(s), gather UMIs above the minimum family size
         (the rest are 'small' families).
      2. For each chunk of UMI families, collect reads, build consensus for each family.
      3. Write the consensus read(s) to a gzipped FASTQ file, plus record stats.

    Parameters
    ----------
    wt_strings : list of str
        A list of "wild type" references. If reads closely match any 
        of these sequences, we treat them as wild type for certain steps.
    min_familysize : int
        Minimum read count required for a UMI to be considered for consensus.
    umi_len : int
        Length (in bases) of the UMI at the start of each read.
    min_umi_qual : int
        Minimum Phred quality for each base of the UMI region.
        Reads with lower quality UMIs are discarded.
    max_umis_per_computation : int, optional
        A limit on how many UMIs we handle in a single chunk (for memory constraints).
    """

    def __init__(
        self,
        wt_strings,
        min_familysize,
        umi_len,
        min_umi_qual,
        max_umis_per_computation=100_000
    ):
        # Wild-type references combined into one '|' separated string
        self._wt_string = '|'.join(wt_strings)
        self._umi_len = umi_len
        self._min_umi_qual = min_umi_qual
        self._min_familysize = min_familysize

        # Constants / thresholds for internal usage
        self._lim = 0.9                # coverage fraction needed for a "dominant" base
        self._min_wildtype_match = 5   # used in certain fallback checks
        self._wt_match_margin = 10     # substring length for wildcard matching
        self._max_n_umis = max_umis_per_computation
        self._min_len_read = 50        # minimum read length after removing UMI & adapter

    # -------------------------------------------------------------------------
    # PHASE 1: READ COLLECTION
    # -------------------------------------------------------------------------
    def _collect_reads(self, fnames_fastq, umi_set):
        """
        From a set of FASTQ files, gather all reads that have a UMI in 'umi_set'.

        This includes checking UMI quality, adapter position, read length, etc.

        Parameters
        ----------
        fnames_fastq : list of str
            Gzipped FASTQ filenames.
        umi_set : set of str
            UMIs of interest.

        Returns
        -------
        (ConsensusStats, dict)
            The consensus stats encountered while parsing reads,
            plus a dict {umi: ([read1, read2, ...], [qual1, qual2, ...])}.
        """
        stats = ConsensusStats()
        # For each UMI, store a tuple: (list_of_reads, list_of_quals)
        umi_reads = {umi: ([], []) for umi in umi_set}

        for fname_fq in fnames_fastq:
            with gzip.open(fname_fq, 'rt') as fq:
                for _, read, _, qual in iter_fastq_simple(fq):
                    stats.reads_raw += 1

                    # Check minimal UMI quality
                    umi_qual_vals = _ascii_list(qual[:self._umi_len])
                    if (min(umi_qual_vals) - PHRED_SHIFT) < self._min_umi_qual:
                        stats.umis_removed += 1
                        continue

                    # Detect adapter to limit read length
                    adapter_start = _adapter_detect(read)
                    umi = read[:self._umi_len]

                    # Truncate read/qual beyond adapter
                    read_seq = read[self._umi_len:adapter_start]
                    qual_seq = qual[self._umi_len:adapter_start]

                    # Check minimal read length (after removing UMI + adapter region)
                    if len(read_seq) < self._min_len_read:
                        stats.reads_short += 1
                        continue

                    # Exclude reads containing 'N'
                    if 'N' in read_seq:
                        stats.reads_removed += 1
                        continue

                    stats.umis_Q += 1

                    # If the read's UMI is in umi_set, store it
                    if umi in umi_set:
                        umi_reads[umi][0].append(read_seq)
                        umi_reads[umi][1].append(qual_seq)

        # For each UMI family, unify read/qual lengths to the min read length in that family
        for umi, (fam_reads, fam_quals) in umi_reads.items():
            if fam_reads:
                minlen = min(len(r) for r in fam_reads)
                fam_reads = [r[:minlen] for r in fam_reads]
                fam_quals = [q[:minlen] for q in fam_quals]
            umi_reads[umi] = (fam_reads, fam_quals)

        return stats, umi_reads

    # -------------------------------------------------------------------------
    # PHASE 2: CONSENSUS COMPUTATION
    # -------------------------------------------------------------------------
    def _calc_consensus(self, reads, quals, lim):
        """
        Build a consensus by examining each base column. If coverage of a 
        single base is >= lim * len(reads), use it. Otherwise, attempt 
        matching to a known wild type substring.

        Parameters
        ----------
        reads : list of str
        quals : list of str
        lim : int
            numeric coverage threshold (or minimal # reads) for a base 
            to be considered "dominant."

        Returns
        -------
        (consensus_seq, consensus_qual) : (str or None, str or None)
            If None, the method was unsuccessful in forming a consensus.
        """
        # Average quality across the entire family
        avg_quality = _avg_qual(quals)

        # Convert reads to a 2D array of ASCII codes
        reads_arr = np.array([_ascii_list(r) for r in reads])
        # For each possible base (A/C/G/T), count how many reads have that base 
        # at each position
        nt_counts = np.array([
            (reads_arr == nt).sum(axis=0) for nt in _DNA_ALPHABET_ARR
        ])
        # The index of the most abundant base at each position
        consens_i = nt_counts.argmax(axis=0)
        # The coverage of the chosen base
        max_counts = nt_counts.max(axis=0)

        # Build the "dominant" consensus array
        consensus_arr = _DNA_ALPHABET_ARR[consens_i]

        # Identify positions that do not meet the coverage threshold (max_counts < lim)
        unclear_positions = np.flatnonzero(max_counts < lim)
        if not unclear_positions.size:
            # Everything is above threshold
            return _decode_ascii(consensus_arr), avg_quality

        # Some positions are below coverage threshold, so try to salvage them
        # by checking partial matches to the known WT strings.
        # If the first unclear position is near the start, we won't have enough 
        # context to do a partial match
        if unclear_positions[0] < self._min_wildtype_match:
            return None, None

        # Build the final consensus in a piecewise manner
        consensus_seq = ''
        for i in unclear_positions:
            # Add up to the unclear position
            consensus_seq += _decode_ascii(consensus_arr[len(consensus_seq):i])

            # Attempt to choose a base among the ones observed
            # sorted by coverage, from most to fewest
            col_counts = nt_counts[:, i]
            valid_mask = col_counts > 0
            # Sort bases in descending order of coverage
            candidate_bases = [
                nt for coverage, nt in 
                sorted(zip(col_counts[valid_mask], _DNA_ALPHABET_ARR[valid_mask]),
                       reverse=True)
            ]
            picked_base = False
            for nt in candidate_bases:
                trial = consensus_seq[-self._wt_match_margin:] + chr(nt)
                if trial in self._wt_string:
                    consensus_seq += chr(nt)
                    picked_base = True
                    break
            if not picked_base:
                # Could not find a suitable base for this position
                return None, None

        # Add any remaining bases to the consensus
        consensus_seq += _decode_ascii(consensus_arr[len(consensus_seq):])
        return consensus_seq, avg_quality

    def _get_consensus(self, reads, quals):
        """
        Compute the consensus for a single UMI family.

        Strategy:
          1) If all reads are identical, that's the consensus.
          2) If the top read is >90% of the family, use that.
          3) If a read matches a known WT string and is sufficiently 
             common, use that.
          4) Otherwise, build a column-by-column consensus. 
             If that fails, pick the best or a "noisy" WT read.

        Returns
        -------
        (consensus_read, consensus_qual, cons_flag)
        """
        read_count = Counter(reads)
        if len(read_count) == 1:
            # all reads are identical
            return reads[0], _avg_qual(quals), _CONS_UNIQUE

        # The top read is the most common read in the set
        read_counts_sorted = read_count.most_common()
        main_read, main_count = read_counts_sorted[0]
        lim_count = math.ceil(len(reads) * self._lim)

        def finalize_result(r, status_flag):
            return r, _avg_qual_filtered(quals, reads, r), status_flag

        # If top read covers >= 90% of the family
        if main_count >= lim_count:
            return finalize_result(main_read, _CONS_FREQUENCY)

        # Check if one of the top few reads is in the WT string 
        # and is "common enough" to overshadow others
        wt_lim = len(reads) - lim_count
        for rd, cnt in read_counts_sorted[:2]:
            if rd in self._wt_string and cnt > wt_lim:
                return finalize_result(rd, _CONS_WT)

        # Attempt a column-based consensus
        consensus_seq, consensus_qual = self._calc_consensus(reads, quals, lim_count)
        if consensus_seq:
            return consensus_seq, consensus_qual, _CONS_CONS

        # Last fallback: if any read is strictly a known WT, pick the most common
        for rd, cnt in read_counts_sorted:
            if rd in self._wt_string:
                return finalize_result(rd, _CONS_WT_NOISY)

        return None, None, _CONS_UNCLEAR

    def _cons_seqs(self, umi_reads):
        """
        Compute consensus sequences for each UMI. Update usage stats.

        Parameters
        ----------
        umi_reads : dict
            Key: UMI
            Value: (list_of_reads, list_of_quals)

        Returns
        -------
        (ConsensusStats, dict)
            The stats for this set of UMI families,
            and a dict {umi: (consensus_seq, consensus_qual, consensus_flag)}.
        """
        stats = ConsensusStats()
        consensus_seqs = {}

        for umi, (reads, quals) in sorted(umi_reads.items()):
            if len(reads) < self._min_familysize:
                stats.umis_small += 1
                stats.consensus_small += len(reads)
                continue

            consensus, cons_qual, cons_flag = self._get_consensus(reads, quals)
            consensus_seqs[umi] = (consensus, cons_qual, cons_flag)

            # Update stats
            stats.consensus_Q += len(reads)
            stats.consensus_raw += 1

            if cons_flag == _CONS_UNCLEAR:
                stats.consensus_unclear += 1
            else:
                stats.consensus_clear += 1

        return stats, consensus_seqs

    # -------------------------------------------------------------------------
    # PHASE 3: OUTPUT WRITING
    # -------------------------------------------------------------------------
    def _write_consensus_seqs(self,
                              umi_reads,
                              consensus_seqs,
                              sequencer,
                              flowcell,
                              index,
                              fname_out):
        """
        Write the final consensus sequences (where available) to a gzipped FASTQ.

        Parameters
        ----------
        umi_reads : dict
            {umi: ([reads], [quals])}
        consensus_seqs : dict
            {umi: (consensus_seq or None, consensus_qual, cons_log)}
        sequencer : str
        flowcell : str
        index : str
        fname_out : str
            The gzipped FASTQ output path.
        """
        with gzip.open(fname_out, 'at') as fq_out:
            for umi, (consensus, qual, cons_flag) in consensus_seqs.items():
                if consensus is None:
                    continue
                umi_log_str = _CONS_LOG_MSGS[cons_flag]
                n_reads = len(umi_reads[umi][0])
                header_obj = ConsensusHeader(
                    sequencer=sequencer,
                    flowcell=flowcell,
                    umi=umi,
                    index=index,
                    umi_log_str=umi_log_str,
                    n_reads=n_reads
                )
                fq_out.write('\n'.join([
                    str(header_obj),
                    consensus,
                    '+',
                    qual,
                    ''
                ]))

    # -------------------------------------------------------------------------
    # Public Methods
    # -------------------------------------------------------------------------
    def compute_consensus_sequences(self, fnames_fastq_in, fname_out):
        """
        Compute consensus sequences for a group of FASTQ files 
        representing reads from a single well (or sample).
        Writes the consensus to a gzipped FASTQ file.

        Parameters
        ----------
        fnames_fastq_in : str or list of str
            Path(s) to input gzipped FASTQ(s).
        fname_out : str
            Path to the gzipped FASTQ for consensus output.

        Returns
        -------
        (str, ConsensusStats)
            The flowcell ID, plus aggregated statistics of the run.
        """
        if isinstance(fnames_fastq_in, str):
            fnames_fastq_in = [fnames_fastq_in]

        # Identify flowcell info from the first FASTQ
        sequencer, flowcell, index = experiment_info_from_fastq(fnames_fastq_in[0])

        # Count all UMIs in the input files, separate "small" families
        umi_counter = _count_umis(fnames_fastq_in, self._umi_len)
        umis_small = {
            umi: ct for umi, ct in umi_counter.items() 
            if ct < self._min_familysize
        }
        # The "big enough" UMIs
        umis = list(umi_counter.keys() - umis_small.keys())
        # Chunk UMIs to limit memory usage
        umis_chunked = _chunk_list(umis, self._max_n_umis)

        # Prepare stats
        stats = ConsensusStats()
        stats.flowcell = flowcell
        stats.fastq = os.path.split(fname_out)[1]
        stats.umis_raw = len(umi_counter)
        # Families below the size threshold
        stats.umis_small += len(umis_small)
        stats.consensus_small += sum(umis_small.values())

        # Ensure the output file exists (start fresh)
        with gzip.open(fname_out, 'at') as fq_out:
            fq_out.write('')

        # Process each chunk of UMIs
        for i, umi_set in enumerate(umis_chunked):
            chunk_stats, umi_reads = self._collect_reads(fnames_fastq_in, umi_set)
            # For the first chunk, incorporate read stats
            if i == 0:
                stats += chunk_stats

            # Build consensus
            stats_cons, consensus_seqs = self._cons_seqs(umi_reads)
            stats += stats_cons

            # Write out results
            self._write_consensus_seqs(
                umi_reads, 
                consensus_seqs,
                sequencer, 
                flowcell, 
                index, 
                fname_out
            )

        return flowcell, stats

    def pooled_consensus_computation(self,
                                     fnames_fastq_in,
                                     fnames_cons_out,
                                     fname_stats_out,
                                     cpu_limit):
        """
        Parallelize consensus computation across multiple input files 
        using a process pool.

        Parameters
        ----------
        fnames_fastq_in : list of (str or list of str)
            Each element is either a single FASTQ path or a list of FASTQ paths 
            that belong to one sample/well.
        fnames_cons_out : list of str
            Output FASTQ paths for each element in fnames_fastq_in.
        fname_stats_out : str
            File where aggregated stats are written (tab-delimited).
        cpu_limit : int
            Number of processes to spawn.
        """
        # Map each input group to compute_consensus_sequences
        with Pool(processes=cpu_limit) as pool:
            stats_list = pool.starmap(
                self.compute_consensus_sequences,
                zip(fnames_fastq_in, fnames_cons_out)
            )
        # stats_list is a list of (flowcell, ConsensusStats)
        # We only need the stats objects
        final_stats = [st for (_, st) in stats_list]
        ConsensusStats.write_stats(fname_stats_out, final_stats)



