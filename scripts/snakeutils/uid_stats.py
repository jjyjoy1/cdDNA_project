import numpy as np
from collections import defaultdict, Counter
from pss.utils import nested_defaultdict


class UidStats:
    """
    Tracks counts of reads in UMI (UID) families across multiple wells and amplicons,
    including quantispike reads.

    Attributes
    ----------
    _fc_dist : defaultdict of {str : Counter}
        For each amplicon (key: amp_id), a Counter storing how many UID families
        have a particular read count.
    _readcounts : defaultdict of {str : list of int}
        For each amplicon (key: amp_id), store a list of read counts for each UID family.
    _amp_ids : set of str
        A set of all non-quantispike amplicon IDs encountered.
    _qs_amp_id_map : dict of {str : str}
        Maps a quantispike ID to its target amplicon ID (the aligned amplicon).
    _readcounts_per_well : nested defaultdict
        2-layer dict keyed by [fastq_filename][amp_id] -> list of read counts.
    """

    def __init__(self):
        """
        Initialize empty data structures.
        """
        # For each amp_id, store a Counter that maps "read_count" -> number_of_families
        self._fc_dist = defaultdict(Counter)

        # For each amp_id, store a list of family sizes
        self._readcounts = defaultdict(list)

        # A set of normal (non-quantispike) amplicon IDs
        self._amp_ids = set()

        # If a read is identified as quantispike, map the QS ID to its amplicon
        self._qs_amp_id_map = {}

        # 2D structure for well-level data: [fastq_file][amp_id] -> list of readcounts
        self._readcounts_per_well = nested_defaultdict(2, list)

    def _add_uid_fam(self, fastq, amp_id, n_reads, qs_id):
        """
        Record one UID family's read count for a specific well and amplicon (or quantispike).

        Parameters
        ----------
        fastq : str
            The FASTQ filename (identifying a well).
        amp_id : str
            Amplicon ID or "None" if we are dealing with quantispike
            reads (in which case we store the QS ID in qs_id).
        n_reads : int
            The size of this UID family (number of reads).
        qs_id : str or None
            If present, indicates a quantispike ID. In that case we record the relationship
            between quantispike and the amplicon ID in `_qs_amp_id_map`.
        """
        if qs_id:
            # This is a quantispike read
            self._qs_amp_id_map[qs_id] = amp_id
            key_id = qs_id
        else:
            # Normal amplicon read
            self._amp_ids.add(amp_id)
            key_id = amp_id

        # At the well level
        self._readcounts_per_well[fastq][key_id].append(n_reads)

        # At the global level
        self._readcounts[key_id].append(n_reads)
        self._fc_dist[key_id].update([n_reads])

    def add(self, fastq, consensus, aligned_reads):
        """
        Incorporate the read count (n_reads) of a single UID family, based on
        the final consensus entry.

        Parameters
        ----------
        fastq : str
            The FASTQ filename from which the read came (represents a well).
        consensus : ConsensusEntry
            The consensus read object that includes `header.n_reads`.
        aligned_reads : dict
            A mapping of read sequence -> AlignedRead. We'll look up the AlignedRead
            to see if it's a normal amplicon read or a quantispike read.
        """
        aligned_read = aligned_reads[consensus.read]
        # We only record amplicon or quantispike reads
        if not (aligned_read.is_amplicon_read() or aligned_read.is_quantispike):
            return

        self._add_uid_fam(
            fastq,
            aligned_read.amplicon_id,
            consensus.header.n_reads,
            aligned_read.quantispike
        )

    def uid_count(self, amp_id):
        """
        How many UID families were observed for this amplicon (or quantispike ID).

        Parameters
        ----------
        amp_id : str
            Amplicon ID or quantispike ID.

        Returns
        -------
        int
            Number of families (unique UMIs) found for that ID.
        """
        return len(self._readcounts[amp_id])

    def well_uid_count(self, fastq, amp_id):
        """
        Number of UID families for a given amplicon in a specific well (FASTQ).

        Parameters
        ----------
        fastq : str
            The FASTQ filename identifying the well.
        amp_id : str
            Amplicon ID or quantispike ID.

        Returns
        -------
        int
            The count of UID families in that well.
        """
        return len(self._readcounts_per_well[fastq][amp_id])

    def get_amplicon_ids(self, include_quantispike=False):
        """
        Retrieve all amplicon IDs for which UID families were found.
        Optionally include quantispike IDs.

        Parameters
        ----------
        include_quantispike : bool, default=False
            If True, also return quantispike IDs.

        Returns
        -------
        set of str
            Amplicon IDs, plus QS IDs if requested.
        """
        if include_quantispike:
            return self._amp_ids.union(self._qs_amp_id_map)
        return self._amp_ids

    def qs_counts(self):
        """
        If quantispike reads are present, find the pair (quantispike_id, target_amp_id)
        with the largest number of UID families.

        Because multiple quantispikes may exist, this picks the most abundant pair.

        Returns
        -------
        (int, int, str, str) or (None, None, None, None)
            (qs_families, target_families, qs_id, target_amp_id)
            If no quantispikes exist, returns all None.
        """
        if not self._qs_amp_id_map:
            return None, None, None, None

        # Build a list of (qs_fam_count, amp_fam_count, qs_id, amp_id) for each QS
        qs_target_counts = [
            (self.uid_count(qs_id), self.uid_count(amp_id), qs_id, amp_id)
            for qs_id, amp_id in self._qs_amp_id_map.items()
        ]
        # Sort and return the highest
        return sorted(qs_target_counts)[-1]

    def max_readcount(self):
        """
        The maximum read count (i.e. the largest UID family size) among all families.

        Returns
        -------
        int
        """
        return max(
            c
            for read_counts in self._readcounts.values()
            for c in read_counts
        )

    def min_readcount(self):
        """
        The minimum read count (smallest UID family size) among all families.

        Returns
        -------
        int
        """
        return min(
            c
            for read_counts in self._readcounts.values()
            for c in read_counts
        )

    def get_fc(self, amp_id, read_count):
        """
        Return how many families (frequency count, FC) have the given read_count for amp_id.

        Parameters
        ----------
        amp_id : str
        read_count : int

        Returns
        -------
        int
            Number of families for which the size == read_count.
        """
        return self._fc_dist[amp_id][read_count]

    def uid_family_size_stats(self, amp_id=None):
        """
        Compute mean and standard deviation of UID family sizes.

        Parameters
        ----------
        amp_id : str or None
            If specified, compute statistics for that single amplicon/quantispike ID.
            If None, aggregate across all amplicons/quantispikes.

        Returns
        -------
        (float, float)
            (mean_size, std_size)
        """
        if amp_id is None:
            # Flatten all readcounts across all IDs
            fam_sizes = np.array([c for counts in self._readcounts.values() for c in counts])
        else:
            fam_sizes = np.array(self._readcounts[amp_id])

        return fam_sizes.mean(), fam_sizes.std()

