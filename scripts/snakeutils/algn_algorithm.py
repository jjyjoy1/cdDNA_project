###############################################################################
# Refactored versions of: ukkonen, single_best_hit, single_best_primer_hit,
# and the Alignment class.
###############################################################################

import itertools
from collections import deque
import numpy as np
import re

###############################################################################
# Constants as defined in alignment_algo.py
###############################################################################
BAM_EQUAL = 7
BAM_DIFF = 8
BAM_INSERTION = 1
BAM_DELETION = 2

_SAM_CHARS = {
    BAM_EQUAL: '=',
    BAM_DIFF: 'X',
    BAM_INSERTION: 'I',
    BAM_DELETION: 'D'
}

_ALIGN_MATCH = 0
_ALIGN_INSERTION = 1
_ALIGN_DELETION = 2


def ukkonen(text, pattern, k):
    """
    Perform an approximate pattern matching using a Ukkonen-like algorithm,
    allowing up to k errors (substitutions, insertions, or deletions).

    Parameters
    ----------
    text : str
        Reference sequence (the 'haystack').
    pattern : str
        Query sequence to match (the 'needle').
    k : int
        Maximum permitted errors (sum of mismatches + indels).

    Returns
    -------
    list of tuples
        A list of (cost, end_position) where:
          - cost is how many edits were required (<= k),
          - end_position is the ending index of the match in 'text'.
        The matches are discovered while scanning through 'text'.

    Notes
    -----
    - If either `text` or `pattern` is empty, no matches are returned.
    - cost[j] is dynamically updated to track the edit distance
      when aligning pattern up to index j with the current character
      in text.
    - This implementation slides along `text` one character at a time and
      updates the cost array for each position in `pattern`.
    """
    if not text or not pattern:
        return []

    # cost[j] will track how many edits are needed to match pattern[:j+1]
    # against the current portion of text.
    # Initialize cost array for length = len(pattern).
    cost = list(range(1, len(pattern) + 1))

    # We'll keep track of “lact” as the rightmost position we’re still
    # allowed to check in the cost array (i.e., cost[lact] <= k).
    # This helps skip positions that exceed the error threshold.
    lact = min(k + 1, len(pattern) - 1)

    hits = []

    # Slide over each character in the text
    for i, char_t in enumerate(text):
        cost_diag = 0   # the cost for the previous diagonal cell
        cost_ij = 0

        # Update cost array from left to right
        for j, char_p in enumerate(pattern):
            # Record the old cost[j] to use as cost_diag in the next iteration
            old_cost_j = cost[j]

            # If the characters match exactly, cost_ij = cost_diag
            # (no additional edits). Otherwise, add 1 for a mismatch.
            if char_t == char_p:
                cost_ij = cost_diag
            else:
                # min(...) + 1 covers substitution, insertion, or deletion
                cost_ij = min(cost_ij, cost_diag, cost[j]) + 1

            # Move cost_ij into cost[j], then store old cost_j in cost_diag
            cost[j] = cost_ij
            cost_diag = old_cost_j

        # Try to adjust lact so we skip positions that can't produce a valid match
        while lact >= 0 and cost[lact] > k:
            lact -= 1
        if lact < len(pattern) - 1:
            lact += 1

        # If we've progressed far enough along text that a match of pattern could
        # fully “fit,” record it. The idea is that once i >= len(pattern)-1,
        # we have a potential alignment that ends around i.
        if i >= len(pattern) - 1 - k:
            # cost[lact] is a known distance, and lact points at the best feasible
            # position or is at len(pattern)-1. So cost[lact] <= k is a valid match.
            hits.append((cost[lact], i))

    return hits


def single_best_hit(hitlist):
    """
    From a list of (cost, position) hits (e.g., from ukkonen), return
    the single “best” hit. The best is defined as:
      1. the lowest cost (fewest errors)
      2. if multiple share the same cost, choose the rightmost position.

    Parameters
    ----------
    hitlist : list of (int, int)
        Each element is (cost, end_position).

    Returns
    -------
    (best_cost, best_pos) : (int, int) or (None, None)
        The best cost (fewest errors) and the rightmost end_position for that cost.
        If hitlist is empty, returns (None, None).
    """
    if not hitlist:
        return (None, None)
    # Sort primarily by cost ascending, secondarily by position descending
    # Then pick the first item in the sorted sequence
    best_cost, best_pos = sorted(hitlist, key=lambda x: (x[0], -x[1]))[0]
    return best_cost, best_pos


def single_best_primer_hit(hitlist, primer_length):
    """
    Pick the single best hit (lowest cost, then rightmost position) from
    a list of (cost, position), and adjust if needed for a known offset issue
    with certain primers (common “KRAS / NRAS117 problem”).

    Parameters
    ----------
    hitlist : list of (int, int)
        Each element is (cost, end_position).
    primer_length : int
        The length of the primer we’re aligning.

    Returns
    -------
    (final_cost, final_pos) : (int, int) or (None, None)
        The “best” cost and position. If `hitlist` is empty, returns (None, None).

    Notes
    -----
    - If the best cost is not None, we check if the best position is strictly
      > primer_length - 1 AND if the next position down is also in the hitlist
      with the same cost. If so, we shift one position to the left.
      This is a known fix for certain corner cases in primer alignment.
    """
    cost, pos = single_best_hit(hitlist)
    if cost is None:
        return (None, None)

    # Attempt to fix a known off-by-1 scenario
    # Only do it if the best position is strictly beyond primer_length - 1
    if pos is not None and pos > primer_length - 1:
        # If that shift is in the list, we shift left by 1
        # E.g., if (cost, pos-1) also in hitlist, we pick pos-1
        if (cost, pos - 1) in hitlist:
            pos = pos - 1
    return (cost, pos)


class Alignment:
    """
    Represents an alignment between two DNA sequences. Provides utilities
    to examine mismatches, insertions, and deletions, as well as to generate
    mutation sets and CIGAR-like strings.

    Parameters
    ----------
    seq_ref : str
        The “reference” sequence.
    seq_query : str
        The “query” sequence to be aligned against seq_ref.
    align : list of int, optional
        Encoded alignment operations, using e.g. BAM_EQUAL, BAM_DIFF, etc.
        If None, the alignment object is considered invalid (Falsey).
    cost : float, optional
        Numerical “cost” or “score” of the alignment as returned by the
        alignment algorithm. Interpreted based on the algorithm.

    Attributes
    ----------
    _seq_ref : str
        Internal storage for reference.
    _seq_query : str
        Internal storage for query.
    _align_arr : np.array of int or None
        Operations describing the alignment (see BAM_* constants).
    _cost : float
        The alignment’s “cost” or “score.”
    """

    def __init__(self, seq_ref, seq_query, align, cost):
        self._seq_ref = seq_ref
        self._seq_query = seq_query
        self._align_arr = np.array(align) if align is not None else None
        self._cost = cost

    def __bool__(self):
        """
        Return True if this alignment is valid (has a non-None alignment array),
        False otherwise.
        """
        return self._align_arr is not None

    def __str__(self):
        """
        Produce a 3-line representation of the alignment with the reference
        on the first line, the query on the last line, and a middle line that
        indicates equality, difference, insertion, or deletion at each position.
        """
        if self._align_arr is None:
            return "(Invalid alignment)"

        alignment_len = len(self._align_arr)
        ref_line = np.array(['-'] * alignment_len)
        query_line = np.array(['-'] * alignment_len)
        mid_line = np.array([' '] * alignment_len)

        # Fill in reference characters except where we have an insertion
        ref_line[self._align_arr != BAM_INSERTION] = list(self._seq_ref)
        # Fill in query characters except where we have a deletion
        query_line[self._align_arr != BAM_DELETION] = list(self._seq_query)

        mid_line[self._align_arr == BAM_EQUAL] = '|'
        mid_line[self._align_arr == BAM_DIFF] = '.'

        return '\n'.join([
            ''.join(ref_line),
            ''.join(mid_line),
            ''.join(query_line),
            ''
        ])

    def __len__(self):
        """Number of alignment operations."""
        if self._align_arr is None:
            return 0
        return len(self._align_arr)

    @property
    def seq_ref(self):
        """str: Return the reference sequence used in this alignment."""
        return self._seq_ref

    @property
    def seq_query(self):
        """str: Return the query sequence used in this alignment."""
        return self._seq_query

    @property
    def cost(self):
        """float: The alignment cost (or score)."""
        return self._cost

    def _cigar_tuples(self):
        """
        Describe the alignment as consecutive runs of the same operation.

        Returns
        -------
        list of (operation, length)
            For example, (BAM_EQUAL, 3), (BAM_INSERTION, 2), ...
        """
        if self._align_arr is None:
            return []
        return [
            (op, sum(1 for _ in group))
            for op, group in itertools.groupby(self._align_arr)
        ]

    def cigar_list(self):
        """
        Get a list of counted operations as used in SAM files.

        Returns
        -------
        list of str
            E.g. ["10=", "1X", "2I"] etc.
        """
        return [
            f"{count}{_SAM_CHARS[op]}"
            for (op, count) in self._cigar_tuples()
        ]

    def has_unclear_indel(self):
        """
        Check if there are multiple sets of consecutive insertions or deletions
        that are interrupted by a different alignment operation.

        Returns
        -------
        bool
            True if we detect separated runs of the same gap operation.
        """
        arr = self._align_arr
        if arr is None or len(arr) < 2:
            return False

        def _gap_is_broken(gap_type):
            # Count how many times we see gap_type appear AFTER a different op.
            # E.g. `X I I X I I` => two separate runs of 'I'
            # We'll look for transitions from non-gap to gap.
            return np.logical_and(arr[:-1] != gap_type, arr[1:] == gap_type).sum()

        # If we have more than 1 separate run of insertions or deletions => unclear
        return _gap_is_broken(BAM_INSERTION) > 1 or _gap_is_broken(BAM_DELETION) > 1

    def count_indels_and_mismatches(self):
        """
        Count all operations that are NOT a perfect match.

        Returns
        -------
        int
        """
        if self._align_arr is None:
            return 0
        return (self._align_arr != BAM_EQUAL).sum()

    def count_indels_and_mismatches_borders(self, margin):
        """
        Count all operations that are NOT a perfect match, but only in the first
        or last `margin` operations of the alignment.

        Parameters
        ----------
        margin : int
            Number of alignment operations at each end to count.

        Returns
        -------
        int
            Sum of mismatches/indels in the first `margin` and last `margin` positions.
        """
        if self._align_arr is None:
            return 0
        arr = self._align_arr
        return (arr[:margin] != BAM_EQUAL).sum() + (arr[-margin:] != BAM_EQUAL).sum()

    def count_insertions_end(self, margin):
        """
        Count how many insertion operations occur in the last `margin` alignment ops.

        Parameters
        ----------
        margin : int

        Returns
        -------
        int
        """
        if self._align_arr is None:
            return 0
        return (self._align_arr[-margin:] == BAM_INSERTION).sum()

    def count_deletions_end(self, margin):
        """
        Count how many deletion operations occur in the last `margin` alignment ops.

        Parameters
        ----------
        margin : int

        Returns
        -------
        int
        """
        if self._align_arr is None:
            return 0
        return (self._align_arr[-margin:] == BAM_DELETION).sum()

    def indel_mismatch_ratio(self):
        """
        Fraction of alignment operations that are not perfect matches.

        Returns
        -------
        float
        """
        total_ops = len(self)
        if total_ops == 0:
            return 0.0
        return self.count_indels_and_mismatches() / total_ops

    def count_matches(self):
        """Number of positions that are perfect matches."""
        if self._align_arr is None:
            return 0
        return (self._align_arr == BAM_EQUAL).sum()

    def count_insertions(self):
        """Number of insertion operations in the alignment."""
        if self._align_arr is None:
            return 0
        return (self._align_arr == BAM_INSERTION).sum()

    def count_deletions(self):
        """Number of deletion operations in the alignment."""
        if self._align_arr is None:
            return 0
        return (self._align_arr == BAM_DELETION).sum()

    def count_mismatches(self):
        """Number of mismatch operations in the alignment."""
        if self._align_arr is None:
            return 0
        return (self._align_arr == BAM_DIFF).sum()

    def get_mutations(self, shift_ref=0, shift_query=0):
        """
        Convert alignment operations into a set of “mutation” records
        (inserts, deletions, substitutions). Large consecutive indels
        become a single mutation event.

        Parameters
        ----------
        shift_ref : int
            Offset applied to reference positions (useful if seq_ref is
            a subregion of a bigger reference).
        shift_query : int
            Offset for the query positions.

        Returns
        -------
        MutationSet
            The collection of observed mutation events.
        """
        # We'll re-use the existing logic from alignment_algo.py:
        from collections import namedtuple

        # Minimal local version of the MutationSet / Mutation structure:
        Mutation = namedtuple("Mutation", ["pos_ref", "pos_query", "seq_ref", "seq_query"])

        class MutationSet:
            def __init__(self):
                self._muts = []

            def add(self, pos_r, pos_q, s_ref, s_query):
                self._muts.append(Mutation(pos_r, pos_q, s_ref, s_query))

            def __repr__(self):
                return f"MutationSet({self._muts})"

        mutations = MutationSet()

        if self._align_arr is None:
            return mutations

        i = 0
        j = 0
        # group by consecutive runs of the same operation
        def cigar_tuples():
            return [
                (op, sum(1 for _ in g)) for op, g in itertools.groupby(self._align_arr)
            ]

        for (op, count) in cigar_tuples():
            if op == BAM_DELETION:
                # (ref segment, empty query)
                ref_str = self._seq_ref[i : i + count]
                mutations.add(
                    i + shift_ref, j + shift_query, ref_str, ""
                )
                i += count
            elif op == BAM_INSERTION:
                # (empty ref, query segment)
                query_str = self._seq_query[j : j + count]
                mutations.add(
                    i + shift_ref, j + shift_query, "", query_str
                )
                j += count
            elif op == BAM_DIFF:
                # break each mismatch into a single-base event
                for _ in range(count):
                    ref_char = self._seq_ref[i]
                    query_char = self._seq_query[j]
                    mutations.add(
                        i + shift_ref, j + shift_query, ref_char, query_char
                    )
                    i += 1
                    j += 1
            else:
                # op == BAM_EQUAL => skip
                i += count
                j += count

        return mutations

    def shift_end_gaps(self):
        """
        Attempt to shift trailing gaps if there’s a more “natural” alignment
        that avoids internal mismatches. If cost is the same,
        the shift can produce a more meaningful alignment.

        Returns
        -------
        Alignment
            Possibly adjusted alignment, or self if no shift was done.
        """
        # Implementation detail is the same as the original code's logic:
        # We look for the last block of gaps, see if we can shift them earlier.
        # If a shift would keep cost the same but align more identical bases,
        # do so. Otherwise, return unmodified.
        if self._align_arr is None or len(self._align_arr) == 0:
            return self

        arr_rev = self._align_arr[::-1]
        # i: the index in reverse array for the first match at the end
        i = np.argmax(arr_rev == BAM_EQUAL)
        # j: index in reverse array for the first gap after i
        j = i + np.argmax(arr_rev[i:] != BAM_EQUAL)

        # If j is beyond array length, no gap found
        if j >= len(arr_rev):
            return self

        indel_op = arr_rev[j]
        if indel_op == BAM_DIFF:
            # The original code just returns self if the op is mismatch
            return self

        # If all the preceding ops up to i are not the same gap, return self
        if not (arr_rev[:i] == indel_op).all():
            return self

        k = j + np.argmax(arr_rev[j:] != indel_op)
        gap_length = k - j
        match_length = j - i

        # If the match region is smaller than the gap region, shifting won't help
        if match_length <= gap_length:
            return self

        # Try to see if we can realign the last match region
        # so the cost remains the same but the alignment is more “natural.”
        alignment_len = len(self._align_arr)
        pos_rev = alignment_len - k
        # Number of insert ops in the region after pos_rev
        ins_count_after = (self._align_arr[k:] == BAM_INSERTION).sum()
        # Number of del ops in the region after pos_rev
        del_count_after = (self._align_arr[k:] == BAM_DELETION).sum()

        ir = pos_rev - ins_count_after
        iq = pos_rev - del_count_after

        # Check if the reference and query substrings match, so shifting is beneficial
        ref_sub = self._seq_ref[ir : ir + match_length]
        query_sub = self._seq_query[iq : iq + match_length]
        if ref_sub == query_sub:
            # Build a new alignment array with the shift
            new_align = np.concatenate([
                self._align_arr[:pos_rev],
                np.array([BAM_EQUAL] * match_length),
                np.array([indel_op] * (alignment_len - pos_rev - match_length))
            ])
            return Alignment(self._seq_ref, self._seq_query, new_align, self._cost)
        return self

    @classmethod
    def from_str(cls, str_ref, str_query):
        """
        Construct an Alignment object from two “aligned” strings:
        reference and query, each possibly containing '-' characters
        to denote gaps.

        Parameters
        ----------
        str_ref : str
            The reference alignment string (with '-').
        str_query : str
            The query alignment string (with '-').

        Returns
        -------
        Alignment
        """
        # They should be the same length (including '-')
        if len(str_ref) != len(str_query):
            raise ValueError("Aligned strings must be the same length")

        arr_len = len(str_ref)
        align_ops = np.empty(arr_len, dtype=int)

        # We'll build the real (gapless) seq_ref, seq_query as we go
        real_ref = []
        real_query = []

        for i in range(arr_len):
            r_char, q_char = str_ref[i], str_query[i]
            if r_char == '-':
                align_ops[i] = BAM_INSERTION
                real_query.append(q_char)
            elif q_char == '-':
                align_ops[i] = BAM_DELETION
                real_ref.append(r_char)
            elif r_char == q_char:
                align_ops[i] = BAM_EQUAL
                real_ref.append(r_char)
                real_query.append(q_char)
            else:
                align_ops[i] = BAM_DIFF
                real_ref.append(r_char)
                real_query.append(q_char)

        cost_val = np.count_nonzero(align_ops != BAM_EQUAL)
        return cls(
            seq_ref=''.join(real_ref),
            seq_query=''.join(real_query),
            align=align_ops,
            cost=cost_val
        )

    # Class method that delegates to a generic alignment function
    @classmethod
    def _align(cls, seq_ref, seq_query, align_func, **params):
        """
        Generic helper to perform an alignment using the given function
        and build an Alignment object.

        Parameters
        ----------
        seq_ref : str
        seq_query : str
        align_func : callable
            Something that returns (cost, array_of_ops)
        **params : dict
            Additional parameters for align_func

        Returns
        -------
        Alignment
        """
        cost_val, ops = align_func(seq_ref, seq_query, **params)
        return cls(seq_ref, seq_query, ops, cost_val)



