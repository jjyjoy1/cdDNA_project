import os
from multiprocessing import Pool

from pss.utils import iter_file_nocomment
from pss.seqtools import convert_sense
from pss.alignment_algo import (
    ukkonen,
    single_best_primer_hit,
    single_best_hit,
    Alignment
)
from pss.amplicons import AmpliconGroup

import collections
from collections import defaultdict

# If your code references Amplicon, AmpliconGroup, or other classes, import them:
# from pss.amplicons import AmpliconGroup


###############################################################################
# Helper Classes for Primer Matching
###############################################################################
class PrimerMatch:
    """
    Stores information about how a single primer matched a read.

    Parameters
    ----------
    amplicon_id : str
        Identifier of the amplicon whose primer matched.
    reverse_oriented : bool
        True if this match corresponds to the reverse primer.
    cost : int
        The alignment or mismatch cost for the primer match.
    start_position : int
        Start position (0-based) of the matched region on the read.
    end_position : int
        End position (0-based, exclusive) of the matched region on the read.

    Attributes
    ----------
    amplicon_id : str
    reverse_oriented : bool
    cost : int
    start_position : int
    end_position : int
    """

    def __init__(self,
                 amplicon_id,
                 reverse_oriented,
                 cost,
                 start_position,
                 end_position):
        self.amplicon_id = amplicon_id
        self.reverse_oriented = reverse_oriented
        self.cost = cost
        self.start_position = start_position
        self.end_position = end_position

    @property
    def is_forward(self):
        """bool: True if this is a forward primer."""
        return not self.reverse_oriented

    @property
    def matchlen(self):
        """int: Length of the matched primer sequence."""
        return self.end_position - self.start_position


class PrimerMatchPair:
    """
    Encapsulates two primer matches (forward and reverse),
    or partial matches if only one primer was found.

    Attributes
    ----------
    first_primer : PrimerMatch or None
    second_primer : PrimerMatch or None
    template_id : str or None
        Used if the read matched a 'template' region instead of standard primers.
    misprime : bool
        True if the read matched two primers from different amplicons (i.e., mismatch).
    """

    def __init__(self):
        self.first_primer = None
        self.second_primer = None
        self.template_id = None
        self.misprime = False

    def __str__(self):
        return "\n".join([
            "Primer Match",
            "status = " + self.status(),
            "amplicon(s) = " + self.get_id(),
            ""
        ])

    @property
    def is_forward(self):
        """bool: True if the 'first_primer' is forward-oriented."""
        return self.first_primer and self.first_primer.is_forward

    @property
    def is_reverse(self):
        """bool: True if the 'first_primer' is reverse-oriented."""
        return self.first_primer and self.first_primer.reverse_oriented

    @property
    def only_one_primer_found(self):
        """bool: True if exactly one primer was matched."""
        return (
            self.first_primer is not None
            and self.second_primer is None
        )

    @property
    def no_primers_found(self):
        """bool: True if no primers were matched."""
        return self.first_primer is None

    @property
    def matching_primers_found(self):
        """bool: True if both primers from the same amplicon were found."""
        return (
            self.first_primer is not None
            and self.second_primer is not None
            and not self.misprime
        )

    def status(self):
        """
        Returns a string describing the pairing status:
         - 'ok' if both primers from the same amplicon are found
         - 'not found' if no primer is found
         - 'template' if a template region (rather than a standard primer) was matched
         - 'misprime' if two different amplicon primers were found
        """
        if self.matching_primers_found:
            return 'ok'
        elif self.first_primer is None:
            return 'not found'
        elif self.template_id:
            return 'template'
        return 'misprime'

    def get_id(self):
        """
        Returns a string describing the amplicon ID(s) or template.

        Returns
        -------
        str
        """
        if self.template_id:
            return self.template_id
        if self.matching_primers_found:
            return self.first_primer.amplicon_id
        if self.first_primer is None:
            return "None"
        if self.second_primer is None:
            return self.first_primer.amplicon_id
        return "+".join([
            self.first_primer.amplicon_id,
            self.second_primer.amplicon_id
        ])

    def read_align_region(self, read, margin=0):
        """
        Extract the portion of the read spanned by first_primer -> second_primer,
        optionally including an extra margin of bases on each side.

        Parameters
        ----------
        read : str
            The full read sequence.
        margin : int
            Number of extra bases to include on each side.

        Returns
        -------
        str
            Substring of read from (first_primer.end_position - margin)
            to (second_primer.start_position + margin).
        """
        if not self.first_primer or not self.second_primer:
            return ""
        start = max(0, self.first_primer.end_position - margin)
        end = min(len(read), self.second_primer.start_position + margin)
        return read[start:end]


def _get_best_primer_match(primer_matches):
    """
    From a list of PrimerMatch objects, pick the best (lowest cost,
    then earliest start).

    Parameters
    ----------
    primer_matches : list of PrimerMatch

    Returns
    -------
    PrimerMatch
    """
    if not primer_matches:
        return None
    # sort by (cost ASC, start_position ASC) and pick the first
    primer_matches_sorted = sorted(primer_matches, key=lambda pm: (pm.cost, pm.start_position))
    return primer_matches_sorted[0]


###############################################################################
# PrimerFinder
###############################################################################
class PrimerFinder:
    """
    Locates forward/reverse primers in a read by approximate string matching
    (using ukkonen and single_best_primer_hit).

    Parameters
    ----------
    amplicons : AmpliconGroup or similar
        Container that provides a list of amplicons, each with primer sequences.

    Attributes
    ----------
    _amplicons : AmpliconGroup
    _primer_max_err : int
        Allowed error for the quick ukkonen-based primer detection.
    _primer_match_len : int
        Subset of the primer (suffix or prefix) used for alignment in _prim_align.
    """

    def __init__(self, amplicons):
        self._amplicons = amplicons
        self._primer_max_err = 1
        self._primer_match_len = 16

    def _match_prim(self, read, primer):
        """
        Use ukkonen + single_best_primer_hit to find the best match of `primer`
        in the initial region of `read`. Return (cost, start, end) or (None, None, None).

        We only search read[:detect_limit] to reduce computational overhead.
        """
        detect_limit = max(len(primer), len(read) // 2)
        hits = ukkonen(read[:detect_limit], primer, self._primer_max_err)
        cost, end = single_best_primer_hit(hits, len(primer))
        if cost is None:
            return None, None, None
        start = end - len(primer) + 1
        return cost, start, end

    def _prim_align(self, read, primer):
        """
        Align the last `_primer_match_len` bases of a primer to the read,
        then infer the full primer's position.

        Returns (cost, start, end) or (None, None, None).
        """
        primer_suffix = primer[-self._primer_match_len:]
        cost, start, end = self._match_prim(read, primer_suffix)
        if start is not None:
            # shift the start to account for the rest of the primer
            offset = len(primer) - len(primer_suffix)
            start = max(0, start - offset)
        return cost, start, end

    def _prim_align_prefix(self, read, primer):
        """
        Align the first `_primer_match_len` bases of a primer to the read,
        then infer the full primer's position.

        Returns (cost, start, end) or (None, None, None).
        """
        primer_prefix = primer[:self._primer_match_len]
        cost, start, end = self._match_prim(read, primer_prefix)
        if end is not None:
            # shift the end to account for the rest of the primer
            offset = len(primer) - len(primer_prefix)
            end += offset
        return cost, start, end

    def _primer_matches_amplicon(self, amp, read):
        """
        Attempt to match both primers (forward + reverse) of `amp` to `read`.
        Returns a list of PrimerMatch instances.

        Parameters
        ----------
        amp : Amplicon
            Holds primer sequences (both_primers).
        read : str

        Returns
        -------
        list of PrimerMatch
        """
        results = []
        for i, primer in enumerate(amp.both_primers):
            cost, start, end = self._prim_align(read, primer)
            if cost is not None:
                # i==0 => forward, i==1 => reverse
                results.append(
                    PrimerMatch(amp.id, bool(i), cost, start, end+1)
                )
        return results

    def detect_primer(self, read):
        """
        Attempt to find a matching forward primer and reverse primer
        for the read across all amplicons. If found, returns a PrimerMatchPair
        with both matches. If not, returns partial or empty matches.

        Parameters
        ----------
        read : str

        Returns
        -------
        PrimerMatchPair
        """
        match_pair = PrimerMatchPair()

        # 1) Find all forward-primer matches
        forward_matches = []
        for amp in self._amplicons:
            fmatches = self._primer_matches_amplicon(amp, read)
            forward_matches.extend(fmatches)

        if not forward_matches:
            return match_pair  # no primers found at all

        # 2) Pick the single best forward primer
        best_fp = _get_best_primer_match(forward_matches)
        match_pair.first_primer = best_fp

        # 3) Identify the associated amplicon for that primer
        amp_fp = self._amplicons[best_fp.amplicon_id]
        reverse_primer = amp_fp.get_primer(forward=best_fp.reverse_oriented)

        # 4) Attempt to match the reverse primer on the reverse complement of the read
        read_rev = convert_sense(read)
        for align_func in [self._prim_align, self._prim_align_prefix]:
            cost, start, end = align_func(read_rev, reverse_primer)
            if cost is not None:
                # success => store second primer, done
                match_pair.second_primer = PrimerMatch(
                    amp_fp.id,
                    best_fp.reverse_oriented,
                    cost,
                    len(read_rev) - end - 1,
                    len(read_rev) - start
                )
                return match_pair

        # 5) Check if this read might align to a "template" region
        #    (some custom or partial amplicon region)
        for template_id, seq, err in amp_fp.iter_templates():
            if ukkonen(read, seq, err):
                match_pair.template_id = template_id
                return match_pair

        # 6) Possibly a misprime event (search for other amplicons' reverse primer in read_rev)
        revp_matches = []
        for amp in self._amplicons:
            rmatches = self._primer_matches_amplicon(amp, read_rev)
            revp_matches.extend(rmatches)

        if revp_matches:
            match_pair.second_primer = _get_best_primer_match(revp_matches)
            match_pair.misprime = True

        return match_pair


###############################################################################
# AlignedRead
###############################################################################
class AlignedRead:
    """
    Describes the result of matching a read to a pair of primers and (if found)
    performing an alignment to the corresponding amplicon sequence.

    Parameters
    ----------
    primer_match_pair : PrimerMatchPair or None
    alignment : Alignment or None
    amplicon : Amplicon or None
    failed : bool
        True if the read was attempted but alignment could not be completed.
    align_margin : int
        Number of bases used as margin in the alignment region.

    Attributes
    ----------
    alignment_failed : bool
        True if alignment was attempted but failed.
    is_pseudogene : bool
        Post-processing flag set if the read is identified as pseudogene.
    quantispike : str or None
        ID of the quantispike if matched, or 'unclear_...' if partial, else None.
    """

    def __init__(self,
                 primer_match_pair=None,
                 alignment=None,
                 amplicon=None,
                 failed=False,
                 align_margin=0):
        self._primer_match_pair = primer_match_pair
        self._alignment = alignment
        self._amplicon = amplicon
        self._mutations = None
        self._mutations_not_filtered = None
        self.alignment_failed = failed
        self._align_margin = align_margin
        self.is_pseudogene = False
        self.quantispike = None
        self._complex_mutations = None
        self._max_mut_dist = -1

        # If we do have an alignment, extract its mutation sets
        if alignment and amplicon and primer_match_pair:
            # The reference shift is the length of the matched primer
            # so that alignment mutation positions line up with the full amplicon
            forward_flag = primer_match_pair.is_forward
            start_ref = len(amplicon.get_primer(forward_flag))
            start_query = primer_match_pair.first_primer.end_position

            # Mutations including primer region
            self._mutations_not_filtered = alignment.get_mutations(
                shift_ref=(start_ref - align_margin),
                shift_query=(start_query - align_margin)
            )
            # Filter out mutations in the primer region
            end_ref = start_ref + len(amplicon.get_ampseq())
            self._mutations = self._mutations_not_filtered.filter_region(
                start_ref, end_ref
            )

    def __str__(self):
        lines = [
            str(self._primer_match_pair) if self._primer_match_pair else "(No PrimerMatchPair)",
            str(self._alignment) if self._alignment else "(No Alignment)",
            f"Amplicon is forward: {self._amplicon.is_forward if self._amplicon else 'N/A'}",
            f"Matched in Amplicon direction: {self.matched_forward}",
            f"Pseudogene: {self.is_pseudogene}",
            f"Quantispike: {self.quantispike}"
        ]
        return "\n".join(lines) + "\n"

    @property
    def primer_match_pair(self):
        """PrimerMatchPair: The matched primer info."""
        return self._primer_match_pair

    @property
    def alignment(self):
        """Alignment: The alignment object, if any."""
        return self._alignment

    @property
    def amplicon(self):
        """Amplicon: The amplicon object, if any."""
        return self._amplicon

    @property
    def chromosome(self):
        """str: Chromosome name, if amplicon is known."""
        return self._amplicon.chromosome if self._amplicon else "unknown"

    @property
    def no_primers_found(self):
        """bool: True if no primers matched."""
        return self._primer_match_pair.no_primers_found if self._primer_match_pair else True

    @property
    def only_one_primer_found(self):
        """bool: True if only one primer was found."""
        return self._primer_match_pair.only_one_primer_found if self._primer_match_pair else False

    @property
    def misprime(self):
        """bool: True if two primers from different amplicons were matched."""
        return self._primer_match_pair.misprime if self._primer_match_pair else False

    @property
    def matching_primers_found(self):
        """bool: True if both primers for the same amplicon were found."""
        return self._primer_match_pair.matching_primers_found if self._primer_match_pair else False

    def is_amplicon_read(self):
        """
        Determine if the read is a valid amplicon read, i.e.:
          - matching primers found
          - alignment didn't fail
          - not flagged as pseudogene
          - not flagged as quantispike

        Returns
        -------
        bool
        """
        return (
            self.matching_primers_found
            and not self.alignment_failed
            and not self.is_pseudogene
            and not self.quantispike
        )

    @property
    def matched_forward(self):
        """bool: True if the read alignment is in the same orientation as the amplicon."""
        if not self._primer_match_pair:
            return False
        return self._primer_match_pair.is_forward

    @property
    def on_complementary_strand(self):
        """
        bool: True if the read alignment is on the opposite strand
        of the reference (amplicon).

        If the amplicon is forward, but the read matched in reverse orientation,
        or vice versa, return True.
        """
        if not self._amplicon or not self._primer_match_pair:
            return False
        return (self._amplicon.is_forward != self._primer_match_pair.is_forward)

    @property
    def query_seq(self):
        """str: The aligned portion of the read, if alignment is available."""
        return self._alignment.seq_query if self._alignment else ""

    @property
    def amplicon_id(self):
        """str: The amplicon identifier from the matched primer(s)."""
        return self._primer_match_pair.get_id() if self._primer_match_pair else "None"

    @property
    def amp_qs_id(self):
        """
        str: If this read is determined to be a quantispike,
        return its ID. If it's an 'unclear' quantispike, return 'None'.
        Otherwise, return the amplicon_id.
        """
        if self.is_unclear_quantispike:
            return "None"
        if self.is_quantispike:
            return self.quantispike
        return self.amplicon_id

    @property
    def mutations(self):
        """
        The subset of alignment mutations outside the primer regions.

        Returns
        -------
        MutationSet or None
        """
        return self._mutations

    @property
    def mutations_not_filtered(self):
        """
        The full set of alignment mutations, including those in primer regions.

        Returns
        -------
        MutationSet or None
        """
        return self._mutations_not_filtered

    @property
    def amplicon_key(self):
        """
        Returns (amplicon_id, matched_forward) for dictionary keys,
        or ("None", False) if no amplicon was matched.

        Returns
        -------
        tuple of (str, bool)
        """
        return (self.amplicon_id, self.matched_forward)

    @property
    def is_unclear_quantispike(self):
        """
        bool: True if the read was detected as a quantispike,
        but labeled 'unclear_...'.
        """
        return self.quantispike and "unclear_" in self.quantispike

    @property
    def is_quantispike(self):
        """
        bool: True if the read was detected as a definite quantispike.
        """
        return self.quantispike and not self.is_unclear_quantispike

    @property
    def primer_ids(self):
        """
        tuple of str: (forward_primer_amplicon_id, reverse_primer_amplicon_id)
        """
        if not self._primer_match_pair or not self._primer_match_pair.second_primer:
            return (None, None)
        return (
            self._primer_match_pair.first_primer.amplicon_id,
            self._primer_match_pair.second_primer.amplicon_id
        )

    def get_quality_substring(self, qual):
        """
        Extract the portion of the quality string that corresponds
        to the aligned portion of the read.

        Reverses it if the primer match is on the reverse orientation.

        Parameters
        ----------
        qual : str
            Full phred-quality string for the read.

        Returns
        -------
        str
        """
        if not self._primer_match_pair or not self._alignment:
            return ""

        start = self._primer_match_pair.first_primer.end_position - self._align_margin
        end = start + len(self._alignment.seq_query)
        qual_sub = qual[start:end]
        if self._primer_match_pair.is_reverse:
            return qual_sub[::-1]
        return qual_sub

    def get_complex_mutations(self, max_mut_dist):
        """
        Group the read's alignment mutations into 'complex' clusters, if applicable.

        Parameters
        ----------
        max_mut_dist : int
            Maximum distance between adjacent mutations to be merged into one group.

        Returns
        -------
        list of MutationSet
        """
        if not self._mutations:
            return []
        if self._complex_mutations and self._max_mut_dist == max_mut_dist:
            return self._complex_mutations
        self._max_mut_dist = max_mut_dist
        self._complex_mutations = self._mutations.group_complex_mutations(max_mut_dist)
        return self._complex_mutations

    def alignment_with_mutations(self):
        """
        Check if the alignment has any indels or mismatches.

        Returns
        -------
        bool
        """
        if not self._alignment:
            return False
        return (self._alignment.count_indels_and_mismatches() > 0)



class ReadAligner:
    """
    Aligns reads against a set of amplicons by first matching primers
    and then performing a sequence alignment with flexible parameters.

    Parameters
    ----------
    amplicons : AmpliconGroup (or similar)
        Container with the reference amplicon sequences and associated data.

    Attributes
    ----------
    _primer_finder : PrimerFinder
        Helper object that detects which primer(s) match a given read.
    _amplicons : AmpliconGroup
        The main amplicon group for the assay.
    _left_margin : int
        How many bases to include upstream of the read region for alignment.
    _total_margin : int
        Overall margin used in alignment truncations.
    _gotoh_par_one : dict
        Parameters for one version of the Gotoh k-band alignment.
    _gotoh_par_match_enhanced : dict
        Enhanced mismatch penalty parameters for tricky alignments.
    _gotoh_par_low_gapext : dict
        Parameters that reduce gap extension cost for especially tricky reads.
    """

    def __init__(self, amplicons):
        from .align_read import PrimerFinder  # or define/import PrimerFinder as needed

        self._primer_finder = PrimerFinder(amplicons)
        self._amplicons = amplicons

        # Parameters that adjust how alignment is performed
        self._left_margin = 2
        self._total_margin = 6

        # Different parameter sets for the (k-band) Gotoh alignments
        self._gotoh_par_one = {
            'gap_start': 4,
            'gap_ext': 1,
            'mismatch': 1,
            'match_score': -1
        }
        self._gotoh_par_match_enhanced = {
            'gap_start': 4,
            'gap_ext': 1,
            'mismatch': 1,
            'match_score': -2
        }
        self._gotoh_par_low_gapext = {
            'gap_start': 4,
            'gap_ext': 0.2,
            'mismatch': 1,
            'match_score': -2
        }

    def align_read(self, read):
        """
        Align a single read to the set of amplicons. The pipeline is:
          1. Identify matching primers.
          2. If both primers are found, align the read to the selected amplicon.
          3. Apply k-band or Gotoh alignment with different parameter sets
             if needed (e.g., if too many mismatches or unclear gaps).

        Parameters
        ----------
        read : str
            The DNA read sequence.

        Returns
        -------
        AlignedRead
            An object describing how/if the read was aligned to an amplicon.
            Contains primer matching info, alignment operations, etc.
        """
        # 1) Match primers
        primer_matches = self._primer_finder.detect_primer(read)
        if not primer_matches.matching_primers_found:
            from .align_read import AlignedRead  # or define/import AlignedRead as needed
            return AlignedRead(primer_match_pair=primer_matches)

        # 2) If matched, pick the best primer info
        fp_match = primer_matches.first_primer
        amp = self._amplicons[fp_match.amplicon_id]

        # Make sure we are in the correct orientation
        forward = fp_match.is_forward
        primer_len = len(amp.get_primer(forward))
        amp_construct = amp.get_construct(forward)
        read_start = fp_match.end_position

        # 3) Attempt a quick approximate alignment using ukkonen to see
        #    how many mismatches might exist
        ukk_ref = amp_construct[primer_len:]
        ukk_query = read[read_start : read_start + len(ukk_ref)]
        cost, pos = single_best_hit(
            ukkonen(ukk_ref, ukk_query, amp.error)
        )

        from .align_read import AlignedRead  # or define/import AlignedRead as needed
        if cost is None:
            return AlignedRead(
                primer_match_pair=primer_matches,
                failed=True
            )

        # k-band width based on the approximate cost
        kband_k = max(30, 3 * cost)

        # Prepare for a more formal alignment
        # The reference for alignment includes a margin
        len_inseq = len(amp.get_ampseq())
        seq_amp = amp_construct[primer_len - self._left_margin:]
        seq_read = read[read_start - self._left_margin:]
        align_len = min(
            len_inseq + self._total_margin,
            len(seq_amp), len(seq_read)
        )

        seq_amp_trunc = seq_amp[:align_len]
        seq_read_trunc = seq_read[:align_len]

        # 4) Perform an initial k-band alignment
        align = Alignment.kband(seq_amp_trunc, seq_read_trunc, kband_k)

        # If alignment is suspicious (unclear indels, high mismatch ratio),
        # then try gotoh-based approach
        if align and (align.has_unclear_indel()
                      or align.indel_mismatch_ratio() >= 0.1
                      or align.count_indels_and_mismatches() > 4):
            align = Alignment.gotoh_kband(
                seq_amp_trunc,
                seq_read_trunc,
                kband_k,
                **self._gotoh_par_one
            )

        if align and align.indel_mismatch_ratio() >= 0.1:
            # Try a narrower alignment ignoring some margin
            seq_read_amp_len = seq_read[: len(seq_amp)]
            align = Alignment.gotoh_kband(
                seq_amp, seq_read_amp_len,
                kband_k,
                **self._gotoh_par_match_enhanced
            )
            if align and align.indel_mismatch_ratio() >= 0.05:
                align = Alignment.gotoh_kband(
                    seq_amp, seq_read_amp_len,
                    kband_k,
                    **self._gotoh_par_low_gapext
                )
        elif align and align.count_indels_and_mismatches_borders(10) >= 4:
            # Try to align ignoring some trailing parts
            align_len_read = min(len(seq_amp), len(seq_read))
            seq_amp_rl = seq_amp[:align_len_read]
            seq_read_rl = seq_read[:align_len_read]
            align = Alignment.gotoh_kband(
                seq_amp_rl, seq_read_rl,
                kband_k,
                **self._gotoh_par_match_enhanced
            )

        if align and (
            align.count_indels_and_mismatches_borders(10) >= 4
            or (
                align.count_insertions_end(10)
                == align.count_deletions_end(10) > 0
            )
        ):
            # Possibly align the full read + full reference with gotoh
            if abs(len(seq_amp) - len(seq_read)) > kband_k:
                return AlignedRead(
                    primer_match_pair=primer_matches,
                    failed=True
                )
            align = Alignment.gotoh_kband(
                seq_amp,
                seq_read,
                kband_k,
                **self._gotoh_par_one
            )

        if not align:
            return AlignedRead(primer_match_pair=primer_matches, failed=True)

        # Tweak the alignment to shift end gaps if there's a simpler layout
        align = align.shift_end_gaps()

        # Filter out extreme mismatch or large-gap alignments
        n_mis = align.count_mismatches()
        if (n_mis + align.count_insertions() > len(read) - 30
                or n_mis + align.count_deletions() > len(read) - 30):
            return AlignedRead(primer_match_pair=primer_matches, failed=True)

        # If we survived all checks, create the final AlignedRead object
        return AlignedRead(
            primer_match_pair=primer_matches,
            alignment=align,
            amplicon=amp,
            align_margin=self._left_margin
        )


##############################################################################
# 1) PseudogeneFilter
##############################################################################

def _read_pseudogene_file(fname_pseudogenes):
    """
    Parse a text/csv file containing pseudogene definitions (ID, sequence).

    Parameters
    ----------
    fname_pseudogenes : str
        Path to a file, each line containing "ID,sequence".

    Returns
    -------
    list of (str, str)
        A list of (pseudogene_id, pseudogene_sequence), all uppercase.
    """
    from pss.utils import iter_file_nocomment  # local import for clarity
    with open(fname_pseudogenes, 'r') as f:
        records = []
        for line in iter_file_nocomment(f):
            # Each line has form: PSEUDOGENE_ID,ATGCAT...
            parts = line.strip().split(",")
            if len(parts) == 2:
                pg_id, pg_seq = parts
                records.append((pg_id.upper(), pg_seq.upper()))
        return records


class PseudogeneFilter:
    """
    Identifies pseudogene reads by comparing the read's alignment mutations
    against known pseudogene mutation signatures.

    Parameters
    ----------
    fname_pseudogenes : str
        Path to a file with pseudogene definitions ("ID,sequence").
    amplicons : AmpliconGroup (or AmpliconSet)
        Container for reference amplicon sequences.

    Attributes
    ----------
    _pseudogenes : list of (str, str)
        List of (pseudogene_id, pseudogene_sequence).
    _amplicons : AmpliconGroup
        The amplicons used in alignment.
    _mutations : dict
        Maps (amplicon_id, orientation) -> list of MutationSet objects
        representing pseudogene signature(s).
    """

    def __init__(self, fname_pseudogenes, amplicons):
        self._pseudogenes = _read_pseudogene_file(fname_pseudogenes)
        self._amplicons = amplicons
        self._mutations = self._build_filter()

    def __str__(self):
        """
        Summarize the stored pseudogene mutation sets for debugging.
        """
        results = []
        for (amp_id, forward_flag), mutation_sets in self._mutations.items():
            results.append(f"Amplicon {amp_id} (is_forward={forward_flag}):")
            for mset in mutation_sets:
                results.append(str(mset) + "\n")
        return "\n".join(results)

    def _build_filter_one_pseudogene(self, pg_id, pg_seq):
        """
        For a single pseudogene, align it (and its reverse complement)
        against all amplicons matching pg_id's prefix. Collect the
        resulting mutation sets.

        Parameters
        ----------
        pg_id : str
            Pseudogene identifier.
        pg_seq : str
            Pseudogene sequence (DNA).

        Returns
        -------
        dict
            Key: (amplicon_id, orientation)
            Value: list of MutationSet objects (the pseudogene "signature").
        """
        # Subset the amplicons to those that might correspond to this pseudogene
        # (based on some naming convention, e.g. matching prefix).
        sub_amplicons = self._amplicons.get_amplicons_filtered(pg_id.split("_")[0])
        result = defaultdict(list)

        if not sub_amplicons:
            return result

        from pss.alignment_algo import Alignment
        from .align_read import ReadAligner  # relative import

        read_aligner = ReadAligner(sub_amplicons)

        # Attempt alignment in forward orientation and reversed orientation
        for seq_candidate in (pg_seq, convert_sense(pg_seq)):
            aligned = read_aligner.align_read(seq_candidate)
            if aligned.matching_primers_found:
                # If alignment succeeded, gather all "not filtered" mutations
                mut_set = aligned.mutations_not_filtered
                if len(mut_set) > 0:
                    result[aligned.amplicon_key].append(mut_set)
        return result

    def _build_filter(self):
        """
        Build the pseudogene filter across all known pseudogenes.

        Returns
        -------
        dict
            Key: (amplicon_id, orientation)
            Value: list of MutationSet objects. Each list is the set of
            known pseudogene "signatures" for that amplicon orientation.
        """
        all_results = defaultdict(list)
        for (pg_id, pg_seq) in self._pseudogenes:
            single_pg_dict = self._build_filter_one_pseudogene(pg_id, pg_seq)
            for key, mut_list in single_pg_dict.items():
                all_results[key].extend(mut_list)
        return all_results

    def detect(self, aligned_read):
        """
        Check whether a given aligned read matches any known pseudogene pattern.

        Parameters
        ----------
        aligned_read : AlignedRead
            A read aligned to an amplicon.

        Returns
        -------
        bool
            True if the read is identified as a pseudogene read, else False.
        """
        # Possibly find relevant sets of signature mutations
        candidate_mut_sets = self._mutations[aligned_read.amplicon_key]

        # Compare the read's mutations against each pseudogene signature
        read_muts = aligned_read.mutations_not_filtered

        for pg_mut_set in candidate_mut_sets:
            # Number of mutations in common
            shared_count = pg_mut_set.count_same(read_muts)
            n_ref = len(pg_mut_set)

            # The checks below replicate the original logic for deciding
            # a match:
            if shared_count == n_ref:
                return True
            if n_ref > 5 and shared_count >= 4:
                return True
            if n_ref > 2 and shared_count >= n_ref - 1:
                return True
        return False

    def get_all_sequences(self, add_reverse=True):
        """
        Return all known pseudogene DNA sequences.

        Parameters
        ----------
        add_reverse : bool
            If True, also include each pseudogene's reverse complement.

        Returns
        -------
        list of str
        """
        seqs = [seq for (_id, seq) in self._pseudogenes]
        if add_reverse:
            seqs += [convert_sense(s) for s in seqs]
        return seqs


##############################################################################
# 2) QuantispikeFilter
##############################################################################

class QuantispikeFilter:
    """
    Detect “quantispike” reads by counting how many signature mutations
    match in a read’s alignment.

    Parameters
    ----------
    fname_quantispikes : str
        Path to a file describing quantispike amplicons.
    amplicons : AmpliconGroup
        The standard amplicon group for the main assay.
    min_match : int
        Minimum number of shared signature mutations to qualify a read
        as a definite quantispike.
    wt_max : int
        If shared mutations exceed this threshold but do not reach
        `min_match`, the read is tagged as "unclear_{quantispike_id}".
    """

    def __init__(self, fname_quantispikes, amplicons, min_match, wt_max):
        self._min_match = min_match
        self._wt_max = wt_max

        # Build an AmpliconGroup (or AmpliconSet) from the quantispike file
        self._quantispikes = AmpliconGroup.from_files(fname_quantispikes)

        # Precompute signature mutations for each quantispike
        self._mutations = self._build_filter(self._quantispikes, amplicons)

    def _build_filter(self, qspikes, amplicons):
        """
        Align each quantispike amplicon to the main amplicons,
        gather the resulting mutation sets.

        Returns
        -------
        defaultdict(list)
            Key: (amplicon_id, orientation)
            Value: list of (quantispike_id, mutation_set).
        """
        from .align_read import ReadAligner
        aligner = ReadAligner(amplicons)

        qs_mutations = defaultdict(list)

        for qs_amp in qspikes:
            for orientation in (True, False):
                aligned = aligner.align_read(qs_amp.get_construct(orientation))
                if aligned.matching_primers_found:
                    mut_set = aligned.mutations_not_filtered
                    if len(mut_set):
                        qs_mutations[aligned.amplicon_key].append(
                            (qs_amp.id, mut_set)
                        )
        return qs_mutations

    def detect(self, aligned_read):
        """
        Check whether the read is a quantispike read by looking
        at shared mutations with known quantispike signatures.

        Parameters
        ----------
        aligned_read : AlignedRead
            A read aligned to an amplicon.

        Returns
        -------
        str or None
            Returns the quantispike ID, possibly with prefix 'unclear_',
            or None if not a quantispike.
        """
        candidate_signatures = self._mutations[aligned_read.amplicon_key]
        read_muts = aligned_read.mutations_not_filtered

        max_shared = 0
        best_qs_id = None

        for (qs_id, qs_mutset) in candidate_signatures:
            shared_count = qs_mutset.count_same(read_muts)
            if shared_count > max_shared:
                max_shared = shared_count
                best_qs_id = qs_id

        # Decide if the read is definitely quantispike, unclear, or not at all
        if max_shared >= self._min_match:
            return best_qs_id
        elif max_shared > self._wt_max:
            return f"unclear_{best_qs_id}" if best_qs_id else None
        return None

    def get_all_sequences(self, add_reverse=True):
        """
        Gather all quantispike sequences.

        Parameters
        ----------
        add_reverse : bool
            If True, add the reverse complement as well.

        Returns
        -------
        list of str
        """
        return self._quantispikes.get_all_sequences(add_reverse)

    def get_amplicon_targets(self):
        """
        Return all (quantispike_id, amplicon_id) pairs found in the filter data.

        Returns
        -------
        list of (str, str)
        """
        pairs = []
        for (amp_id, orientation), entries in self._mutations.items():
            for (qs_id, _mutset) in entries:
                pairs.append((qs_id, amp_id))
        return sorted(set(pairs))


##############################################################################
# 3) AmpliconReadAlign
##############################################################################

class AmpliconReadAlign:
    """
    High-level read aligner that detects amplicon, pseudogene, and quantispike
    statuses for each read.

    Parameters
    ----------
    amplicons : AmpliconGroup
        The main set of reference amplicons for the assay.
    pseudogene_filter : PseudogeneFilter
        A filter that detects pseudogene reads.
    quantispike_filter : QuantispikeFilter
        A filter that detects quantispike reads.
    """

    def __init__(self, amplicons, pseudogene_filter, quantispike_filter):
        from .align_read import ReadAligner  # relative import to keep it consistent

        self._read_aligner = ReadAligner(amplicons)
        self._pseudogene_filter = pseudogene_filter
        self._quantispike_filter = quantispike_filter

    def align_read(self, read):
        """
        Align a single read to the amplicons. If it matches fully,
        check for pseudogene and quantispike signatures.

        Parameters
        ----------
        read : str
            Nucleotide sequence.

        Returns
        -------
        AlignedRead
            Describes the alignment, plus flags for pseudogene or quantispike.
        """
        aligned = self._read_aligner.align_read(read)

        # If the read does not match the amplicon or has other issues, just return
        if not aligned.is_amplicon_read():
            return aligned

        # Check for pseudogene
        if self._pseudogene_filter is not None:
            aligned.is_pseudogene = self._pseudogene_filter.detect(aligned)

        # Check for quantispike
        if self._quantispike_filter is not None:
            aligned.quantispike = self._quantispike_filter.detect(aligned)

        return aligned

    def align_reads(self, reads, cpu_limit):
        """
        Align a list of reads in parallel using multiple processes.

        Parameters
        ----------
        reads : list of str
            Nucleotide reads.
        cpu_limit : int
            Number of worker processes to use.

        Returns
        -------
        list of AlignedRead
        """
        with Pool(processes=cpu_limit) as pool:
            results = pool.map(self.align_read, reads)
        return results

