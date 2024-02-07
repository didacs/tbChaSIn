from typing import Optional

import pytest
from attr import frozen
from Bio import Align

from pytomebio.core.aligner import get_glocal_aligner
from pytomebio.core.aligner import get_query_prefix_aligner


@frozen
class AlignmentForTest:
    """Stores the output of aligner for testing

    Attributes:
        query: sequence of query
        target: sequence of target
        query_start_index: index of first base of the query (read) with a match to the
                           target (attachment site).
        target_start_index: index of the first base of the target (read) with a match to the query.
        best_score: score of the best alignment
    """

    query: str
    target: str
    query_start_index: int
    target_start_index: int
    best_score: int

    @classmethod
    def _build(
        cls, alignments: Align.PairwiseAlignments, query: str, target: str
    ) -> Optional["AlignmentForTest"]:
        """Builds object to store output of aligner
        Args:
            aligner: pairwise aligner to test
            query: sequence of query
            target: sequence of target

        Returns:
            TestAlignment object
        """
        # keep the best alignment, which we assume is the first alignment
        # we don't care about ties
        return AlignmentForTest(
            query=query,
            target=target,
            query_start_index=alignments[0].aligned[0][0][0],
            target_start_index=alignments[0].aligned[1][0][0],
            best_score=alignments[0].score,
        )

    @classmethod
    def build_query_prefix(cls, query: str, target: str) -> Optional["AlignmentForTest"]:
        """Builds object to store output of query_prefix_aligner (prefix of query to full target)

        Args:
            query: sequence of query
            target: sequence of target

        Returns:
            TestAlignment object
        """
        aligner = get_query_prefix_aligner()
        alignments = aligner.align(query, target)
        if alignments[0].aligned.size == 0:
            return None
        return cls._build(alignments=alignments, query=query, target=target)

    @classmethod
    def build_glocal(cls, query: str, target: str) -> Optional["AlignmentForTest"]:
        """Builds object to store output of glocal_aligner (subset of query to full target)

        Args:
            query: sequence of query
            target: sequence of target

        Returns:
            TestAlignment object
        """
        aligner = get_glocal_aligner()
        alignments = aligner.align(query, target)
        if len(alignments) == 0:
            return None
        return cls._build(alignments=alignments, query=query, target=target)


@pytest.mark.parametrize(
    "query, target,expected_score, expected_query_index, expected_target_index",
    [
        ("AAAGGGG", "AAA", 3, 0, 0),  # 3=
        ("AAAAGGGG", "AATA", -1, 0, 0),  # 2=1X1=
        ("AAAAGGGG", "ATTA", -6, 0, 0),  # 1=2X1=
        ("AAAAGGGG", "ATTT", -8, 0, 0),  # 1=3D
        ("AGGGTTTGGGG", "ATTT", -5, 0, 0),  # 1=3I3=
        ("AATGCAAAA", "AAG", -2, 0, 0),  # 2=1X
        ("AATGCAAAA", "AAC", -2, 0, 0),  # 2=1X
        ("AGCATTAGGGGGGGG", "ATGCATTAG", 1, 0, 0),  # 1=1D7=
        ("ATTGCATTAGGGGGGGG", "ATGCATTAG", 2, 0, 0),  # 2=1I7=
        ("ACCATTAGGGGGGGG", "ATGCATTAG", -4, 0, 0),  # 1=1D1X6=
        ("AATTAGGGGGGGGGG", "ATGCATTAG", -3, 0, 0),  # 1=3D5=
    ],
)
def test_query_prefix_alignment_scores_and_indices(
    query: str,
    target: str,
    expected_score: int,
    expected_query_index: int,
    expected_target_index: int,
) -> None:
    assert AlignmentForTest.build_query_prefix(target=target, query=query) == AlignmentForTest(
        query=query,
        target=target,
        query_start_index=expected_query_index,
        target_start_index=expected_target_index,
        best_score=expected_score,
    )


# check alignments that should return None
@pytest.mark.parametrize("query, target", [("AAA", "GGG"), ("CCCCCCCCC", "GCAAA")])
def test_no_alignment_prefix_query(query: str, target: str) -> None:
    assert AlignmentForTest.build_query_prefix(target=target, query=query) is None


@pytest.mark.parametrize(
    "query, target,expected_score, expected_query_index, expected_target_index",
    [
        # section of query aligns perfectly to target
        ("AAAGGGG", "AAA", 3, 0, 0),  # 3=
        ("CCCCAAA", "AAA", 3, 4, 0),  # 3=
        ("CCCCAAAGGGG", "AAA", 3, 4, 0),  # 3=
        # section of query two matches (not three matches and one mismatch) with target
        # soft-clips added to CIGAR comments because I was confused at first by the alignments
        ("AAAAGGGG", "AATA", 2, 0, 0),  # 2=6S
        ("CCCCAAAAGGGG", "AATA", 2, 4, 0),  # 4S2=6S
        ("CCCCAAAA", "AATA", 2, 4, 0),  # 4S2=2S
        # query has one mismatch with target (not four matches, length is tie-breaker)
        ("CCAAAAGGGGGCC", "AAAATGGGG", 4, 2, 0),  # 4=1X4=
        # query has one extra base
        ("CCGGAAGGATCAACCAACC", "GGAAGGACAACCAA", 7, 2, 0),  # 7=1I7=
        # query has two extra bases
        ("CCAGGAAGGATTCAACCAACCC", "AGGAAGGACAACCAAC", 8, 2, 0),  # 8=2I8=
        # query has one missing base
        ("CCGGAAGGACAACCAACC", "GGAAGGATCAACCAA", 7, 2, 0),  # 7=1D7=
        # query has two missing bases
        ("CCAGGAAGGACAACCAACCC", "AGGAAGGATTCAACCAAC", 8, 2, 0),  # 8=2D8=
    ],
)
def test_glocal_alignment_scores_and_indices(
    query: str,
    target: str,
    expected_score: int,
    expected_query_index: int,
    expected_target_index: int,
) -> None:
    assert AlignmentForTest.build_glocal(target=target, query=query) == AlignmentForTest(
        query=query,
        target=target,
        query_start_index=expected_query_index,
        target_start_index=expected_target_index,
        best_score=expected_score,
    )


@pytest.mark.parametrize("query, target", [("AA", "GG"), ("CC", "TT"), ("AAA", "GTC")])
def test_no_alignment_glocal(query: str, target: str) -> None:
    assert AlignmentForTest.build_glocal(target=target, query=query) is None
