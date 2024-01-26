import pytest
from typing import List
from Bio import Align
from pytomebio.tools.change_seq.fastq_to_bam import gapped_alignment
from pytomebio.tools.change_seq.fastq_to_bam import perform_gapped_rescue
from pytomebio.core.attachment_site import AttachmentSite
from pytomebio.core.attachment_site import AttachmentSiteMatch
from pytomebio.core.aligner import get_query_prefix_aligner


@pytest.fixture(scope="session")
def query_prefix_aligner() -> Align.PairwiseAligner:
    yield get_query_prefix_aligner()


@pytest.mark.parametrize(
    "full_read, attachment_site, is_left, min_score, expected_alignment_length, "
    "expected_alignment_score",
    (
        # perfect match on left
        [
            "AA",
            AttachmentSite(name="test_site", left="AA", overhang="T", right="GG"),
            True,
            0,
            2,
            2,
        ],
        # perfect match on right
        [
            "GG",
            AttachmentSite(name="test_site", left="AA", overhang="T", right="GG"),
            False,
            0,
            2,
            2,
        ],
    ),
)
def test_gapped_aligner(
    query_prefix_aligner: Align.PairwiseAligner,
    full_read: str,
    attachment_site: AttachmentSite,
    is_left: bool,
    min_score: int,
    expected_alignment_length: int,
    expected_alignment_score: int,
) -> None:
    expected_attachment_site_match = AttachmentSiteMatch(
        site=attachment_site,
        is_left=is_left,
        score=expected_alignment_score,
        alignment_length=expected_alignment_length,
    )
    assert expected_attachment_site_match == gapped_alignment(
        full_read=full_read,
        site=attachment_site,
        is_left=is_left,
        aligner=query_prefix_aligner,
        min_score=min_score,
    )


@pytest.mark.parametrize(
    "full_read, attachment_site, is_left, min_score",
    (
        # no alignments produced, attempting to align left
        [
            "CC",
            AttachmentSite(name="test_site", left="AA", overhang="T", right="GG"),
            True,
            0,
        ],
        # no alignments produced, attempting to align right
        [
            "CC",
            AttachmentSite(name="test_site", left="AA", overhang="T", right="GG"),
            False,
            0,
        ],
        # produces alignment, but it's less than min score
        [
            "GGG",
            AttachmentSite(name="test_site", left="AAA", overhang="T", right="GG"),
            True,
            0,
        ],
        # produces alignment, but it's empty
        [
            "GGG",
            AttachmentSite(name="test_site", left="AAA", overhang="T", right="GG"),
            True,
            -9,
        ],
    ),
)
def test_gapped_aligner_returns_none(
    query_prefix_aligner: Align.PairwiseAligner,
    full_read: str,
    attachment_site: AttachmentSite,
    is_left: bool,
    min_score: int,
) -> None:
    assert (
        gapped_alignment(
            full_read=full_read,
            site=attachment_site,
            is_left=is_left,
            aligner=query_prefix_aligner,
            min_score=min_score,
        )
        is None
    )


@pytest.mark.parametrize(
    "full_read, attachment_site_list, rescue_prefix_len, expected_is_left, "
    "expected_alignment_length, expected_alignment_score, match_idx",
    (
        # Prefix matches left, aligns with bad score (aligned bases mismatch)
        [
            "GGGAAA",
            [AttachmentSite(name="match_left_site", left="GGGTTT", overhang="T", right="CCC")],
            3,
            True,
            3,
            -6,
            0,
        ],
        # Prefix matches right, aligns with bad score (aligned bases mismatch)
        [
            "CCCAAA",
            [AttachmentSite(name="match_right_site", left="GGG", overhang="T", right="CCCTTT")],
            3,
            False,
            3,
            -6,
            0,
        ],
        # Prefix matches 1st site left, aligns with bad score (aligned bases mismatch)
        [
            "GGGAAA",
            [
                AttachmentSite(
                    name="match_first_left_site", left="GGGTTT", overhang="T", right="CCC"
                ),
                AttachmentSite(name="test_site_2", left="CCCTTT", overhang="T", right="GGGTTT"),
            ],
            3,
            True,
            3,
            -6,
            0,
        ],
        [
            "CCCAAA",
            [
                AttachmentSite(
                    name="match_first_right_site", left="GGGTTT", overhang="T", right="CCCTTT"
                ),
                AttachmentSite(name="test_site_2", left="CCCTTT", overhang="T", right="GGGTTT"),
            ],
            3,
            False,
            3,
            -6,
            0,
        ],
        # Prefix matches 1st site right, aligns with bad score (aligned bases mismatch)
        [
            "CCCAAA",
            [
                AttachmentSite(
                    name="match_first_right_site", left="GGGTTT", overhang="T", right="CCCTTT"
                ),
                AttachmentSite(name="test_site_2", left="CCCTTT", overhang="T", right="GGGTTT"),
            ],
            3,
            False,
            3,
            -6,
            0,
        ],
        # Prefix matches 2nd site left, aligns with bad score (aligned bases mismatch)
        [
            "GGGAAA",
            [
                AttachmentSite(name="test_site_1", left="CCCTTT", overhang="T", right="CCC"),
                AttachmentSite(
                    name="match_second_left_site",
                    left="GGGTTT",
                    overhang="T",
                    right="AAATTT",
                ),
            ],
            3,
            True,
            3,
            -6,
            1,
        ],
        # Prefix matches 2nd site right, aligns with bad score (aligned bases mismatch)
        [
            "AAACCC",
            [
                AttachmentSite(name="test_site_1", left="CCCTTT", overhang="T", right="CCC"),
                AttachmentSite(
                    name="match_second_right_site",
                    left="GGGTTT",
                    overhang="T",
                    right="AAATTT",
                ),
            ],
            3,
            False,
            3,
            -6,
            1,
        ],
    ),
)
def test_rescue_gapped_alignment(
    query_prefix_aligner: Align.PairwiseAligner,
    full_read: str,
    attachment_site_list: List[AttachmentSite],
    rescue_prefix_len: int,
    expected_is_left: bool,
    expected_alignment_length: int,
    expected_alignment_score: int,
    match_idx: int,
) -> None:
    expected_attachment_site_match = AttachmentSiteMatch(
        site=attachment_site_list[match_idx],
        is_left=expected_is_left,
        score=expected_alignment_score,
        alignment_length=expected_alignment_length,
    )
    assert (
        perform_gapped_rescue(
            full_read=full_read,
            sites=attachment_site_list,
            aligner=query_prefix_aligner,
            rescue_prefix_len=rescue_prefix_len,
            min_score=-10,
        )
        == expected_attachment_site_match
    )


@pytest.mark.parametrize(
    "full_read, attachment_site_list, rescue_prefix_len, min_score",
    (
        # no attachment sites
        [
            "GGGAAA",
            [],
            3,
            0,
        ],
        # 1 attachment site, doesn't match
        [
            "GGGAAA",
            [AttachmentSite(name="prefix_mismatch", left="GGCTTT", overhang="T", right="CCC")],
            3,
            0,
        ],
        # 1 attachment site, prefix doesn't match
        [
            "GGGAAA",
            [AttachmentSite(name="prefix_mismatch", left="GGCTTT", overhang="T", right="CCC")],
            3,
            0,
        ],
        # 1 attachment site, score too low
        [
            "GGGAAA",
            [
                AttachmentSite(
                    name="score_less_than_min_score", left="GGGTTT", overhang="T", right="CCC"
                )
            ],
            3,
            0,
        ],
        # combos, just to test that None will be returned for lists longer than 1
        [
            "GGGAAA",
            [
                AttachmentSite(name="prefix_mismatch_1", left="GGCTTT", overhang="T", right="CCC"),
                AttachmentSite(
                    name="prefix_mismatch_2", left="AAATTT", overhang="T", right="TTTTTT"
                ),
            ],
            3,
            0,
        ],
        # combos, just to test that None will be returned for lists longer than 1
        # neither prefix matches
        [
            "GGGAAA",
            [
                AttachmentSite(name="prefix_mismatch_1", left="GGCTTT", overhang="T", right="CCC"),
                AttachmentSite(
                    name="prefix_mismatch_2", left="AAATTT", overhang="T", right="TTTTTT"
                ),
            ],
            3,
            0,
        ],
        # first prefix doesn't match, second has score less than min score
        [
            "GGGAAA",
            [
                AttachmentSite(name="prefix_mismatch_1", left="GGCTTT", overhang="T", right="CCC"),
                AttachmentSite(
                    name="score_less_than_min_score_2", left="AAA", overhang="T", right="GGGCCC"
                ),
            ],
            3,
            0,
        ],
        # first has score less than min socre, second prefix doesn't match
        [
            "GGGAAA",
            [
                AttachmentSite(
                    name="score_less_than_min_score_1", left="GGGTTT", overhang="T", right="CCC"
                ),
                AttachmentSite(
                    name="prefix_mismatch_2", left="AAATTT", overhang="T", right="TTTTTT"
                ),
            ],
            3,
            0,
        ],
        # both prefixes match, but score less than min score
        [
            "GGGAAA",
            [
                AttachmentSite(
                    name="score_less_than_min_score_1", left="GGGTTT", overhang="T", right="CCC"
                ),
                AttachmentSite(
                    name="score_less_than_min_score_2", left="AAA", overhang="T", right="GGGCCC"
                ),
            ],
            3,
            0,
        ],
    ),
)
def test_rescue_gapped_alignment_returns_none(
    query_prefix_aligner: Align.PairwiseAligner,
    full_read: str,
    attachment_site_list: List[AttachmentSite],
    rescue_prefix_len: int,
    min_score: int,
) -> None:
    assert (
        perform_gapped_rescue(
            full_read=full_read,
            sites=attachment_site_list,
            aligner=query_prefix_aligner,
            rescue_prefix_len=rescue_prefix_len,
            min_score=min_score,
        )
        is None
    )


def test_zero_rescue_prefix_len_raises_exception(
    query_prefix_aligner: Align.PairwiseAligner,
) -> None:
    """raises Exception (failed assert) because rescue_prefix_len == 0"""
    with pytest.raises(Exception):
        perform_gapped_rescue(
            full_read="A",
            sites=[AttachmentSite(name="test_site", left="T", overhang="A", right="G")],
            aligner=query_prefix_aligner,
            rescue_prefix_len=0,
            min_score=-10,
        )
