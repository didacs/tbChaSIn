from pathlib import Path
from typing import ClassVar
from typing import Dict
from typing import Iterator
from typing import List
from typing import Optional
from typing import Tuple

import attr
import numpy
import pysam
import pytest
from attr import frozen
from pysam import AlignedSegment
from pysam import AlignmentFile
from pysam import FastxRecord

from pytomebio.core.attachment_site import AttachmentSite
from pytomebio.core.attachment_site import AttachmentSiteMatch
from pytomebio.tools.change_seq.fastq_to_bam import MatchingMetric
from pytomebio.tools.change_seq.fastq_to_bam import UngappedAligner
from pytomebio.tools.change_seq.fastq_to_bam import fastq_to_bam


def _align(ungapped_aligner: UngappedAligner, full_read: str) -> Optional[AttachmentSiteMatch]:
    """Method for testing. Runs align on a string rather than a FastxRecord.

    Args:
        full_read: string of bases to align to attachment site
    """
    assert all(base in {"A", "T", "G", "C"} for base in full_read), "Read has illegal bases"

    rec = FastxRecord(sequence=full_read)
    return ungapped_aligner.align(rec)


@pytest.mark.parametrize(
    "s1, s2, expected_score", (["AA", "AA", 2], ["AA", "AT", 1], ["AA", "TT", 0])
)
def test_ungapped_aligner_score_method(s1: str, s2: str, expected_score: int) -> None:
    ungapped_aligner_obj = UngappedAligner(
        max_read_start_offset=0,
        max_target_start_offset=0,
        min_score=-10,
        match_score=1,
        mismatch_score=0,
        sites=[AttachmentSite(name="dummy site", left="A", overhang="T", right="A")],
    )
    assert ungapped_aligner_obj.score(s1=s1, s2=s2) == expected_score


@pytest.mark.parametrize(
    "max_read_start_offset, max_target_start_offset, attachment_site, "
    "seq_to_align, is_left, expected_score, alignment_length, read_start_offset",
    (  # perfect match
        [
            0,
            0,
            AttachmentSite("test_site", left="AA", overhang="T", right="TT"),
            "AA",
            True,
            2,
            2,
            0,
        ],
        [
            0,
            0,
            AttachmentSite("test_site", left="AA", overhang="T", right="TT"),
            "TT",
            False,
            2,
            2,
            0,
        ],
        # single mismatch
        [
            0,
            0,
            AttachmentSite("test_site", left="AAT", overhang="T", right="GG"),
            "AT",
            True,
            0,
            2,
            0,
        ],
        [
            0,
            0,
            AttachmentSite("test_site", left="AA", overhang="T", right="CGG"),
            "GG",
            False,
            0,
            2,
            0,
        ],
        # max_target_start_offset > 1
        [
            0,
            1,
            AttachmentSite("test_site", left="AAT", overhang="T", right="GG"),
            "AT",
            True,
            2,
            2,
            0,
        ],
        [
            0,
            1,
            AttachmentSite("test_site", left="AA", overhang="T", right="CGG"),
            "GG",
            False,
            2,
            2,
            0,
        ],
        # max_read_start_offset > 1
        [
            1,
            0,
            AttachmentSite("test_site", left="AA", overhang="T", right="CGG"),
            "TAA",
            True,
            2,
            2,
            1,
        ],
        [
            1,
            0,
            AttachmentSite("test_site", left="AA", overhang="T", right="CGA"),
            "TCG",
            False,
            2,
            2,
            1,
        ],
        # both max_target_offset and max_read_start_offset > 1
        [
            1,
            1,
            AttachmentSite("test_site", left="GAAA", overhang="T", right="CGA"),
            "CAAA",
            True,
            3,
            3,
            1,
        ],
        [
            1,
            1,
            AttachmentSite("test_site", left="GAAA", overhang="T", right="CGTA"),
            "TGTA",
            False,
            3,
            3,
            1,
        ],
        # choose longer of two alignments with the same score
        [
            0,
            0,
            AttachmentSite("test_site", left="ATAAAAA", overhang="T", right="AAAAA"),
            "AAAAAAA",
            True,
            5,
            7,
            0,
        ],
        [
            0,
            0,
            AttachmentSite("test_site", left="AAAAA", overhang="T", right="ATAAAAA"),
            "AAAAAAA",
            False,
            5,
            7,
            0,
        ],
    ),
)
def test_ungapped_aligner_align_method_single_attachment_site(
    max_read_start_offset: int,
    max_target_start_offset: int,
    attachment_site: AttachmentSite,
    seq_to_align: str,
    is_left: bool,
    expected_score: int,
    alignment_length: int,
    read_start_offset: int,
) -> None:
    test_ungapped_aligner = UngappedAligner(
        max_read_start_offset=max_read_start_offset,
        max_target_start_offset=max_target_start_offset,
        min_score=-10,
        match_score=1,
        mismatch_score=-1,
        sites=list([attachment_site]),
    )
    test_attachment_site_match = AttachmentSiteMatch(
        site=attachment_site,
        is_left=is_left,
        score=expected_score,
        alignment_length=alignment_length,
        read_start_offset=read_start_offset,
    )
    assert _align(test_ungapped_aligner, seq_to_align) == test_attachment_site_match


@pytest.mark.parametrize(
    "max_read_start_offset, max_target_start_offset, attachment_site_list, "
    "best_attachment_site_idx, seq_to_align, is_left, expected_score, alignment_length",
    (
        [
            0,
            0,
            list(
                [
                    AttachmentSite("test_site_1", left="AA", overhang="T", right="GG"),
                    AttachmentSite("test_site_2", left="TT", overhang="T", right="CC"),
                ]
            ),
            1,
            "TT",
            True,
            2,
            2,
        ],
        [
            0,
            0,
            list(
                [
                    AttachmentSite("test_site_1", left="AA", overhang="T", right="GG"),
                    AttachmentSite("test_site_2", left="TT", overhang="T", right="CC"),
                ]
            ),
            0,
            "AA",
            True,
            2,
            2,
        ],
        [
            0,
            0,
            list(
                [
                    AttachmentSite("test_site_1", left="AA", overhang="T", right="GG"),
                    AttachmentSite("test_site_2", left="TT", overhang="T", right="CC"),
                ]
            ),
            1,
            "CC",
            False,
            2,
            2,
        ],
        [
            0,
            0,
            list(
                [
                    AttachmentSite("test_site_1", left="AA", overhang="T", right="GG"),
                    AttachmentSite("test_site_2", left="TT", overhang="T", right="CC"),
                ]
            ),
            0,
            "GG",
            False,
            2,
            2,
        ],
    ),
)
def test_ungapped_aligner_multiple_attachment_sites(
    max_read_start_offset: int,
    max_target_start_offset: int,
    attachment_site_list: List[AttachmentSite],
    best_attachment_site_idx: int,
    seq_to_align: str,
    is_left: bool,
    expected_score: int,
    alignment_length: int,
) -> None:
    test_ungapped_aligner = UngappedAligner(
        max_read_start_offset=max_read_start_offset,
        max_target_start_offset=max_target_start_offset,
        min_score=-10,
        match_score=1,
        mismatch_score=-1,
        sites=attachment_site_list,
    )
    test_attachment_site_match = AttachmentSiteMatch(
        site=attachment_site_list[best_attachment_site_idx],
        is_left=is_left,
        score=expected_score,
        alignment_length=alignment_length,
    )
    assert _align(test_ungapped_aligner, seq_to_align) == test_attachment_site_match


@frozen
class QueryAlignment:
    """Class to most basic alignment info for a ReadPairTestCase. Also has code to assert two
    QueryAlignments are equal"""

    target: str
    is_left: bool

    def __str__(self) -> str:
        """String encoding used for writing fastq files"""
        return f"{self.target};{self.is_left}"

    @staticmethod
    def assert_alignments_equal(
        expected_alignment: Optional["QueryAlignment"],
        actual_alignment: Optional["QueryAlignment"],
    ) -> None:
        if expected_alignment is None:
            assert actual_alignment is None, (
                f"Expected no alignment but one was found {actual_alignment}",
            )
        else:
            assert actual_alignment is not None, (
                f"Expected alignment {actual_alignment} but no alignment was found",
            )
            assert expected_alignment.target == actual_alignment.target, (
                "aligned to wrong target",
            )
            assert expected_alignment.is_left == actual_alignment.is_left, (
                "aligned to wrong side",
            )


@frozen
class ReadPairTestCase:
    """
    Class to hold an individual test case (corresponding to a read pair and expected alignment and
    filter outcomes.
    Also has code to:
        -write to fastq files
        -reconstitute ReadPairTestCase objects from the actual outcomes in the output bam files
        -assert two ReadPairTestCases are equal
    """

    description: str = attr.field()
    query1: str = attr.ib()
    query2: str = attr.ib()
    accept: bool = attr.ib()
    rescue_attempted: bool = attr.ib()
    query_alignment1: Optional[QueryAlignment] = attr.ib()
    query_alignment2: Optional[QueryAlignment] = attr.ib()

    # Note underscore is disallowed because we convert spaces to underscore
    _disallowed_characters: ClassVar[str] = "\n\t_"

    def __attrs_post_init__(self) -> None:
        """Ensure description has SAM-compatible characters"""
        if any(
            disallowed in self.description
            for disallowed in ReadPairTestCase._disallowed_characters
        ):
            raise ValueError(f"description contains disallowed characters: {self.description}")

    @staticmethod
    def assert_test_cases_equal(
        expected_test_case: "ReadPairTestCase", actual_test_case: "ReadPairTestCase"
    ) -> None:
        assert expected_test_case.description == actual_test_case.description, (
            f"TestCase '{expected_test_case.description}' was missing from output bams or "
            f"'{actual_test_case.description}' was added"
        )
        assert expected_test_case.accept == actual_test_case.accept
        assert expected_test_case.rescue_attempted == actual_test_case.rescue_attempted
        # test that alignments go as expected
        QueryAlignment.assert_alignments_equal(
            expected_test_case.query_alignment1,
            actual_test_case.query_alignment1,
        )
        QueryAlignment.assert_alignments_equal(
            expected_test_case.query_alignment2,
            actual_test_case.query_alignment2,
        )

    def _get_fastq_comment(self, query_alignment: Optional[QueryAlignment]) -> str:
        """Get comment line for fastq record, using specified QueryAlignment"""
        front = f"{self.accept};{self.rescue_attempted}"
        return f"{front};None;None" if query_alignment is None else f"{front};{query_alignment}"

    @property
    def fastq_comment1(self) -> str:
        """Get comment line for fastq record for read1"""
        return self._get_fastq_comment(self.query_alignment1)

    @property
    def fastq_comment2(self) -> str:
        """Get comment line for fastq record for read2"""
        return self._get_fastq_comment(self.query_alignment2)

    @staticmethod
    def from_bams(keep_bam: Path, reject_bam: Path) -> List["ReadPairTestCase"]:
        """Given keep and reject BAMs, return list of ReadPairTestCases encoding the actual outcome
        of processing all the input fastqs"""
        return [
            ReadPairTestCase._from_sam_records(accept=True, rec1=rec1, rec2=rec2)
            for (rec1, rec2) in ReadPairTestCase._bam_to_rec_pairs(bam=keep_bam)
        ] + [
            ReadPairTestCase._from_sam_records(accept=False, rec1=rec1, rec2=rec2)
            for (rec1, rec2) in ReadPairTestCase._bam_to_rec_pairs(bam=reject_bam)
        ]

    @staticmethod
    def _from_sam_records(
        accept: bool, rec1: AlignedSegment, rec2: AlignedSegment
    ) -> "ReadPairTestCase":
        """
        Given two SAM records from a read pair, return a ReadPairTestCase that encodes the actual
        output of fastq_to_bam
        Args:
            accept: True if the record was in the keep BAM, false if it was in the reject bam.
            rec1: SAM record corresponding to read1
            rec2: SAM record corresponding to read2

        Returns:

        """
        query_name = rec1.query_name
        assert rec2.query_name == query_name
        description = query_name.replace("_", " ")
        assert rec1.is_read1
        assert rec2.is_read2
        sam_tag_1 = f"{rec1.get_tag(AttachmentSiteMatch.READ_MATCH_TAG)}"
        sam_tag_2 = f"{rec2.get_tag(AttachmentSiteMatch.READ_MATCH_TAG)}"
        assert rec1.get_tag(AttachmentSiteMatch.MATE_MATCH_TAG) == sam_tag_2
        assert rec2.get_tag(AttachmentSiteMatch.MATE_MATCH_TAG) == sam_tag_1
        target1, is_left1, read_start_offset1 = ReadPairTestCase._sam_tag_to_alignment_info(
            sam_tag_1
        )
        target2, is_left2, read_start_offset2 = ReadPairTestCase._sam_tag_to_alignment_info(
            sam_tag_2
        )
        return ReadPairTestCase(
            description=description,
            query1=rec1.query_sequence,
            query2=rec2.query_sequence,
            accept=accept,
            rescue_attempted=not (read_start_offset1 == 1 and read_start_offset2 == 1),
            query_alignment1=(
                None if target1 is None else QueryAlignment(target=target1, is_left=is_left1)
            ),
            query_alignment2=(
                None if target2 is None else QueryAlignment(target=target2, is_left=is_left2)
            ),
        )

    @staticmethod
    def _bam_to_rec_pairs(bam: Path) -> Iterator[Tuple[AlignedSegment, AlignedSegment]]:
        """Iterate over all read1, read2 read pairs in a BAM file"""
        with AlignmentFile(str(bam), "rb", threads=1, check_sq=False) as fh_in:
            rec1 = next(fh_in, None)
            while rec1 is not None:
                rec2 = next(fh_in)
                yield rec1, rec2
                rec1 = next(fh_in, None)

    @staticmethod
    def _sam_tag_to_alignment_info(
        sam_tag: str,
    ) -> Tuple[Optional[str], Optional[bool], Optional[int]]:
        """From a SAM tag extract (target, is_left, read_start_offset). If no alignment was found
        return None,None,None"""
        if sam_tag == "None":
            return None, None, None
        else:
            target, side, read_start_offset, __ = sam_tag.split(";", maxsplit=3)
            is_left = side == "left"
            return target, is_left, int(read_start_offset)


@frozen
class EndToEndFastqToBamData:
    """
    Class to manage overall test data:
    -hold test data raw inputs (attachment_sites, test_cases)
    -write input fastq files
    -set path for output bam files and metric TSV file
    -generate the expected metric based on input test cases
    """

    attachment_sites: list[AttachmentSite]
    test_cases: list[ReadPairTestCase]
    read1_fastq: Path
    read2_fastq: Path
    keep_bam: Path
    reject_bam: Path
    metric_tsv: Path

    @staticmethod
    def build(
        temp_path: Path, attachment_sites: list[AttachmentSite], test_cases: list[ReadPairTestCase]
    ) -> "EndToEndFastqToBamData":
        """Generate and store data for end-to-end test of fastq_to_bam inside temp_path"""
        input_path = temp_path / "input"
        input_path.mkdir(parents=True, exist_ok=False)
        output_path = temp_path / "output"
        output_path.mkdir(parents=False, exist_ok=False)
        read1_fastq = input_path / "read1.fastq.gz"
        read2_fastq = input_path / "read2.fastq.gz"

        names = [test_case.description for test_case in test_cases]
        if len(names) != len(set(names)):
            # raise error, we want TestCases to have unique descriptions, both because we
            # don't need to duplicate tests, and because we want to compare results to inputs
            name_counts: Dict[str, int] = {}
            for name in names:
                name_counts[name] = name_counts.get(name, 0) + 1
            for name, counts in name_counts.items():
                if counts > 1:
                    raise ValueError(
                        f"TestCase description '{name}' is duplicated ({counts} occurences)"
                    )

        EndToEndFastqToBamData._write_fastqs(
            fastq=read1_fastq,
            names=names,
            queries=[test_case.query1 for test_case in test_cases],
            comments=[test_case.fastq_comment1 for test_case in test_cases],
        )
        EndToEndFastqToBamData._write_fastqs(
            fastq=read2_fastq,
            names=names,
            queries=[test_case.query2 for test_case in test_cases],
            comments=[test_case.fastq_comment2 for test_case in test_cases],
        )

        return EndToEndFastqToBamData(
            attachment_sites=attachment_sites,
            test_cases=test_cases,
            read1_fastq=read1_fastq,
            read2_fastq=read2_fastq,
            keep_bam=output_path / "keep.bam",
            reject_bam=output_path / "reject.bam",
            metric_tsv=output_path / "metric.tsv",
        )

    @staticmethod
    def _write_fastqs(
        fastq: Path, names: List[str], queries: List[str], comments: List[str]
    ) -> None:
        with pysam.BGZFile(filename=f"{fastq}", mode="wb", index=None) as f_out:
            for name, query, comment in zip(names, queries, comments):
                f_out.write(
                    EndToEndFastqToBamData._query_to_fastq_bytes(
                        name=name, bases=query, comment=comment
                    )
                )

    @staticmethod
    def _query_to_fastq_bytes(name: str, bases: str, comment: Optional[str]) -> bytes:
        fastx_record = EndToEndFastqToBamData._query_to_fastx_record(
            name=name, bases=bases, comment=comment
        )
        return f"{fastx_record}\n".encode("utf-8")

    @staticmethod
    def _query_to_fastx_record(name: str, bases: str, comment: Optional[str]) -> FastxRecord:
        fastx_record = FastxRecord(
            name=name.replace(" ", "_"), sequence=bases, comment=comment, quality="#" * len(bases)
        )
        return fastx_record

    @property
    def expected_matching_metric(self) -> MatchingMetric:
        metric = MatchingMetric()
        for test_case in self.test_cases:
            metric.update(
                r1_match=EndToEndFastqToBamData._query_alignment_to_attachment_site_match(
                    test_case.query_alignment1
                ),
                r2_match=EndToEndFastqToBamData._query_alignment_to_attachment_site_match(
                    test_case.query_alignment2
                ),
                rescued_attempted=test_case.rescue_attempted,
            )
        return metric

    @staticmethod
    def _query_alignment_to_attachment_site_match(
        query_alignment: Optional[QueryAlignment],
    ) -> Optional[AttachmentSiteMatch]:
        """Return corresponding AttachmentSiteMatch. Most fields are unknown, this is intended for
        use only for updating MatchingMetric
        """
        return (
            None
            if query_alignment is None
            else AttachmentSiteMatch(
                site=AttachmentSite(name=query_alignment.target, left="", overhang="", right=""),
                is_left=query_alignment.is_left,
                score=numpy.NaN,
                alignment_length=-1,
                read_start_offset=-1,
            )
        )


def test_fastq_to_bam(tmp_path: Path) -> None:
    # not testing alignment, just integration into tool. So make a few basic sequences for reuse in
    # combos for different code paths
    left1 = "AAAAAAAAAACC"
    right1 = "CCCCCCCCCCAA"
    left2 = "GGGGGGGGGGTT"
    right2 = "TTTTTTTTTTGG"
    overhang = "GT"
    entropy = "GACTGATC"
    # mismatch by any aligner
    mismatch = entropy + entropy
    # perfect matches to front of target with offset = 1
    ungapped_left1 = right1[0] + left1 + entropy
    ungapped_right1 = left1[0] + right1 + entropy
    ungapped_left2 = right2[0] + left2 + entropy
    ungapped_right2 = left2[0] + right2 + entropy
    # have insertion that's too large for ungapped alignment. Front doesn't match any target prefix
    # so gapped rescue doesn't work
    gapped_left1 = left2[:2] + left1 + entropy
    gapped_right1 = right2[:2] + right1 + entropy
    gapped_left2 = left1[:2] + left2 + entropy
    gapped_right2 = right1[:2] + right2 + entropy
    # match perfectly at the front but have an insertion in middle that's too large for ungapped
    # alignment
    gapped_rescue_left1 = left1[:3] + left2[:3] + left1[3:] + entropy
    gapped_rescue_right1 = right1[:3] + right2[:3] + right1[3:] + entropy
    gapped_rescue_left2 = left2[:3] + left1[:3] + left2[3:] + entropy
    gapped_rescue_right2 = right2[:3] + right1[:3] + right2[3:] + entropy

    end_to_end_fastq_to_bam_data = EndToEndFastqToBamData.build(
        temp_path=tmp_path,
        attachment_sites=[
            AttachmentSite(name="1", left=left1, overhang=overhang, right=right1),
            AttachmentSite(name="2", left=left2, overhang=overhang, right=right2),
        ],
        test_cases=[
            # UNGAPPED TEST CASES
            # accept: ungapped alignment match 1/left, 1/right
            ReadPairTestCase(
                description="accept: ungapped alignment match 1/left, 1/right",
                query1=ungapped_left1,
                query2=ungapped_right1,
                accept=True,
                rescue_attempted=False,
                query_alignment1=QueryAlignment(target="1", is_left=True),
                query_alignment2=QueryAlignment(target="1", is_left=False),
            ),
            # accept: ungapped alignment match 2/right, 2/left
            ReadPairTestCase(
                description="accept: ungapped alignment match 2/right, 2/left",
                query1=ungapped_right2,
                query2=ungapped_left2,
                accept=True,
                rescue_attempted=False,
                query_alignment1=QueryAlignment(target="2", is_left=False),
                query_alignment2=QueryAlignment(target="2", is_left=True),
            ),
            # reject: (same side) ungapped alignment match 1/right, 1/right
            ReadPairTestCase(
                description="reject: (same side) ungapped alignment match 1/right, 1/right",
                query1=ungapped_right1,
                query2=ungapped_right1,
                accept=False,
                rescue_attempted=False,
                query_alignment1=QueryAlignment(target="1", is_left=False),
                query_alignment2=QueryAlignment(target="1", is_left=False),
            ),
            # reject: (different site) ungapped alignment match 2/right, 2/left
            ReadPairTestCase(
                description="reject: (different site) ungapped alignment match 2/right, 2/left",
                query1=ungapped_right2,
                query2=ungapped_left1,
                accept=False,
                rescue_attempted=False,
                query_alignment1=QueryAlignment(target="2", is_left=False),
                query_alignment2=QueryAlignment(target="1", is_left=True),
            ),
            # reject: (same side and different site) ungapped alignment match 2/left, 1/left
            ReadPairTestCase(
                description="reject: (same side and different site) ungapped alignment match "
                "2/left, 1/left",
                query1=ungapped_left2,
                query2=ungapped_left1,
                accept=False,
                rescue_attempted=False,
                query_alignment1=QueryAlignment(target="2", is_left=True),
                query_alignment2=QueryAlignment(target="1", is_left=True),
            ),
            # GAPPED (NOT RESCUE) TEST CASES
            # accept: ungapped query 1, gapped alignment for query 2
            ReadPairTestCase(
                description="accept: ungapped query 1, gapped alignment for query 2",
                query1=ungapped_left1,
                query2=gapped_right1,
                accept=True,
                rescue_attempted=True,
                query_alignment1=QueryAlignment(target="1", is_left=True),
                query_alignment2=QueryAlignment(target="1", is_left=False),
            ),
            # accept: gapped alignment for query 1, ungapped for query 2
            ReadPairTestCase(
                description="accept: gapped alignment for query 1, ungapped for query 2",
                query1=gapped_right2,
                query2=ungapped_left2,
                accept=True,
                rescue_attempted=True,
                query_alignment1=QueryAlignment(target="2", is_left=False),
                query_alignment2=QueryAlignment(target="2", is_left=True),
            ),
            # reject: ungapped query 1,
            #         would have gapped alignment for query 2, but it's not the correct side so
            #         gapped alignment fails. No rescue is possible because prefix doesn't match
            ReadPairTestCase(
                description="reject: ungapped query 1, gapped same-side query 2",
                query1=ungapped_left1,
                query2=gapped_left1,
                accept=False,
                rescue_attempted=True,
                query_alignment1=QueryAlignment(target="1", is_left=True),
                query_alignment2=None,
            ),
            # reject: ungapped query 1,
            #         would have gapped alignment for query 2, but it's the wrong site so
            #         gapped alignment fails. No rescue is possible because prefix doesn't match
            ReadPairTestCase(
                description="reject: ungapped query 1, gapped different site query 2",
                query1=ungapped_left1,
                query2=gapped_right2,
                accept=False,
                rescue_attempted=True,
                query_alignment1=QueryAlignment(target="1", is_left=True),
                query_alignment2=None,
            ),
            # reject: ungapped query 2,
            #         would have gapped alignment for query 1, but it's not the correct side so
            #         gapped alignment fails. No rescue is possible because prefix doesn't match
            ReadPairTestCase(
                description="reject: ungapped query 2, gapped same-side query 1",
                query1=gapped_left2,
                query2=ungapped_left2,
                accept=False,
                rescue_attempted=True,
                query_alignment1=None,
                query_alignment2=QueryAlignment(target="2", is_left=True),
            ),
            # reject: ungapped query 2,
            #         would have gapped alignment for query 1, but it's the wrong site so
            #         gapped alignment fails. No rescue is possible because prefix doesn't match
            ReadPairTestCase(
                description="reject: ungapped query 2, gapped different site query 1",
                query1=ungapped_right1,
                query2=gapped_left2,
                accept=False,
                rescue_attempted=True,
                query_alignment1=QueryAlignment(target="1", is_left=False),
                query_alignment2=None,
            ),
            # 1 UNGAPPED x 1 GAPPED RESCUE TEST CASES (I THINK THIS ARRANGEMENT MUST ALWAYS FAIL)
            # reject: ungapped query 1,
            #         gapped (rescue) alignment for query 2, but it's not the correct side
            ReadPairTestCase(
                description="reject: ungapped query 1, rescue query 2 on same site",
                query1=ungapped_left2,
                query2=gapped_rescue_left2,
                accept=False,
                rescue_attempted=True,
                query_alignment1=QueryAlignment(target="2", is_left=True),
                query_alignment2=QueryAlignment(target="2", is_left=True),
            ),
            # reject: ungapped query 1,
            #         gapped (rescue) alignment for query 2, but it's the wrong site
            ReadPairTestCase(
                description="reject: ungapped query 1, rescue query 2 on different site",
                query1=ungapped_left1,
                query2=gapped_rescue_right2,
                accept=False,
                rescue_attempted=True,
                query_alignment1=QueryAlignment(target="1", is_left=True),
                query_alignment2=QueryAlignment(target="2", is_left=False),
            ),
            # reject: ungapped query 2,
            #         gapped (rescue) alignment for query 1, but it's not the correct side
            ReadPairTestCase(
                description="reject: ungapped query 2, rescue query 1 on same side",
                query1=gapped_rescue_right2,
                query2=ungapped_right2,
                accept=False,
                rescue_attempted=True,
                query_alignment1=QueryAlignment(target="2", is_left=False),
                query_alignment2=QueryAlignment(target="2", is_left=False),
            ),
            # reject: ungapped query 1,
            #         gapped (rescue) alignment for query 1, but it's the wrong site
            ReadPairTestCase(
                description="reject: ungapped query 1, query 2 on different site",
                query1=ungapped_right1,
                query2=gapped_rescue_left2,
                accept=False,
                rescue_attempted=True,
                query_alignment1=QueryAlignment(target="1", is_left=False),
                query_alignment2=QueryAlignment(target="2", is_left=True),
            ),
            # TWO GAPPED RESCUE TEST CASES
            # accept: gapped rescue alignment match 1/right, 1/left
            ReadPairTestCase(
                description="accept: gapped rescue alignment match 1/right, 1/left",
                query1=gapped_rescue_right1,
                query2=gapped_rescue_left1,
                accept=True,
                rescue_attempted=True,
                query_alignment1=QueryAlignment(target="1", is_left=False),
                query_alignment2=QueryAlignment(target="1", is_left=True),
            ),
            # accept: gapped rescue alignment match 2/left, 2/right
            ReadPairTestCase(
                description="accept: gapped rescue alignment match 2/left, 2/right",
                query1=gapped_rescue_left2,
                query2=gapped_rescue_right2,
                accept=True,
                rescue_attempted=True,
                query_alignment1=QueryAlignment(target="2", is_left=True),
                query_alignment2=QueryAlignment(target="2", is_left=False),
            ),
            # reject: (same side) gapped rescue alignment match 1/right, 1/right
            ReadPairTestCase(
                description="reject: (same side) gapped rescue alignment match 1/right, 1/right",
                query1=gapped_rescue_right1,
                query2=gapped_rescue_right1,
                accept=False,
                rescue_attempted=True,
                query_alignment1=QueryAlignment(target="1", is_left=False),
                query_alignment2=QueryAlignment(target="1", is_left=False),
            ),
            # reject: (different site) gapped rescue alignment match 2/right, 2/left
            ReadPairTestCase(
                description="reject: (different site) gapped rescue alignment match 2/right, "
                "2/left",
                query1=gapped_rescue_right2,
                query2=gapped_rescue_left1,
                accept=False,
                rescue_attempted=True,
                query_alignment1=QueryAlignment(target="2", is_left=False),
                query_alignment2=QueryAlignment(target="1", is_left=True),
            ),
            # reject: (same side and different site) gapped rescue alignment match 2/left, 1/left
            ReadPairTestCase(
                description="reject: (same side and different site) gapped rescue alignment match "
                "2/left, 1/left",
                query1=gapped_rescue_left2,
                query2=gapped_rescue_left1,
                accept=False,
                rescue_attempted=True,
                query_alignment1=QueryAlignment(target="2", is_left=True),
                query_alignment2=QueryAlignment(target="1", is_left=True),
            ),
            # 1 MISMATCH 1 UNGAPPED ALIGNMENT TEST CASES
            # reject: query1 mismatch, query2 ungapped match 1/left
            ReadPairTestCase(
                description="reject: query1 mismatch, query2 ungapped match 1/left",
                query1=mismatch,
                query2=ungapped_left1,
                accept=False,
                rescue_attempted=True,
                query_alignment1=None,
                query_alignment2=QueryAlignment(target="1", is_left=True),
            ),
            # reject: query1 mismatch, query 2 ungapped match 2/right
            ReadPairTestCase(
                description="reject: query1 mismatch, query2 ungapped match 1/right",
                query1=mismatch,
                query2=ungapped_right2,
                accept=False,
                rescue_attempted=True,
                query_alignment1=None,
                query_alignment2=QueryAlignment(target="2", is_left=False),
            ),
            # reject: query2 mismatch, query1 ungapped match 1/right
            ReadPairTestCase(
                description="reject: query2 mismatch, query1 ungapped match 1/right",
                query1=ungapped_right1,
                query2=mismatch,
                accept=False,
                rescue_attempted=True,
                query_alignment1=QueryAlignment(target="1", is_left=False),
                query_alignment2=None,
            ),
            # reject: query2 mismatch, query 1 ungapped match 2/left
            ReadPairTestCase(
                description="reject: query2 mismatch, query 1 ungapped match 2/left",
                query1=ungapped_left2,
                query2=mismatch,
                accept=False,
                rescue_attempted=True,
                query_alignment1=QueryAlignment(target="2", is_left=True),
                query_alignment2=None,
            ),
            # 1 MISMATCH 1 GAPPED RESCUE ALIGNMENT TEST CASES
            # reject: query1 mismatch, query2 gapped rescue match 1/left
            ReadPairTestCase(
                description="reject: query1 mismatch, query2 gapped rescue match 1/left",
                query1=mismatch,
                query2=gapped_rescue_left1,
                accept=False,
                rescue_attempted=True,
                query_alignment1=None,
                query_alignment2=QueryAlignment(target="1", is_left=True),
            ),
            # reject: query1 mismatch, query 2 gapped rescue match 2/right
            ReadPairTestCase(
                description="reject: query1 mismatch, query 2 gapped rescue match 2/right",
                query1=mismatch,
                query2=gapped_rescue_right2,
                accept=False,
                rescue_attempted=True,
                query_alignment1=None,
                query_alignment2=QueryAlignment(target="2", is_left=False),
            ),
            # reject: query2 mismatch, query1 gapped rescue match 1/right
            ReadPairTestCase(
                description="reject: query2 mismatch, query1 gapped rescue match 1/right",
                query1=gapped_rescue_right1,
                query2=mismatch,
                accept=False,
                rescue_attempted=True,
                query_alignment1=QueryAlignment(target="1", is_left=False),
                query_alignment2=None,
            ),
            # reject: query2 mismatch, query 1 gapped rescue match 2/left
            ReadPairTestCase(
                description="reject: query2 mismatch, query 1 gapped rescue match 2/left",
                query1=gapped_rescue_left2,
                query2=mismatch,
                accept=False,
                rescue_attempted=True,
                query_alignment1=QueryAlignment(target="2", is_left=True),
                query_alignment2=None,
            ),
            # 2 MISMATCH TEST CASES
            ReadPairTestCase(
                description="reject: two mismatches",
                query1=mismatch,
                query2=mismatch,
                accept=False,
                rescue_attempted=True,
                query_alignment1=None,
                query_alignment2=None,
            ),
        ],
    )

    matching_metric: MatchingMetric = fastq_to_bam(
        r1_fq=end_to_end_fastq_to_bam_data.read1_fastq,
        r2_fq=end_to_end_fastq_to_bam_data.read2_fastq,
        keep_bam=end_to_end_fastq_to_bam_data.keep_bam,
        reject_bam=end_to_end_fastq_to_bam_data.reject_bam,
        metric_tsv=end_to_end_fastq_to_bam_data.metric_tsv,
        read_group="test_fastq_to_bam",
        attachment_site=end_to_end_fastq_to_bam_data.attachment_sites,
        _min_score=3,
        _rescue_prefix_len=3,
        _max_read_start_offset=1,
        _max_target_start_offset=1,
    )

    # get the input test cases, sorted by description
    original_test_cases = sorted(
        end_to_end_fastq_to_bam_data.test_cases,
        key=lambda test_case: test_case.description,
    )

    # get the actual outcomes of the test cases, sorted by description
    processed_test_cases = sorted(
        ReadPairTestCase.from_bams(
            keep_bam=end_to_end_fastq_to_bam_data.keep_bam,
            reject_bam=end_to_end_fastq_to_bam_data.reject_bam,
        ),
        key=lambda test_case: test_case.description,
    )

    # assert they're equal
    for original_test_case, processed_test_case in zip(original_test_cases, processed_test_cases):
        ReadPairTestCase.assert_test_cases_equal(
            expected_test_case=original_test_case, actual_test_case=processed_test_case
        )

    # assert matching metric is equal
    assert matching_metric == end_to_end_fastq_to_bam_data.expected_matching_metric
