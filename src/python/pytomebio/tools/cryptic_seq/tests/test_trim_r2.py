import collections
import typing
from pathlib import Path
from typing import Dict
from typing import Iterator
from typing import List

from attr import frozen
from fgpyo import sequence
from fgpyo.sam import Template
from fgpyo.sam.builder import SamBuilder
from pysam import AlignmentFile

from pytomebio.tools.common import TN5_ME
from pytomebio.tools.cryptic_seq.trim_r2 import ScoreCountsMetric
from pytomebio.tools.cryptic_seq.trim_r2 import trim_r2


@frozen
class TrimR2TestCase:
    """
    Holds an individual test case and expected/actual outcome from the trim_r2 tool
    Also has code to:
        - reconstitute TrimR2TestCase objects from the actual outcomes in the output bam
          files
        - assert two TrimR2TestCase are equal

    Attributes:
        name: description of test case (no whitespace other than spaces)
        bases1: bases on read1. Should not be examined by the tool so defaults to a simple default
                string. In a few test cases we will override this to test that bases1 are ignored
                even if they match the adapter sequence
        aln_score: alignment score when aligned to the adapter sequence
        trimmed_len: number of bases trimmed from R2
        keep: whether the read pair is accepted or rejected
        trim: whether R2 is trimmed
        bases2: bases on read2
    """

    name: str
    bases2: str
    aln_score: int
    trimmed_len: int
    keep: bool
    trim: bool
    bases1: str = "GATCGATCGATC"

    @property
    def sam_name(self) -> str:
        """return a SAM-compatible name"""
        return self.name.replace(" ", "_")

    @staticmethod
    def assert_equal(expected: "TrimR2TestCase", actual: "TrimR2TestCase") -> None:
        """To be equivalent, everything must be equal, except that num-mismatches may differ if
        the test case is not accepted (because the mismatches tag isn't added to the reject bam,
        so that information is lost)
        """
        assert expected.name == actual.name
        assert expected.bases1 == actual.bases1, f"'{expected.name}' has wrong bases1"
        assert expected.bases2 == actual.bases2, f"'{expected.name}' has wrong bases2"
        assert expected.keep == actual.keep, f"'{expected.name}' has wrong keep"
        if expected.keep:
            # if the test case is not accepted, it's not possible to reconstruct the number of
            # mismatches. However, this is still tested (indirectly) by comparing the metrics
            assert expected.aln_score == actual.aln_score, f"'{expected.name}' has wrong aln_score"
            assert (
                expected.trimmed_len == actual.trimmed_len
            ), f"'{expected.name}' has wrong trimmed_len"

    @staticmethod
    def from_bams(
        keep_bam: Path,
        reject_bam: Path,
        adapter_len: int,
    ) -> List["TrimR2TestCase"]:
        """Reconstitute test cases from keep and reject BAMs output by trim_for_tn5me"""
        return [
            TrimR2TestCase._from_template(keep=True, template=template, adapter_len=adapter_len)
            for template in TrimR2TestCase._bam_to_template_iterator(keep_bam)
        ] + [
            TrimR2TestCase._from_template(keep=False, template=template, adapter_len=adapter_len)
            for template in TrimR2TestCase._bam_to_template_iterator(reject_bam)
        ]

    @staticmethod
    def _bam_to_template_iterator(bam: Path) -> Iterator[Template]:
        """Use TemplateIterator to iterate over all read1, read2 read pairs in a BAM file"""
        with AlignmentFile(f"{bam}", check_sq=False) as fh_in:
            yield from Template.iterator(fh_in)

    @staticmethod
    def _from_template(
        keep: bool,
        template: Template,
        adapter_len: int,
    ) -> "TrimR2TestCase":
        """Reconstruct a  TrimR2TestCase from read pair, and validate tag correctness"""
        template.validate()
        human_readable_name = f"{template.name}".replace("_", " ")
        assert template.r1 is not None
        assert template.r1.is_read1
        assert template.r2 is not None
        assert template.r2.is_read2
        # read 1 should never be tagged
        assert not template.r1.has_tag("t3")
        assert not template.r1.has_tag("tm")
        assert not template.r1.has_tag("ts")
        if keep:
            trimmed = template.r2.has_tag("t3")
            if trimmed:
                assert template.r2.has_tag("tm")
                assert template.r2.has_tag("ts")
                t3_tag = typing.cast(str, template.r2.get_tag("t3"))
                bases2 = template.r2.query_sequence + t3_tag
                aln_score = int(template.r2.get_tag("ts"))
                trimmed_len = len(t3_tag)
            else:
                bases2 = template.r2.query_sequence
                aln_score = -1
                trimmed_len = len(bases2)
            return TrimR2TestCase(
                name=human_readable_name,
                bases1=template.r1.query_sequence,
                bases2=bases2,
                aln_score=aln_score,
                trimmed_len=trimmed_len,
                keep=True,
                trim=trimmed,
            )
        else:
            # read 1 is neither tagged nor trimmed if it is accepted
            assert not template.r2.has_tag("t3")
            assert not template.r2.has_tag("tm")
            assert not template.r2.has_tag("ts")
            return TrimR2TestCase(
                name=human_readable_name,
                bases1=template.r1.query_sequence,
                bases2=template.r2.query_sequence,
                aln_score=-1,
                trimmed_len=-1,
                keep=False,
                trim=False,
            )


@frozen
class TrimR2TestBuilder:
    """
    Container for test cases.
    - hold test data raw inputs (adapter_seq, min_score, test_cases)
    - write input bam file
    - set path for output bam files and metric TSV file
    - generate the expected metric based on input test cases
    """

    adapter_seq: str
    min_score: int
    test_cases: List[TrimR2TestCase]
    in_bam: Path
    keep_bam: Path
    reject_bam: Path
    metrics_tsv: Path

    @staticmethod
    def from_test_cases(
        test_dir: Path,
        adapter_seq: str,
        min_score: int,
        test_cases: List[TrimR2TestCase],
    ) -> "TrimR2TestBuilder":
        """Generate and store data for end-to-end test of cryptic_seq/trim_r2 inside
        test_dir

        Args:
            test_dir: path to directory that should contain all the temporary files for test
            adapter_seq: the adapter sequence used in the test
            min_score: the minimum alignment score to allow for passing test cases
            test_cases: list of TrimR2TestCase that will comprehensively test logic flow of
                        trim_r2 tool
        """
        TrimR2TestBuilder.assert_names_unique(test_cases)

        in_bam = test_dir / "in.bam"
        input_bam_builder = SamBuilder()
        for test_case in test_cases:
            # bases2 are never examined, so generate them randomly
            input_bam_builder.add_pair(
                name=test_case.sam_name, bases1=test_case.bases1, bases2=test_case.bases2
            )
        input_bam_builder.to_path(path=in_bam)
        return TrimR2TestBuilder(
            adapter_seq=adapter_seq,
            min_score=min_score,
            test_cases=test_cases,
            in_bam=in_bam,
            keep_bam=test_dir / "keep.bam",
            reject_bam=test_dir / "reject.bam",
            metrics_tsv=test_dir / "mismatch_metrics.tsv",
        )

    @staticmethod
    def assert_names_unique(test_cases: List[TrimR2TestCase]) -> None:
        """Raise ValueError if any test cases have duplicate names"""
        names = [test_case.name for test_case in test_cases]
        if len(names) != len(set(names)):
            # raise error, we want TestCases to have unique descriptions, both because we
            # don't need to duplicate tests, and because we want to compare results to inputs
            name_counts: Dict[str, int] = {}
            for name in names:
                name_counts[name] = name_counts.get(name, 0) + 1
            for name, counts in name_counts.items():
                if counts > 1:
                    raise ValueError(
                        f"TestCase name '{name}' is duplicated ({counts} occurrences)"
                    )

    @property
    def expected_mismatch_metric(self) -> Dict[int, ScoreCountsMetric]:
        """Generate expected Dict of Tn5ReadOneMatchingMetric based on input test_cases"""
        counter = collections.Counter({-1: 0})
        for test_case in self.test_cases:
            if test_case.trim:
                counter[test_case.aln_score] += 1
            else:
                counter[-1] += 1
        return ScoreCountsMetric.from_counter(counter=counter)


def test_trim_r2_adapter(
    tmp_path: Path,
    adapter_seq: str = TN5_ME,
    min_score: int = 14,
    min_match_length: int = 18,
) -> None:
    SEQ_LEN = 12
    keep_seq = "A" * SEQ_LEN
    test_builder = TrimR2TestBuilder.from_test_cases(
        adapter_seq=adapter_seq,
        min_score=min_score,
        test_dir=tmp_path,
        test_cases=[
            # accept 0 mismatches
            TrimR2TestCase(
                name="Trim:perfect match",
                bases2=keep_seq + sequence.reverse_complement(adapter_seq),
                aln_score=19,
                trimmed_len=19,
                keep=True,
                trim=True,
            ),
            # accept 1 mismatch at beginning
            TrimR2TestCase(
                name="Trim:1 mismatch at beginning",
                bases2=keep_seq + sequence.reverse_complement("C" + adapter_seq[1:]),
                aln_score=14,
                trimmed_len=19,
                keep=True,
                trim=True,
            ),
            # accept 1 mismatch in middle
            TrimR2TestCase(
                name="Trim:1 mismatch in middle",
                bases2=keep_seq
                + sequence.reverse_complement(adapter_seq[:4] + "C" + adapter_seq[5:]),
                aln_score=14,
                trimmed_len=19,
                keep=True,
                trim=True,
            ),
            # accept 1 mismatch at end
            TrimR2TestCase(
                name="Trim:1 mismatch at end",
                bases2=keep_seq + sequence.reverse_complement(adapter_seq[:-1] + "A"),
                aln_score=14,
                trimmed_len=19,
                keep=True,
                trim=True,
            ),
            # don't trim 2 mismatches at beginning
            TrimR2TestCase(
                name="NoTrim:2 mismatches at beginning",
                bases2=keep_seq + sequence.reverse_complement("CC" + adapter_seq[2:]),
                aln_score=-1,
                trimmed_len=31,
                keep=True,
                trim=False,
            ),
            # don't trim 2 mismatches in middle
            TrimR2TestCase(
                name="NoTrim:2 mismatches in middle",
                bases2=keep_seq
                + sequence.reverse_complement(adapter_seq[:4] + "CC" + adapter_seq[6:]),
                aln_score=-1,
                trimmed_len=31,
                keep=True,
                trim=False,
            ),
            # don't trim 2 mismatches at end
            TrimR2TestCase(
                name="NoTrim:2 mismatches at end",
                bases2=keep_seq + sequence.reverse_complement(adapter_seq[:-2] + "CC"),
                aln_score=-1,
                trimmed_len=31,
                keep=True,
                trim=False,
            ),
            # don't trim inserted base at beginning
            TrimR2TestCase(
                name="Trim:inserted base at beginning",
                bases2=keep_seq + sequence.reverse_complement("C" + adapter_seq),
                aln_score=19,
                trimmed_len=20,
                keep=True,
                trim=True,
            ),
            # don't trim deleted base at beginning
            TrimR2TestCase(
                name="Trim:deleted base at beginning",
                bases2=keep_seq + sequence.reverse_complement(adapter_seq[1:]),
                aln_score=18,
                trimmed_len=18,
                keep=True,
                trim=True,
            ),
            # don't trim reverse complement of mosaic end
            TrimR2TestCase(
                name="NoTrim:reverse complement of mosaic end",
                bases2=keep_seq + adapter_seq,
                aln_score=-1,
                trimmed_len=31,
                keep=True,
                trim=False,
            ),
            # don't trim read2 with tn5 mosaic end (ignore, and reject because read1 doesn't match)
            TrimR2TestCase(
                name="NoTrim:read1 starts with mosaic end",
                bases1=adapter_seq + keep_seq,
                bases2=keep_seq
                + sequence.reverse_complement(adapter_seq[:4] + "CC" + adapter_seq[6:]),
                aln_score=-1,
                trimmed_len=31,
                keep=True,
                trim=False,
            ),
            # don't trim read2 with reverse compliment of mosaic end (ignore, and reject because
            # read1 doesn't match)
            TrimR2TestCase(
                name="NoTrim:read2 starts with reverse-complement mosaic end",
                bases1=sequence.reverse_complement(adapter_seq) + keep_seq,
                bases2=keep_seq
                + sequence.reverse_complement(adapter_seq[:4] + "CC" + adapter_seq[6:]),
                aln_score=-1,
                trimmed_len=31,
                keep=True,
                trim=False,
            ),
            # reject read2 that is too short after trimming
            TrimR2TestCase(
                name="Reject:read2 too short after trimming",
                bases2=keep_seq[:-1] + sequence.reverse_complement(adapter_seq),
                aln_score=-1,
                trimmed_len=-1,
                keep=False,
                trim=False,
            ),
        ],
    )

    actual_mismatch_metrics = trim_r2(
        in_bam=test_builder.in_bam,
        keep_bam=test_builder.keep_bam,
        reject_bam=test_builder.reject_bam,
        out_metrics=test_builder.metrics_tsv,
        sequence=adapter_seq,
        min_score=min_score,
        min_match_length=min_match_length,
        min_keep_length=SEQ_LEN,
    )

    expected_mismatch_metrics = test_builder.expected_mismatch_metric
    if expected_mismatch_metrics != actual_mismatch_metrics:
        # to make the error easier to track down, decrease the size of the assert error by removing
        # common stuff first
        for num_mismatches, expected_metric in list(expected_mismatch_metrics.items()):
            if (
                num_mismatches in actual_mismatch_metrics
                and actual_mismatch_metrics[num_mismatches] == expected_metric
            ):
                expected_mismatch_metrics.pop(num_mismatches)
                actual_mismatch_metrics.pop(num_mismatches)
        assert expected_mismatch_metrics == actual_mismatch_metrics

    expected_test_cases = sorted(test_builder.test_cases, key=lambda test_case: test_case.name)
    actual_test_cases = sorted(
        TrimR2TestCase.from_bams(
            keep_bam=test_builder.keep_bam,
            reject_bam=test_builder.reject_bam,
            adapter_len=len(test_builder.adapter_seq),
        ),
        key=lambda test_case: test_case.name,
    )
    assert len(expected_test_cases) == len(actual_test_cases)
    for expected_test_case, actual_test_case in zip(expected_test_cases, actual_test_cases):
        TrimR2TestCase.assert_equal(expected=expected_test_case, actual=actual_test_case)
