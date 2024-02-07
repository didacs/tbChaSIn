import collections
import typing
from pathlib import Path
from typing import Dict
from typing import Iterator
from typing import List
from typing import Optional

from attr import frozen
from fgpyo import sequence
from fgpyo.sam import Template
from fgpyo.sam.builder import SamBuilder
from pysam import AlignmentFile

from pytomebio.tools.cryptic_seq.trim_for_tn5me import Tn5ReadOneMatchingMetric
from pytomebio.tools.cryptic_seq.trim_for_tn5me import trim_for_tn5me


@frozen
class TrimForTn5meTestCase:
    """
    Holds an individual test case and expected/actual outcome from the trim_for_tn5me tool
    Also has code to:
        - reconstitute TrimForTn5meTestCase objects from the actual outcomes in the output bam
          files
        - assert two TrimForTn5meTestCase are equal

    Attributes:
        name: description of test case (no whitespace other than spaces)
        bases1: bases on read1
        num_mismatches: number of mismatches when aligned to Tn5 mosaic end
        keep: whether the read pair is accepted or rejected
        bases2: bases on read2. Should not be examined by the tool so defaults to a simple default
                string. In a few test cases we will override this to test that bases2 are ignored
                even if they match the Tn5 mosaic end
    """

    name: str
    bases1: str
    num_mismatches: int
    keep: bool
    bases2: str = "GATCGATCGATC"

    @property
    def sam_name(self) -> str:
        """return a SAM-compatible name"""
        return self.name.replace(" ", "_")

    @staticmethod
    def assert_equal(expected: "TrimForTn5meTestCase", actual: "TrimForTn5meTestCase") -> None:
        """To be equivalent, everything must be equal, except that num-mismatches may differ if
        the test case is not accepted (because the mismatches tag isn't added to the reject bam,
        so that information is lost)
        """
        assert expected.name == actual.name
        assert expected.bases1 == actual.bases1, f"'{expected.name}' has wrong bases1"
        assert expected.keep == actual.keep, f"'{expected.name}' has wrong keep"
        if expected.keep:
            # if the test case is not accepted, it's not possible to reconstruct the number of
            # mismatches. However, this is still tested (indirectly) by comparing the metrics
            assert (
                expected.num_mismatches == actual.num_mismatches
            ), f"'{expected.name}' has wrong num_mismatches"
        assert expected.bases2 == actual.bases2, f"'{expected.name}' has wrong bases2"

    @staticmethod
    def from_bams(
        keep_bam: Path,
        reject_bam: Path,
        tn5_mosaic_end_len: int,
    ) -> List["TrimForTn5meTestCase"]:
        """Reconstitute test cases from keep and reject BAMs output by trim_for_tn5me"""
        return [
            TrimForTn5meTestCase._from_template(
                keep=True, template=template, tn5_mosaic_end_len=tn5_mosaic_end_len
            )
            for template in TrimForTn5meTestCase._bam_to_template_iterator(keep_bam)
        ] + [
            TrimForTn5meTestCase._from_template(keep=False, template=template)
            for template in TrimForTn5meTestCase._bam_to_template_iterator(reject_bam)
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
        tn5_mosaic_end_len: Optional[int] = None,
    ) -> "TrimForTn5meTestCase":
        """Reconstruct a  TrimForTn5meTestCase from read pair, and validate tag correctness"""
        template.validate()
        human_readable_name = f"{template.name}".replace("_", " ")
        assert template.r1 is not None
        assert template.r1.is_read1
        assert template.r2 is not None
        assert template.r2.is_read2
        # read 2 should never be tagged
        assert not template.r2.has_tag("t5")
        assert not template.r2.has_tag("tm")
        if keep:
            # read 1 is tagged and trimmed if it is accepted
            assert template.r1.has_tag("t5")
            assert template.r1.has_tag("tm")
            t5_tag = typing.cast(str, template.r1.get_tag("t5"))
            assert len(t5_tag) == tn5_mosaic_end_len
            return TrimForTn5meTestCase(
                name=human_readable_name,
                bases1=t5_tag + template.r1.query_sequence,
                bases2=template.r2.query_sequence,
                num_mismatches=int(template.r1.get_tag("tm")),
                keep=True,
            )
        else:
            # read 1 is neither tagged nor trimmed if it is accepted
            assert not template.r1.has_tag("t5")
            assert not template.r1.has_tag("tm")
            return TrimForTn5meTestCase(
                name=human_readable_name,
                bases1=template.r1.query_sequence,
                bases2=template.r2.query_sequence,
                num_mismatches=-1,
                keep=False,
            )


@frozen
class TrimForTn5meTestBuilder:
    """
    Container for test cases.
    - hold test data raw inputs (tn5_mosaic_end, max_mismatches, test_cases)
    - write input bam file
    - set path for output bam files and metric TSV file
    - generate the expected metric based on input test cases
    """

    tn5_mosaic_end: str
    max_mismatches: int
    test_cases: List[TrimForTn5meTestCase]
    in_bam: Path
    keep_bam: Path
    reject_bam: Path
    metrics_tsv: Path

    @staticmethod
    def from_test_cases(
        test_dir: Path,
        tn5_mosaic_end: str,
        max_mismatches: int,
        test_cases: List[TrimForTn5meTestCase],
    ) -> "TrimForTn5meTestBuilder":
        """Generate and store data for end-to-end test of cryptic_seq/trim_for_tn5me inside
        test_dir

        Args:
            test_dir: path to directory that should contain all the temporary files for test
            tn5_mosaic_end: the Tn5 mosaic end used for tagmentation in the test
            max_mismatches: the maximum number of mismatches to allow for passing test cases
            test_cases: list of TrimForTn5meTestCase that will comprehensively test logic flow of
                        trim_for_tn5me tool
        """
        TrimForTn5meTestBuilder.assert_names_unique(test_cases)

        in_bam = test_dir / "in.bam"
        input_bam_builder = SamBuilder()
        for test_case in test_cases:
            # bases2 are never examined, so generate them randomly
            input_bam_builder.add_pair(
                name=test_case.sam_name, bases1=test_case.bases1, bases2=test_case.bases2
            )
        input_bam_builder.to_path(path=in_bam)
        return TrimForTn5meTestBuilder(
            tn5_mosaic_end=tn5_mosaic_end,
            max_mismatches=max_mismatches,
            test_cases=test_cases,
            in_bam=in_bam,
            keep_bam=test_dir / "keep.bam",
            reject_bam=test_dir / "reject.bam",
            metrics_tsv=test_dir / "mismatch_metrics.tsv",
        )

    @staticmethod
    def assert_names_unique(test_cases: List[TrimForTn5meTestCase]) -> None:
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
    def expected_mismatch_metric(self) -> Dict[int, Tn5ReadOneMatchingMetric]:
        """Generate expected Dict of Tn5ReadOneMatchingMetric based on input test_cases"""
        counter = collections.Counter({0: 0})
        for test_case in self.test_cases:
            counter[test_case.num_mismatches] += 1
        return Tn5ReadOneMatchingMetric.from_counter(
            counter=counter, max_mismatches=self.max_mismatches
        )


def test_trim_for_tn5me(
    tmp_path: Path,
    tn5_mosaic_end: str = "AGATGTGTATAAGAGACAG",
    max_mismatches: int = 1,
) -> None:
    test_builder = TrimForTn5meTestBuilder.from_test_cases(
        tn5_mosaic_end=tn5_mosaic_end,
        max_mismatches=max_mismatches,
        test_dir=tmp_path,
        test_cases=[
            # accept 0 mismatches
            TrimForTn5meTestCase(
                name="Accept:perfect match",
                bases1=tn5_mosaic_end + "AAAAAAAA",
                num_mismatches=0,
                keep=True,
            ),
            # accept 1 mismatch at beginning
            TrimForTn5meTestCase(
                name="Accept:1 mismatch at beginning",
                bases1="C" + tn5_mosaic_end[1:] + "AAAAAAAA",
                num_mismatches=1,
                keep=True,
            ),
            # accept 1 mismatch in middle
            TrimForTn5meTestCase(
                name="Accept:1 mismatch in middle",
                bases1=tn5_mosaic_end[:4] + "C" + tn5_mosaic_end[5:] + "AAAAAAAA",
                num_mismatches=1,
                keep=True,
            ),
            # accept 1 mismatch at end
            TrimForTn5meTestCase(
                name="Accept:1 mismatch at end",
                bases1=tn5_mosaic_end[:-1] + "C" + "AAAAAAAA",
                num_mismatches=1,
                keep=True,
            ),
            # reject 2 mismatches at beginning
            TrimForTn5meTestCase(
                name="Reject:2 mismatches at beginning",
                bases1="CC" + tn5_mosaic_end[2:] + "AAAAAAAA",
                num_mismatches=2,
                keep=False,
            ),
            # reject 2 mismatches in middle
            TrimForTn5meTestCase(
                name="Reject:2 mismatches in middle",
                bases1=tn5_mosaic_end[:4] + "CC" + tn5_mosaic_end[6:] + "AAAAAAAA",
                num_mismatches=2,
                keep=False,
            ),
            # reject 2 mismatches at end
            TrimForTn5meTestCase(
                name="Reject:2 mismatches at end",
                bases1=tn5_mosaic_end[:-2] + "CC" + "AAAAAAAA",
                num_mismatches=2,
                keep=False,
            ),
            # reject inserted base at beginning
            TrimForTn5meTestCase(
                name="Reject:inserted base at beginning",
                bases1="C" + tn5_mosaic_end + "AAAAAAAA",
                num_mismatches=18,
                keep=False,
            ),
            # reject deleted base at beginning
            TrimForTn5meTestCase(
                name="Reject:deleted base at beginning",
                bases1=tn5_mosaic_end[1:] + "AAAAAAAA",
                num_mismatches=18,
                keep=False,
            ),
            # reject reverse complement of mosaic end
            TrimForTn5meTestCase(
                name="Reject:reverse complement of mosaic end",
                bases1=sequence.reverse_complement(tn5_mosaic_end) + "AAAAAAAA",
                num_mismatches=13,
                keep=False,
            ),
            # reject read2 with tn5 mosaic end (ignore, and reject because read1 doesn't match)
            TrimForTn5meTestCase(
                name="Reject:read2 starts with mosaic end",
                bases1=tn5_mosaic_end[:-2] + "CC" + "AAAAAAAA",
                bases2=tn5_mosaic_end + "AAAAAAAA",
                num_mismatches=2,
                keep=False,
            ),
            # reject read2 with reverse compliment of mosaic end (ignore, and reject because read1
            # doesn't match)
            TrimForTn5meTestCase(
                name="Reject:read2 starts with reverse-complement mosaic end",
                bases1=tn5_mosaic_end[:-2] + "CC" + "AAAAAAAA",
                bases2=sequence.reverse_complement(tn5_mosaic_end) + "AAAAAAAA",
                num_mismatches=2,
                keep=False,
            ),
        ],
    )

    actual_mismatch_metrics = trim_for_tn5me(
        in_bam=test_builder.in_bam,
        keep_bam=test_builder.keep_bam,
        reject_bam=test_builder.reject_bam,
        out_metrics=test_builder.metrics_tsv,
        tn5_mosaic_end=tn5_mosaic_end,
        max_mismatches=max_mismatches,
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
        TrimForTn5meTestCase.from_bams(
            keep_bam=test_builder.keep_bam,
            reject_bam=test_builder.reject_bam,
            tn5_mosaic_end_len=len(test_builder.tn5_mosaic_end),
        ),
        key=lambda test_case: test_case.name,
    )
    assert len(expected_test_cases) == len(actual_test_cases)
    for expected_test_case, actual_test_case in zip(expected_test_cases, actual_test_cases):
        TrimForTn5meTestCase.assert_equal(expected=expected_test_case, actual=actual_test_case)
