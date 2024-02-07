import logging
from collections import Counter
from pathlib import Path
from typing import Dict

from attr import frozen
from Bio import Align
from fgpyo import sam
from fgpyo.sequence import reverse_complement
from fgpyo.util.metric import Metric
from pysam import AlignmentFile

from pytomebio.core.aligner import get_glocal_aligner


@frozen
class ScoreCountsMetric(Metric["ScoreCountsMetric"]):
    """Stores counts for how reads matched the attachment sites.

    Attributes:
        alignment_score: An alignment score observed for one or more SAM records
        count: the number of SAM records that were aligned with this score
        frac_at_score: the proportion of SAM records that aligned with this score
        frac_ge_score: the proportion of SAM records with this score or better
    """

    alignment_score: int
    count: int
    frac_at_score: float
    frac_ge_score: float

    @staticmethod
    def from_counter(counter: Counter[int]) -> Dict[int, "ScoreCountsMetric"]:
        """Return a dict from alignment score to corresponding ScoreCountsMetric

        Args:
            counter: Counter that stores number of occurrences for a given alignment score
        """
        metrics: Dict[int, ScoreCountsMetric] = {}
        total = sum(1 for _ in counter.elements())
        running_sum = 0
        for score, count in sorted(counter.items(), key=lambda tup: -tup[0]):
            running_sum += count
            metrics[score] = ScoreCountsMetric(
                alignment_score=score,
                count=count,
                frac_at_score=count / total if total > 0 else 0.0,
                frac_ge_score=running_sum / total if total > 0 else 0.0,
            )
        return metrics


def trim_for_tn5me(
    *,
    in_bam: Path,
    out_bam: Path,
    out_metrics: Path,
    circularization_palindrome: str = "ACGT",
    tn5_mosaic_end: str = "AGATGTGTATAAGAGACAG",
    min_score: int = 19,
    is_circularized: bool = True,
) -> Dict[int, ScoreCountsMetric]:
    """Finds the Tn5 mosaic end, trims it and subsequent bases.

    In some cases, the reads will sequence into (and perhaps across) the circularization
    breakpoint.  This tool will attempt to find evidence for this, and trim the sequence in the
    breakpoint and any subsequent bases.

    The `t5` SAM tag will store the number of bases trimmed from the 3' end of the read and
    alignment score, semicolon delimited.

    Args:
        in_bam: the input unmapped BAM
        out_bam: the output trimmed BAM
        out_metrics: the output trimming metrics TSV
        circularization_palindrome: the circular palindrome used to circularize by ligation
        tn5_mosaic_end: the Tn5 mosaic end used for tagmentation
        min_score: the minimum alignment score for trimming
        is_circularized: true to search for the full circularized sequence (revcomp(Tn5ME) +
                         palindrome + Tn5ME), for example for Change-Seq, false to just search for
                         revcomp(Tn5ME), for example for Cryptic-Seq.
    """
    logger = logging.getLogger(__name__)

    min_leading_to_keep = 12

    # assert min_score >= 0 since we count unmatched reads as -1
    assert min_score >= 0, "--min-score must be >= 0"
    assert circularization_palindrome == reverse_complement(
        circularization_palindrome
    ), f"Not palindromic: {circularization_palindrome}"

    # The full sequence we're searching for.
    full_target: str
    if is_circularized:
        full_target = (
            reverse_complement(tn5_mosaic_end) + circularization_palindrome + tn5_mosaic_end
        )
    else:
        full_target = reverse_complement(tn5_mosaic_end)

    aligner: Align.PairwiseAligner = get_glocal_aligner()

    counter = Counter({-1: 0})
    with AlignmentFile(str(in_bam), check_sq=False) as reader, sam.writer(
        out_bam, header=reader.header
    ) as writer:
        record_number: int = 1
        for record in reader:
            assert record.is_unmapped, f"Record was mapped: {record.query_name}"

            alignments = aligner.align(record.query_sequence, full_target)
            if len(alignments) == 0 or alignments[0].score < min_score:
                counter[-1] += 1
            else:
                counter[alignments[0].score] += 1
                # the first aligned base in the read/query is the number of bases we should keep
                num_leading_to_keep = alignments[0].aligned[0][0][0]
                if min_leading_to_keep <= num_leading_to_keep:
                    record.set_tag(
                        "t5", f"{record.query_length - num_leading_to_keep};{alignments[0].score}"
                    )
                    query_qualities = record.query_qualities
                    record.query_sequence = record.query_sequence[:num_leading_to_keep]
                    record.query_qualities = query_qualities[:num_leading_to_keep]
            writer.write(record)

            if record_number % 10000 == 0:
                logger.info(f"Processed {record_number:,d} records")
            record_number += 1
        logger.info(f"Processed {record_number:,d} records")

    # convert raw Counter into formatted ScoreCountsMetric
    score_counts_metrics = ScoreCountsMetric.from_counter(counter)
    # write TSV. Note that ScoreCountsMetrics are inserted in order of increasing alignment score
    # (and dictionary view objects preserve insertion order), so there is no need to sort here
    ScoreCountsMetric.write(out_metrics, *score_counts_metrics.values())
    return score_counts_metrics
