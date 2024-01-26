from attr import frozen
import logging
from collections import Counter
from fgpyo.sam import Template
from fgpyo.sam import TemplateIterator
from fgpyo.util.metric import Metric
from pathlib import Path
from pysam import AlignmentFile
from samwell import sam
from typing import Dict


def _num_mismatches(s1: str, s2: str) -> int:
    """Returns the number of mismatches between two sequences"""
    mm: int = 0
    for b1, b2 in zip(s1, s2):
        if b1 != b2:
            mm += 1
    return mm


@frozen
class Tn5ReadOneMatchingMetric(Metric["Tn5ReadOneMatchingMetric"]):
    """Stores counts for records with specified number of mismatches between R1 and the Tn5
    transposase mosaic end

    Attributes:
        num_mismatches: The number of mismatches for an alignment between read 1 and Tn5 mosaic end
                        for one or more SAM records
        count: the number of SAM records that were aligned with this number of mismatches
        frac_at_num_mismatches: the proportion of SAM records with this many mismatches
        frac_ge_num_mismatches: the proportion of SAM records with this many mismatches or greater
        frac_le_num_mismatches: the proportion of SAM records with this many mismatches or fewer
        kept: true if this number of mismatches was kept (<= max_mismatches), false if rejected
    """

    num_mismatches: int
    count: int
    frac_at_num_mismatches: float
    frac_le_num_mismatches: float
    frac_ge_num_mismatches: float
    kept: bool

    @staticmethod
    def from_counter(
        counter: Counter[int],
        max_mismatches: int,
    ) -> Dict[int, "Tn5ReadOneMatchingMetric"]:
        """Return a dict from number of mismatches to corresponding Tn5ReadOneMatchingMetric

        Args:
            counter: Counter that stores number of occurrences for a given number of mismatches
            max_mismatches: the maximum number of mismatches allowed
        """
        metrics: Dict[int, "Tn5ReadOneMatchingMetric"] = {}
        total = sum(1 for _ in counter.elements())
        running_sum_up = 0
        running_sum_down = total
        # ensure every number of mismatches is present in range [0, max_observed_mismatches + 1]
        for num_mismatches in range(0, max(counter.keys()) + 1):
            count = counter.get(num_mismatches, 0)
            running_sum_up += count
            metrics[num_mismatches] = Tn5ReadOneMatchingMetric(
                num_mismatches=num_mismatches,
                count=count,
                frac_at_num_mismatches=count / total if total > 0 else 0.0,
                frac_le_num_mismatches=running_sum_up / total if total > 0 else 0.0,
                frac_ge_num_mismatches=running_sum_down / total if total > 0 else 0.0,
                kept=num_mismatches <= max_mismatches,
            )
            running_sum_down -= count
        return metrics


def trim_for_tn5me(
    *,
    in_bam: Path,
    keep_bam: Path,
    reject_bam: Path,
    out_metrics: Path,
    tn5_mosaic_end: str = "AGATGTGTATAAGAGACAG",
    max_mismatches: int = 1,
) -> Dict[int, Tn5ReadOneMatchingMetric]:
    """Finds the Tn5 mosaic end at the start of R1 and trims it.

    Any read pair that where the Tn5 mosaic end is not found at the start of R1 is rejected.

    The `t5` SAM tag will store the raw bases trimmed, while the `tm` SAM tag will store the number
    of mismatches in the alignment.

    Args:
        in_bam: the input unmapped BAM
        keep_bam: the output BAM with reads with the Tn5 mosaic end trimmed in R1
        reject_bam: the output BAM with reads where the leading Tn5 mosaic end in R1 was not found
        out_metrics: the output trimming metrics
        tn5_mosaic_end: the Tn5 mosaic end used for tagmentation
        max_mismatches: the maximum number of mismatches to allow
    """
    logger = logging.getLogger(__name__)

    tn5_mosaic_end = tn5_mosaic_end.upper()
    tn5_mosaic_end_len = len(tn5_mosaic_end)

    # assert min_score >= 0 since we count unmatched reads as -1
    assert max_mismatches >= 0, "--max-mismatches must be >= 0"

    counter = Counter({0: 0})
    with AlignmentFile(str(in_bam), check_sq=False) as reader, sam.writer(
        keep_bam, header=reader.header
    ) as keep_writer, sam.writer(reject_bam, header=reader.header) as reject_writer:
        record_number: int = 1

        template_iter: TemplateIterator = Template.iterator(reader)

        for template in template_iter:
            assert template.r1 is not None, f"Missing r1 for {template.name}"
            assert template.r2 is not None, f"Missing r2 for {template.name}"
            template.validate()

            r1_bases = template.r1.query_sequence
            num_mismatches = _num_mismatches(tn5_mosaic_end, r1_bases)
            counter[num_mismatches] += 1
            if num_mismatches <= max_mismatches:
                template.r1.set_tag("t5", r1_bases[:tn5_mosaic_end_len])
                template.r1.set_tag("tm", num_mismatches)
                r1_quals = template.r1.query_qualities
                template.r1.query_sequence = r1_bases[tn5_mosaic_end_len:]
                template.r1.query_qualities = r1_quals[tn5_mosaic_end_len:]
                for rec in template.all_recs():
                    keep_writer.write(rec)
            else:
                for rec in template.all_recs():
                    reject_writer.write(rec)

            for _ in template.all_recs():
                if record_number % 100000 == 0:
                    logger.info(f"Processed {record_number:,d} records")
                record_number += 1
        logger.info(f"Processed {record_number:,d} records")

    # convert raw Counter into formatted MismatchMetric
    mismatch_metrics = Tn5ReadOneMatchingMetric.from_counter(
        counter=counter, max_mismatches=max_mismatches
    )
    # write TSV. Note that Tn5ReadOneMatchingMetrics are inserted in order of increasing number
    # of mismatches (and dictionary view objects preserve insertion order), so there is no need to
    # sort in this step
    Tn5ReadOneMatchingMetric.write(out_metrics, *mismatch_metrics.values())
    return mismatch_metrics
