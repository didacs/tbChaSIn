import logging
from collections import Counter
from pathlib import Path
from typing import Dict
from typing import Optional

from attr import frozen
from Bio import Align
from fgpyo import sam
from fgpyo.sam import Template
from fgpyo.sequence import reverse_complement
from fgpyo.util.metric import Metric
from pysam import AlignmentFile

from pytomebio.core.aligner import get_end_adapter_aligner
from pytomebio.tools.common import TN5_ME


@frozen
class ScoreCountsMetric(Metric["ScoreCountsMetric"]):
    """Stores counts for how reads matched the sequences.

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


def trim_r2(
    *,
    in_bam: Path,
    keep_bam: Path,
    reject_bam: Path,
    out_metrics: Path,
    sequence: str = TN5_ME,
    min_score: Optional[int] = None,
    min_match_length: Optional[int] = None,
    min_keep_length: Optional[int] = None,
) -> Dict[int, ScoreCountsMetric]:
    """Finds the specified sequence (the Tn5 mosaic end by default) in R2 and trims it and
    subsequent bases. It is strongly recommended to specify at least one of `--min-score` or
    `--min-match-length`.

    For reads that are trimmed, the following SAM tags will be added:
    * `t3`: the trimmed sequence
    * `tm`: the number of mismatches between the target sequence and the read
    * `ts`: the score of the best alignment between the target sequence and the read

    Args:
        in_bam: the input unmapped BAM
        out_bam: the output trimmed BAM
        out_metrics: the output trimming metrics TSV
        sequence: the sequence to trim; defaults to the Tn5 mosaic end
        anchored: whether the sequence is anchored at the end of the read (True) or gaps are
        allowed at the beginning of the query or target (False)
        min_score: optional; the minimum alignment score for trimming; by default, the read is
        trimmed if the sequence aligns
        min_match_length: optional; the minimum number of bases that must match between the target
        sequence and the read (alignment length - mismatches) for trimming; by default, the read is
        trimmed if the sequence aligns
        min_keep_length: optional; the minimum number of remaining R2 bases after trimming required
        to keep the read pair; by default, the trimmed read is always kept
    """
    logger = logging.getLogger(__name__)

    # assert min_score >= 0 since we count unmatched reads as -1
    if min_score is not None:
        assert min_score >= 0, "--min-score must be >= 0"
    if min_match_length is not None:
        assert min_match_length >= 0, "--min-match-length must be >= 0"
    if min_keep_length is not None:
        assert min_keep_length >= 0, "--min-keep-length must be >= 0"

    aligner: Align.PairwiseAligner = get_end_adapter_aligner()
    # The target sequence we're searching for - reverse complement since we're matching to R2
    target: str = reverse_complement(sequence)
    # count the number of times we observe each alignment score
    score_counter = Counter({-1: 0})
    # total number of records processed
    record_number: int = 0

    with (
        AlignmentFile(str(in_bam), check_sq=False) as reader,
        sam.writer(keep_bam, header=reader.header) as keep_writer,
        sam.writer(reject_bam, header=reader.header) as reject_writer,
    ):
        for template in Template.iterator(reader):
            assert template.r1 is not None, f"Missing r1 for {template.name}"
            assert template.r2 is not None, f"Missing r2 for {template.name}"
            template.validate()

            keep = True
            count_score = -1

            alignments = aligner.align(template.r2.query_sequence, target)
            if len(alignments) > 0 and len(alignments[0].aligned[0]) > 0:
                aln = alignments[0]
                # print(aln)
                score = round(aln.score)
                start, end = aln.aligned[0][0]
                mismatches = aln.counts().mismatches
                # print(f"{start}:{end} {score} {mismatches}")
                trim = True
                if min_score is not None:
                    # only trim the read if the alignment score is at least the min score
                    trim = score >= min_score
                if trim and min_match_length is not None:
                    # only trim the read if the at least the min number of sequence bases match
                    trim = end - start - mismatches >= min_match_length
                if trim and min_keep_length is not None:
                    # only keep the read pair if at least the min number of R2 bases will be kept
                    trim = keep = start >= min_keep_length
                if trim:
                    count_score = score
                    template.r2.set_tag("t3", template.r2.query_sequence[start:])
                    template.r2.set_tag("ts", score)
                    template.r2.set_tag("tm", mismatches)
                    query_qualities = template.r2.query_qualities
                    template.r2.query_sequence = template.r2.query_sequence[:start]
                    template.r2.query_qualities = query_qualities[:start]

            score_counter[count_score] += 1

            for rec in template.all_recs():
                record_number += 1
                if record_number % 10000 == 0:
                    logger.info(f"Processed {record_number:,d} records")
                if keep:
                    keep_writer.write(rec)
                else:
                    reject_writer.write(rec)

    logger.info(f"Processed {record_number:,d} records")

    # convert raw Counter into formatted ScoreCountsMetric
    score_counts_metrics = ScoreCountsMetric.from_counter(score_counter)
    # write TSV. Note that ScoreCountsMetrics are inserted in order of increasing alignment score
    # (and dictionary view objects preserve insertion order), so there is no need to sort here
    ScoreCountsMetric.write(out_metrics, *score_counts_metrics.values())
    return score_counts_metrics
