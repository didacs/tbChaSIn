import logging
from Bio import Align
from collections import Counter
from pathlib import Path
from pysam import AlignmentFile
from samwell import sam
from samwell.dnautils import reverse_complement


def setup_aligner() -> Align.PairwiseAligner:
    """Creates a glocal aligner (partial query, full target) using BWA mem inspired scoring
    parameters"""
    # Gapped alignment parameters
    # NB: setting zero for the target_end_* options effectively makes this glocal, allowing gaps
    # after the end of the target (can open gaps at the end of the query)
    aligner: Align.PairwiseAligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 1
    aligner.mismatch_score = -4
    aligner.open_gap_score = -7
    aligner.extend_gap_score = -1
    aligner.query_internal_open_gap_score = -7
    aligner.query_internal_extend_gap_score = -1
    aligner.query_end_open_gap_score = 0
    aligner.query_end_extend_gap_score = 0
    aligner.target_internal_open_gap_score = -7
    aligner.target_internal_extend_gap_score = -1
    aligner.target_end_open_gap_score = -7
    aligner.target_end_extend_gap_score = -1
    return aligner


def trim_for_tn5me(
    *,
    in_bam: Path,
    out_bam: Path,
    out_metrics: Path,
    circularization_palindrome: str = "ACGT",
    tn5_mosaic_end: str = "AGATGTGTATAAGAGACAG",
    min_score: int = 19,
) -> None:
    """Finds the Tn5 mosaic end, trims it and subsequent bases.

    In some cases, the reads will sequence into (and perhaps across) the circularization
    breakpoint.  This tool will attempt to find evidence for this, and trim the sequence in the
    breakpoint and any subsequent bases.

    The `t5` SAM tag will store the number of bases trimmed from the 3' end of the read and
    alignment score, semicolon delimited.

    Args:
        in_bam: the input unmapped BAM
        out_bam: the output trimmed BAM
        out_metrics: the output trimming metrics
        circularization_palindrome: the circular palindrome used to circularize by ligation
        tn5_mosaic_end: the Tn5 mosaic end used for tagmentation
        min_score: the minimum alignment score for trimming
    """
    logger = logging.getLogger(__name__)

    min_leading_to_keep = 12

    # assert min_score >= 0 since we count unmatched reads as -1
    assert min_score >= 0, "--min-score must be >= 0"
    assert circularization_palindrome == reverse_complement(
        circularization_palindrome
    ), f"Not palindromic: {circularization_palindrome}"

    # The full sequence we're searching for.
    full_target: str = (
        reverse_complement(tn5_mosaic_end) + circularization_palindrome + tn5_mosaic_end
    )

    aligner: Align.PairwiseAligner = setup_aligner()

    counter = Counter({-1: 0})
    with AlignmentFile(str(in_bam), check_sq=False) as reader, sam.writer(
        out_bam, header=reader.header
    ) as writer:
        record_number: int = 1
        for record in reader:
            assert record.is_unmapped, f"Record was mapped: {record.name}"

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

    # write out the metrics
    with out_metrics.open("w") as writer:
        total = sum(1 for _ in counter.elements())
        running_sum = 0
        writer.write("alignment_score\tcount\tfrac_at_score\tfrac_ge_score\n")
        for score, count in sorted(counter.items(), key=lambda tup: -tup[0]):
            running_sum += count
            frac_at_score = float(count) / total
            frac_ge_score = float(running_sum) / total
            writer.write(f"{score}\t{count:,d}\t{frac_at_score:.4f}\t{frac_ge_score:.4f}\n")
