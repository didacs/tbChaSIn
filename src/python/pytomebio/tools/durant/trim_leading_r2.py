import gzip
import logging
from collections import Counter
from pathlib import Path

from pysam import FastxFile


def _num_mismatches(s1: str, s2: str) -> int:
    """Returns the number of mismatches between two sequences"""
    mm: int = 0
    for b1, b2 in zip(s1, s2):
        if b1 != b2:
            mm += 1
    return mm


def trim_leading_r2(
    *,
    in_r1_fq: Path,
    in_r2_fq: Path,
    out_prefix: Path,
    stagger: str,
    donor_inner_primer: str,
    max_mismatches: int = 1,
) -> None:
    """Finds R2s with the leading stagger and donor inner primer and trims it.

    Any read pair where the leading stagger and donor inner primer are not found at the start of R2
    are rejected.

    The following output files will be written:
    - <output>.keep.R1.fq.gz: the R1s kept
    - <output>.keep.R2.fq.gz: the R1s kept
    - <output>.reject.R1.fq.gz: the R1s rejected
    - <output>.reject.R2.fq.gz: the R1s rejected
    - <output>.metrics.txt: the output metrics file

    Args:
        in_r1_fq: the input R1 FASTQ
        in_r2_fq: the input R1 FASTQ
        out_prefix: the output prefix
        stagger: the stagger sequence
        donor_inner_primer: the R2 donor inner primer
        max_mismatches: the maximum number of mismatches to allow
    """
    logger = logging.getLogger(__name__)

    # assert min_score >= 0 since we count unmatched reads as -1
    assert max_mismatches >= 0, "--max-mismatches must be >= 0"
    assert len(donor_inner_primer) > 0, "donor inner primer must be given"

    prefix_seq: str = (stagger + donor_inner_primer).upper()
    prefix_seq_len: int = len(prefix_seq)

    counter = Counter({0: 0})

    out_keep_r1_fq = f"{out_prefix}.keep.R1.fq.gz"
    out_keep_r2_fq = f"{out_prefix}.keep.R2.fq.gz"
    out_reject_r1_fq = f"{out_prefix}.reject.R1.fq.gz"
    out_reject_r2_fq = f"{out_prefix}.reject.R2.fq.gz"
    out_metrics_txt = Path(f"{out_prefix}.metrics.txt")

    with (
        # open read1 and read2 fastq files for reading
        FastxFile(str(in_r1_fq)) as in_r1_fh,
        FastxFile(str(in_r2_fq)) as in_r2_fh,
        # open read1 and read2 fastq files for writing
        gzip.open(str(out_keep_r1_fq), mode="wt") as out_keep_r1_fh,
        gzip.open(str(out_keep_r2_fq), mode="wt") as out_keep_r2_fh,
        gzip.open(str(out_reject_r1_fq), mode="wt") as out_reject_r1_fh,
        gzip.open(str(out_reject_r2_fq), mode="wt") as out_reject_r2_fh,
    ):
        record_number: int = 1

        for r1, r2 in zip(in_r1_fh, in_r2_fh):
            assert r1.name == r2.name, (
                f"Error: FASTQ names disagreed for the {record_number}th entry"
                f" '{r1.name}' != '{r2.name}'"
            )

            r2_bases = r2.sequence
            num_mismatches = _num_mismatches(prefix_seq, r2_bases.upper())
            counter[num_mismatches] += 1

            if num_mismatches <= max_mismatches:
                r2.set_sequence(
                    sequence=r2_bases[prefix_seq_len:], quality=r2.quality[prefix_seq_len:]
                )
                out_keep_r1_fh.write(f"{r1}\n")
                out_keep_r2_fh.write(f"{r2}\n")
            else:
                out_reject_r1_fh.write(f"{r1}\n")
                out_reject_r2_fh.write(f"{r2}\n")

            for _ in (r1, r2):
                if record_number % 100000 == 0:
                    logger.info(f"Processed {record_number:,d} records")
                record_number += 1
        logger.info(f"Processed {record_number:,d} records")

    for i in range(max(counter.elements()) + 1):
        counter[i] += 0

    # write out the metrics
    with out_metrics_txt.open("w") as writer:
        total = sum(1 for _ in counter.elements())
        running_sum = 0
        writer.write(
            "\t".join(
                [
                    "num_mismatches",
                    "count",
                    "frac_at_num_mismatches",
                    "frac_le_num_mismatches",
                    "frac_ge_num_mismatches",
                    "kept",
                ]
            )
            + "\n"
        )
        for num_mismatches, count in sorted(counter.items(), key=lambda tup: tup[0]):
            running_sum += count
            frac_at_num_mismatches = float(count) / total
            frac_ge_num_mismatches = float(total - running_sum) / total
            frac_le_num_mismatches = float(running_sum) / total
            kept = num_mismatches <= max_mismatches
            writer.write(
                "\t".join(
                    [
                        f"{num_mismatches}",
                        f"{count:,d}",
                        f"{frac_at_num_mismatches:.4f}",
                        f"{frac_le_num_mismatches:.4}",
                        f"{frac_ge_num_mismatches:.4f}",
                        f"{kept}",
                    ]
                )
                + "\n"
            )
