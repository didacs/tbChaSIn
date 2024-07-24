import json
import logging
from collections import Counter
from pathlib import Path
from typing import Iterator
from typing import List
from typing import Optional

from Bio import Align
from fgpyo import sam
from fgpyo.sam import Template
from fgpyo.sequence import reverse_complement
from pysam import AlignmentFile

from pytomebio.core.aligner import get_query_prefix_aligner
from pytomebio.core.attachment_site import AttachmentSite
from pytomebio.core.attachment_site import AttachmentSiteMatch


def get_best_alignment(alignments: List[Align.PairwiseAlignments]) -> Align.PairwiseAlignments:
    scores = [alignments[idx].score for idx in range(0, len(alignments))]
    best_alignment = [alignment for alignment in alignments if alignment.score == max(scores)]
    return best_alignment[0]


def _num_mismatches(s1: str, s2: str) -> int:
    """Returns the number of mismatches between two sequences"""
    return sum(1 for b1, b2 in zip(s1, s2) if b1 != b2)


def aligns_better_to_site(
    query_prefix_aligner: Align.PairwiseAligner,
    r2_bases: str,
    best_match: AttachmentSiteMatch,
    min_score: int,
) -> bool:
    site: AttachmentSite = best_match.site
    site_bases = site.left + site.overhang + reverse_complement(site.right)
    if not best_match.is_left:
        site_bases = reverse_complement(site_bases)
    alignments = query_prefix_aligner.align(r2_bases, site_bases)
    if len(alignments) == 0:
        return False
    alignment = get_best_alignment(alignments=alignments)
    return alignment.score >= min_score


def trim_leading_attachment_site(
    *,
    in_bam: Path,
    keep_bam: Path,
    reject_bam: Optional[Path] = None,
    full_site_bam: Optional[Path] = None,
    keep_full_sites: bool = False,
    out_metrics: Path,
    out_summary: Optional[Path] = None,
    attachment_site: List[AttachmentSite],
    max_mismatches: int = 1,
) -> None:
    """Finds a leading left (right) side of an attachment site at the end (start) of R2 and
    trims it.

    Any read pair that where the attachment site left/right is not found at the start of R2 is
    rejected.  If an attachment site left/right is found, then full attachment site
    (left+overhang+right) is aligned to the start of the read, and if found, the read is kept but
    neither trimmed nor has the attachment-site-specific SAM tags added (see below).

    The "rm" and "mm" SAM tags store information about which attachment site, if any, the given
    read and mate matched.

    Args:
        in_bam: the input unmapped BAM
        keep_bam: the output BAM with trimmed R2s.
        reject_bam: the output BAM with reads where the leading attachment site in R2 was not found
        full_site_bam: optional; the output BAM with templates where R2 contains the full
                       attachment site
        keep_full_sites: whether to keep reads where R2 contains the full attachment site; if True,
                         full-site reads are written to '--full-site-bam' if specified, otherwise
                         '--keep-bam'.
        out_metrics: the path of the output file with trimming metrics
        out_summary: optional; the path of the output JSON file with summary statistics
        attachment_site: one or more attachment sites (e.g. attB, attP).  Each
                         attachment site is a colon-delimited quadruple of
                         `<name>:<left_seq>:<overhang>:<right_seq>` where
                         `<left_seq><overhang><right_seq>`
                         is the full attachment site sequence.
        max_mismatches: the maximum number of mismatches to allow
    """
    logger = logging.getLogger(__name__)

    # assert min_score >= 0 since we count unmatched reads as -1
    assert max_mismatches >= 0, "--max-mismatches must be >= 0"
    assert len(attachment_site) > 0, "no attachment sites given"

    aligner = get_query_prefix_aligner()

    num_mismatches_counter = Counter({0: 0})
    record_number: int = 1
    num_templates = 0
    num_full_sites = 0
    num_rejected = 0

    with (
        AlignmentFile(str(in_bam), check_sq=False) as reader,
        sam.writer(keep_bam, header=reader.header) as keep_writer,
    ):
        # open optional writers
        reject_writer = None
        close_reject_writer = False
        if reject_bam is not None:
            reject_writer = sam.writer(reject_bam, header=reader.header)
            close_reject_writer = True

        full_site_writer = None
        close_full_site_writer = False
        if keep_full_sites:
            if full_site_bam is not None:
                full_site_writer = sam.writer(full_site_bam, header=reader.header)
                close_full_site_writer = True
            else:
                full_site_writer = keep_writer

        template_iter: Iterator[Template] = Template.iterator(reader)

        for template in template_iter:
            assert template.r1 is not None, f"Missing r1 for {template.name}"
            assert template.r2 is not None, f"Missing r2 for {template.name}"
            template.validate()

            num_templates += 1

            r2_bases = template.r2.query_sequence

            best_match: Optional[AttachmentSiteMatch] = None
            best_num_mismatches: Optional[int] = None
            for site in attachment_site:
                num_mismatches = _num_mismatches(site.left, r2_bases)
                if best_num_mismatches is None or num_mismatches < best_num_mismatches:
                    best_match = AttachmentSiteMatch(
                        site=site,
                        is_left=True,
                        score=len(site.left) - _num_mismatches(site.left, r2_bases),
                        alignment_length=len(site.left),
                    )
                    best_num_mismatches = num_mismatches
                num_mismatches = _num_mismatches(site.right, r2_bases)
                if best_num_mismatches is None or num_mismatches < best_num_mismatches:
                    best_match = AttachmentSiteMatch(
                        site=site,
                        is_left=False,
                        score=len(site.right) - _num_mismatches(site.right, r2_bases),
                        alignment_length=len(site.right),
                    )
                    best_num_mismatches = num_mismatches

            assert best_match is not None, "a best match should always be found"

            # If we not only align to the leading attachment site, but also the full site (e.g.
            # full attP or attB), we may keep the record (if `keep_full_sites is True`) but skip
            # adding the TAG
            aligns_to_full_site = False
            if best_num_mismatches <= max_mismatches:
                best_score = best_match.alignment_length - best_num_mismatches
                aligns_to_full_site = aligns_better_to_site(
                    query_prefix_aligner=aligner,
                    r2_bases=r2_bases,
                    best_match=best_match,
                    min_score=best_score,
                )

            num_mismatches_counter[best_num_mismatches] += 1

            # If R2 aligns fully to the full attachment site, we may still keep the read (if
            # `keep_full_sites is True`), but if so we do not add the "attachment site match"
            # tags. The `find-sites` tool will later filter out reads without this tag (as a
            # proxy for those that match the leading but not full attachment site). Such reads
            # are kept in case we wish to know the number of reads that later align to the full
            # attB or attP.

            template.r1.set_tag(AttachmentSiteMatch.READ_MATCH_TAG, "None")
            template.r1.set_tag(AttachmentSiteMatch.MATE_MATCH_TAG, "None")
            template.r2.set_tag(AttachmentSiteMatch.READ_MATCH_TAG, "None")
            template.r2.set_tag(AttachmentSiteMatch.MATE_MATCH_TAG, "None")

            if aligns_to_full_site:
                num_full_sites += 1
                writer = full_site_writer
            elif best_num_mismatches > max_mismatches:
                num_rejected += 1
                writer = reject_writer
            else:
                template.r1.set_tag(AttachmentSiteMatch.MATE_MATCH_TAG, best_match.to_sam_tag())
                template.r2.set_tag(AttachmentSiteMatch.READ_MATCH_TAG, best_match.to_sam_tag())
                r2_quals = template.r2.query_qualities
                template.r2.query_sequence = r2_bases[best_match.alignment_length :]
                template.r2.query_qualities = r2_quals[best_match.alignment_length :]
                writer = keep_writer

            for rec in template.all_recs():
                if writer is not None:
                    writer.write(rec)
                if record_number % 100000 == 0:
                    logger.info(f"Processed {record_number:,d} records")
                record_number += 1

        # close optional writers
        if close_reject_writer:
            reject_writer.close()
        if close_full_site_writer:
            full_site_writer.close()

        logger.info(f"Processed {record_number:,d} records")

    # write out the metrics
    with out_metrics.open("w") as metrics_writer:
        total = sum(num_mismatches_counter.values())
        running_sum = 0
        metrics_writer.write(
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
        for num_mismatches in range(max(num_mismatches_counter.keys())):
            count = num_mismatches_counter.get(num_mismatches, 0)
            running_sum += count
            frac_at_num_mismatches = float(count) / total
            frac_ge_num_mismatches = float(total - running_sum) / total
            frac_le_num_mismatches = float(running_sum) / total
            kept = num_mismatches <= max_mismatches
            metrics_writer.write(
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

    if out_summary is not None:
        with out_summary.open("w") as summary_writer:
            json.dump(
                {
                    "total_records": record_number,
                    "total_templates": num_templates,
                    "full_sites": num_full_sites,
                    "rejected": num_rejected,
                },
                summary_writer,
            )
