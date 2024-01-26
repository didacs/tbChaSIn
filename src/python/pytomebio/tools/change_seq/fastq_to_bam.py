import logging
import pysam
from Bio import Align
from attr import define
from attr import frozen
from fgpyo.util.metric import Metric
from pathlib import Path
from pysam import AlignedSegment
from pysam import AlignmentFile
from pysam import AlignmentHeader
from pysam import FastxFile
from pysam import FastxRecord
from typing import Any
from typing import Dict
from typing import List
from typing import Optional
from typing import Iterator

from pytomebio.core.attachment_site import AttachmentSite
from pytomebio.core.attachment_site import AttachmentSiteMatch
from pytomebio.core.aligner import get_query_prefix_aligner


@frozen
class UngappedAligner:
    """Utility class for ungapped alignment.
    Attributes:
        sites: the list of sites to align to
        min_score: the minimum alignment score to keep an alignment
        match_score: the match score (positive)
        mismatch_score: the mismatch score (negative)
        max_read_start_offset: the maximum offset from the start of the read to start the alignment
        max_target_start_offset: the maximum offset from the start of the target to start the
                                 alignment
    """

    sites: List[AttachmentSite]
    min_score: int
    match_score: int
    mismatch_score: int
    max_read_start_offset: int
    max_target_start_offset: int

    def _iter_attachment_site_matches(self, rec: FastxRecord) -> Iterator[AttachmentSiteMatch]:
        """Iterate over all allowed ungapped AttachmentSiteMatches so the optimal alignement can be
        chosen from among them.

        Args:
            rec: the read to align
        """
        # Iterate over all allowed ungapped AttachmentSiteMatches
        full_read = rec.sequence.upper()
        for site in self.sites:
            for read_start_offset in range(0, self.max_read_start_offset + 1):
                read = full_read[read_start_offset:]
                for target_start_offset in range(0, self.max_target_start_offset + 1):
                    # left side of the side
                    target_left = site.left[target_start_offset:]
                    yield AttachmentSiteMatch(
                        site=site,
                        is_left=True,
                        score=self.score(read, target_left),
                        alignment_length=min(len(read), len(target_left)),
                        read_start_offset=read_start_offset,
                    )

                    # right side of the site
                    target_right = site.right[target_start_offset:]
                    yield AttachmentSiteMatch(
                        site=site,
                        is_left=False,
                        score=self.score(read, target_right),
                        alignment_length=min(len(read), len(target_right)),
                        read_start_offset=read_start_offset,
                    )

    def align(self, rec: FastxRecord) -> Optional[AttachmentSiteMatch]:
        """Aligns the given read to both the left and right sides of each site, as well as up to
        the maximum read and target start offset respectively.  Returns the best alignment with
        at least the minimum score.

        Args:
            rec: the read to align
        """
        # find best match by maximizing score, using alignment_length as a tie-breaker:
        best = max(
            self._iter_attachment_site_matches(rec=rec),
            default=None,
            key=lambda attachment_site_match: (
                attachment_site_match.score,
                attachment_site_match.alignment_length,
            ),
        )

        return None if best is None or best.score < self.min_score else best

    def score(self, s1: str, s2: str) -> int:
        """Returns the ungapped alignment score of the two sequences"""
        score = 0
        for b1, b2 in zip(s1, s2):
            score += self.match_score if b1 == b2 else self.mismatch_score
        return score


def to_aligned_segment(
    fq: FastxRecord,
    header: AlignmentHeader,
    read_group: str,
    is_r1: bool,
    read_match: Optional[AttachmentSiteMatch],
    mate_match: Optional[AttachmentSiteMatch],
) -> AlignedSegment:
    """Converts the given FASTQ record into SAM, adding appropriate SAM tags.

    Args:
        fq: the FASTQ record
        header: the output SAM header
        read_group: the read group identifier to use for the RG SAM tag
        is_r1: true if the read is read one in a pair, false if read two
        read_match: the match for the current read
        mate_match: the match for the mate of the current read
    """
    rec = AlignedSegment(header=header)

    bases = fq.sequence
    quals = fq.quality
    if read_match is not None:
        bases = bases[read_match.alignment_length :]
        quals = quals[read_match.alignment_length :]

    # flags
    rec.is_paired = True
    rec.is_read1 = is_r1
    rec.is_read2 = not is_r1
    rec.is_qcfail = False
    rec.is_duplicate = False
    rec.is_secondary = False
    rec.is_supplementary = False
    rec.is_reverse = False
    rec.is_unmapped = True

    # Basic attributes
    rec.query_name = fq.name
    rec.reference_name = None
    rec.reference_start = -1
    rec.mapping_quality = 0
    rec.query_sequence = bases
    rec.query_qualities = pysam.qualitystring_to_array(quals)

    # Tags
    attrs = dict()
    attrs["RG"] = read_group
    attrs[AttachmentSiteMatch.READ_MATCH_TAG] = (
        "None" if read_match is None else read_match.to_sam_tag()
    )
    attrs[AttachmentSiteMatch.MATE_MATCH_TAG] = (
        "None" if mate_match is None else mate_match.to_sam_tag()
    )
    rec.set_tags(list(attrs.items()))

    return rec


def perform_gapped_rescue(
    full_read: str,
    sites: List[AttachmentSite],
    aligner: Align.PairwiseAligner,
    rescue_prefix_len: int,
    min_score: int,
) -> Optional[AttachmentSiteMatch]:
    """Performs gapped alignment of the read across the given sites.

    Only performs the alignment if the prefixes of the read and site sequence match up to the given
    length.

    Args:
        full_read: the full read sequence
        sites: the list of sites to align the read
        aligner: the glocal pairwise aligner
        rescue_prefix_len: the # of bases at the start of the read and site sequence that must
                           exactly match to perform gapped alignment
        min_score: the minimum alignment score for a returned alignment
    """
    match: Optional[AttachmentSiteMatch] = None
    assert rescue_prefix_len > 0, "rescue_prefix_len must be greater than 0"
    for site in sites:
        if full_read[:rescue_prefix_len] == site.left[:rescue_prefix_len]:
            match = gapped_alignment(
                full_read=full_read,
                site=site,
                is_left=True,
                aligner=aligner,
                min_score=min_score,
            )
        elif full_read[:rescue_prefix_len] == site.right[:rescue_prefix_len]:
            match = gapped_alignment(
                full_read=full_read, site=site, is_left=False, aligner=aligner, min_score=min_score
            )
        if match is not None:
            return match
    return match


def gapped_alignment(
    full_read: str,
    site: AttachmentSite,
    is_left: bool,
    aligner: Align.PairwiseAligner,
    min_score: int,
) -> Optional[AttachmentSiteMatch]:
    """Performs pairwise gapped alignment returning a match if found.

    Args:
        full_read: the full read sequence
        site: the site to align to
        is_left: true if aligning to the left side of the site, false if on the right
        aligner: the glocal pairwise aligner
        min_score: the minimum alignment score for a returned alignment
    """
    target = site.left if is_left else site.right
    alignments = aligner.align(full_read, target)
    if len(alignments) == 0 or alignments[0].score < min_score:
        return None
    alignment = alignments[0]
    aligned = alignment.aligned
    if len(aligned[0]) == 0:
        return None  # got an empty alignment
    # the last aligned base in the read/query is the number of bases we have to trim off
    alignment_length = aligned[0][-1][-1]
    match: AttachmentSiteMatch = AttachmentSiteMatch(
        site=site, is_left=is_left, score=alignment.score, alignment_length=alignment_length
    )
    return match


@define
class MatchingMetric(Metric["MatchingMetric"]):
    """Stores counts for how reads matched the attachment sites.

    Attributes:
        valid: the number of read pairs that matched different sides of the same attachment site
        rescued: the number of valid read pairs found due to gapped rescue
        different_site: the number of read pairs that matched different sites
        same_end: the number of read pairs that matched the same side of the same attachment site
        one_end: the number of read pairs where only one read matched an attachment site
        no_match: the number of read pairs that matched no attachment site
        rejected: the number of read pairs that are not valid
    """

    valid: int = 0
    rescued: int = 0
    different_site: int = 0
    same_end: int = 0
    one_end: int = 0
    no_match: int = 0
    rejected: int = 0

    def update(
        self,
        r1_match: Optional[AttachmentSiteMatch],
        r2_match: Optional[AttachmentSiteMatch],
        rescued_attempted: bool,
    ) -> bool:
        """Updates the counts for read matching.

        Args:
            r1_match: the match for read one
            r2_match: the match for read two
            rescued_attempted: true if gapped rescue was attempted, false if not
        """
        if r1_match is None and r2_match is None:
            self.no_match += 1
        elif r1_match is None or r2_match is None:
            self.one_end += 1
        elif r1_match.site.name != r2_match.site.name:
            self.different_site += 1
        elif r1_match.is_left == r2_match.is_left:
            self.same_end += 1
        else:
            self.valid += 1
            if rescued_attempted:
                self.rescued += 1
            return True
        self.rejected += 1
        return False


def fastq_to_bam(
    *,
    r1_fq: Path,
    r2_fq: Path,
    keep_bam: Path,
    reject_bam: Path,
    metric_tsv: Path,
    read_group: str,
    attachment_site: List[AttachmentSite],
    threads: int = 1,
    _min_score: int = 16,
    _rescue_prefix_len: int = 16,
    _max_read_start_offset: int = 3,
    _max_target_start_offset: int = 3,
    _gapped_rescue: bool = True,
) -> MatchingMetric:
    """Convert FASTQ to BAM for CHANGE-Seq.

    Matches the start of each read to each side of each attachment site using ungapped alignment.
    If a match is found for only one read in a pair, attempts to align the unmatched read pair to
    the opposite end of the matched attachment site using gapped alignment.  If still reads are not
    matched, attempts gapped alignment independent of read end to all sites when the start of the
    read matches exactly the given sequence of the attachment site's side a minimum # of bases.

    The read pairs where both reads match the same attachment site, but on different sides, are
    written to the "keep" BAM.  All other reads are written to the "reject" BAM.

    The "rm" and "mm" SAM tags store information about which attachment site, if any, the given
    read and mate matched.

    The output metric file stores counts of read pairs to describe how they matched the attachment
    sites.

    Args:
        r1_fq: the path to the input read one (R1) FASTQ.gz
        r2_fq: the path to the input read two (R2) FASTQ.gz
        keep_bam: the path to the output BAM with kept reads
        reject_bam: the path to the output BAM with discarded reads
        metric_tsv: the path to the output TSV with matching metrics
        read_group: the read group identifier to use
        attachment_site: one or more attachment sites (e.g. attB, attP).  Each
                         attachment site is a colon-delimited quadruple of
                         `<name>:<left_seq>:<overhang>:<right_seq>` where
                         `<left_seq><overhang><right_seq>`
                         is the full attachment site sequence.
        threads: threads to use for compressing the BAMs
        _min_score: the minimum alignment score to accept
        _rescue_prefix_len: require this number of leading bases between the read and site sequence
                            before trying gapped rescue
        _max_read_start_offset: the maximum offset from the start of the read for the ungapped
                                alignment to start
        _max_target_start_offset: the maximum offset from the start of the target for the ungapped
                                  alignment to start
        _gapped_rescue: true to use gapped rescue, false otherwise
    """
    # Args with underscores are fixed tool parameters not exposed to the command line, intended to
    # be changed in testing.

    logger = logging.getLogger(__name__)

    ##############################################################################
    # Definitions of the attachment sites.
    ##############################################################################
    # name,left,overhang,right
    # attB:CACCACGCGTGGCCGGCTTGTCGACGACGGCG:GT:CTCCGTCGTCAGGATCATCCGGGGATCCCGGG
    # attP:GCCGCTAGCGGTGGTTTGTCTGGTCAACCACCGCG:GT:CTCAGTGGTGTACGGTACAAACCCAGCTACCGGTC

    # check that have at least one attachment site that was successfully added to the list
    assert len(attachment_site) >= 1, "Must specify at least one attachment site"

    ##############################################################################
    # build the class to match attachment sites to the starts of reads
    ##############################################################################
    aligner = get_query_prefix_aligner()
    matcher = UngappedAligner(
        sites=attachment_site,
        min_score=_min_score,
        match_score=aligner.match_score,
        mismatch_score=aligner.mismatch_score,
        max_read_start_offset=_max_read_start_offset,
        max_target_start_offset=_max_target_start_offset,
    )

    ##############################################################################
    # Set up the output SAM/BAM
    ##############################################################################
    sam_header_dict: Dict[str, Any] = {
        "HD": {"VN": "1.6", "SO": "unsorted", "GO": "query"},
        "RG": [{"ID": read_group, "SM": read_group, "PL": "ILLUMINA"}],
    }
    sam_header: AlignmentHeader = AlignmentHeader.from_dict(sam_header_dict)

    ##############################################################################
    # main loop
    ##############################################################################
    record_number: int = 1
    metric: MatchingMetric = MatchingMetric()
    with (
        # open read1 and read2 fastq files for reading
        FastxFile(str(r1_fq)) as r1_fh,
        FastxFile(str(r2_fq)) as r2_fh,
        # open keep / reject bams for writing
        AlignmentFile(str(keep_bam), "wb", header=sam_header, threads=threads) as fh_out_keep,
        AlignmentFile(str(reject_bam), "wb", header=sam_header, threads=threads) as fh_out_reject,
    ):
        for r1, r2 in zip(r1_fh, r2_fh):
            assert r1.name == r2.name, (
                f"Error: FASTQ names disagreed for the {record_number}th entry"
                f" '{r1.name}' != '{r2.name}'"
            )

            # ungapped alignment, allowing the alignment to start offset from the start of the read
            # or target
            r1_match: Optional[AttachmentSiteMatch] = matcher.align(rec=r1)
            r2_match: Optional[AttachmentSiteMatch] = matcher.align(rec=r2)

            # Gapped alignment is expensive, so we try not do that above. Nonetheless, there are
            # some reads that would have been perfectly fine with gapped alignment
            rescued_attempted = False
            if _gapped_rescue and (r1_match is None or r2_match is None):
                rescued_attempted = True
                # If one read has a match and the other one doesn't, we could try a gapped
                # alignment to the other end of that site
                if r1_match is not None and r2_match is None:
                    r2_match = gapped_alignment(
                        full_read=r2.sequence,
                        site=r1_match.site,
                        is_left=(not r1_match.is_left),
                        aligner=aligner,
                        min_score=_min_score,
                    )
                elif r1_match is None and r2_match is not None:
                    r1_match = gapped_alignment(
                        full_read=r1.sequence,
                        site=r2_match.site,
                        is_left=(not r2_match.is_left),
                        aligner=aligner,
                        min_score=_min_score,
                    )

                # If the read matches a shorter sequence (e.g. 16bp), then we could try a gapped
                # alignment
                if r1_match is None:
                    r1_match = perform_gapped_rescue(
                        full_read=r1.sequence,
                        sites=matcher.sites,
                        aligner=aligner,
                        rescue_prefix_len=_rescue_prefix_len,
                        min_score=_min_score,
                    )
                if r2_match is None:
                    r2_match = perform_gapped_rescue(
                        full_read=r2.sequence,
                        sites=matcher.sites,
                        aligner=aligner,
                        rescue_prefix_len=_rescue_prefix_len,
                        min_score=_min_score,
                    )

            # Create the output BAM records
            r1_rec = to_aligned_segment(
                fq=r1,
                header=sam_header,
                read_group=read_group,
                is_r1=True,
                read_match=r1_match,
                mate_match=r2_match,
            )
            r2_rec = to_aligned_segment(
                fq=r2,
                header=sam_header,
                read_group=read_group,
                is_r1=False,
                read_match=r2_match,
                mate_match=r1_match,
            )

            metric.update(
                r1_match=r1_match, r2_match=r2_match, rescued_attempted=rescued_attempted
            )

            if (
                r1_match is not None
                and r2_match is not None
                and r1_match.site == r2_match.site
                and r1_match.is_left != r2_match.is_left
            ):
                fh_out_keep.write(r1_rec)
                fh_out_keep.write(r2_rec)
            else:
                fh_out_reject.write(r1_rec)
                fh_out_reject.write(r2_rec)

            if record_number % 10000 == 0:
                logger.info(
                    f"Processed {record_number:,d} records ({metric.valid:,d}"
                    f" kept, {metric.rejected:,d} rejected)"
                )
            record_number += 1

    MatchingMetric.write(metric_tsv, metric)

    logger.info(f"Read {record_number:,d} paired end reads.")
    logger.info(f"Wrote {metric.valid:,d} paired end reads to {keep_bam}")
    logger.info(f"Wrote {metric.rejected:,d} paired end reads to {reject_bam}")

    return metric
