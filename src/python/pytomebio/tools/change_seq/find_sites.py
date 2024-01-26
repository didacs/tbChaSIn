from pathlib import Path
from typing import List
from typing import Optional

from fgpyo.sam import Template
from pysam import AlignedSegment

from pytomebio.core.attachment_site import AttachmentSiteMatch
from pytomebio.core.sites import FindSitesMetric
from pytomebio.core.sites import Side
from pytomebio.core.sites import SiteKey
from pytomebio.core.sites import SitesGenerator


class ChangeSeqSitesGenerator(SitesGenerator):
    def __init__(self, intra_read_slop: int) -> None:
        self.intra_read_slop = intra_read_slop

    def get_key(self, template: Template) -> Optional[SiteKey]:
        read1: AlignedSegment = template.r1
        read2: AlignedSegment = template.r2
        read1_pos = read1.reference_end if read1.is_reverse else read1.reference_start
        read2_pos = read2.reference_end if read2.is_reverse else read2.reference_start
        pos_diff = abs(read1_pos - read2_pos)
        if pos_diff > self.intra_read_slop:
            return None
        tag_value: str = str(read1.get_tag(AttachmentSiteMatch.READ_MATCH_TAG))
        site_match_r1 = AttachmentSiteMatch.from_sam_tag(tag_value)
        position = (read1_pos + read2_pos) // 2
        # whichever read pair has the "left side" (e.g. P) defines the strand
        positive_strand = not (read1.is_reverse if site_match_r1.is_left else read2.is_reverse)

        return SiteKey(
            reference_name=read1.reference_name,
            position=position,
            attachment_site=site_match_r1.site.name,
            positive_strand=positive_strand,
        )

    def get_side(self, template: Template) -> Optional[Side]:
        # for Change-Seq, both sides of the integration are observed in one template, so always
        # return `Side.Both`.
        return Side.Both


def find_sites(
    *,
    in_bam: Path,
    out_txt: Path,
    intra_read_slop: int = 5,
    inter_site_slop: int = 5,
    min_mapq: int = 20,
) -> None:
    """Finds integration sites by examining read pair alignments.

    The read pairs are assumed to have the leading attachment site trimmed, without the overhang
    trimmed (see tomebio-tools change-seq fastq-to-bam).

    The integration site should occur at the first sequenced base in each read in a read pair. Thus
    each read pair is examined and a site is called if the genomic position of the first sequenced
    base in each read pair is within the given intra-read-slop.  Next, the counts for each
    integration site in the genome is aggregated.  Finally, sites within inter-site-slop distance
    from each other are aggregated.

    For sites that are aggregated, the position is the site with the highest read count.

    Args:
        in_bam: path to the input BAM either sorted by queryname or grouped by query.
        out_txt: the list of sites found
        intra_read_slop: the maximum absolute difference between the first base sequenced within
                         a read pair.  Typically, at least 2bp to account for the overhang in the
                         attachment site, which is not trimmed.
        inter_site_slop: aggregate sites that are within this distance of each other.  The counts
                         are summed, and the coordinate is the site with the largest count.  Set to
                         zero to not aggregate counts.
        min_mapq: the minimum mapping quality to consider for a read pair
    """

    site_generator = ChangeSeqSitesGenerator(intra_read_slop=intra_read_slop)

    all_sites: List[FindSitesMetric] = site_generator.generate(
        in_bam=in_bam,
        inter_site_slop=inter_site_slop,
        min_mapq=min_mapq,
        use_duplicates=True,
        same_reference=True,
    )

    FindSitesMetric.write(out_txt, *all_sites)
