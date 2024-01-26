from pathlib import Path
from typing import List
from typing import Optional

from fgpyo.sam import Cigar
from fgpyo.sam import CigarElement
from fgpyo.sam import CigarOp
from fgpyo.sam import Template
from pysam import AlignedSegment

from pytomebio.core.attachment_site import AttachmentSiteMatch
from pytomebio.core.sites import FindSitesMetric
from pytomebio.core.sites import Side
from pytomebio.core.sites import SiteKey
from pytomebio.core.sites import SitesGenerator


class CrypticSeqSitesGenerator(SitesGenerator):
    def __init__(self, r2_max_start_soft_clip: int = 0) -> None:
        self.r2_max_start_soft_clip = r2_max_start_soft_clip

    def get_key(self, template: Template) -> Optional[SiteKey]:
        read2: AlignedSegment = template.r2

        # Do not consider reads where R2 does not have a leading attachment site mathc
        if not read2.has_tag(AttachmentSiteMatch.READ_MATCH_TAG):
            return None

        # Do not consider reads that have too many soft-clipped bases at the start of the read
        cigar: Cigar = Cigar.from_cigartuples(read2.cigartuples)
        elem: CigarElement = cigar.elements[-1] if read2.is_reverse else cigar.elements[0]
        r2_start_soft_clip: int = elem.length if elem.operator == CigarOp.S else 0
        if r2_start_soft_clip > self.r2_max_start_soft_clip:
            return None

        tag_value: str = str(read2.get_tag(AttachmentSiteMatch.READ_MATCH_TAG))
        site_match = AttachmentSiteMatch.from_sam_tag(tag_value)

        position = read2.reference_end if read2.is_reverse else read2.reference_start

        # if R2 has "left side" (e.g. P) then use alignment strand, otherwise negate the strand
        positive_strand = site_match.is_left == (not read2.is_reverse)

        return SiteKey(
            reference_name=read2.reference_name,
            position=position,
            attachment_site=site_match.site.name,
            positive_strand=positive_strand,
        )

    def get_side(self, template: Template) -> Optional[Side]:
        read2: AlignedSegment = template.r2
        # Only one side of the integration is observed in Cryptic-Seq (at the start of R2)
        tag_value: str = str(read2.get_tag(AttachmentSiteMatch.READ_MATCH_TAG))
        site_match = AttachmentSiteMatch.from_sam_tag(tag_value)
        return Side.Left if site_match.is_left else Side.Right


def find_sites(
    *,
    in_bam: Path,
    out_txt: Path,
    r2_max_start_soft_clip: int = 0,
    inter_site_slop: int = 5,
    min_mapq: int = 20,
) -> None:
    """Finds integration sites by examining read pair alignments.

    Reads are assumed to have UMIs removed.

    The R2s are assumed to have the leading attachment site trimmed, without the overhang
    trimmed (see tomebio-tools cryptic-seq trim-leading-attachment-site), as well as having
    sequence trimmed including and after any Tn5 mosaic end found.

    R1s are assumed to have the Tn5 mosaic end trimmed from the start of the read.

    Since R1s are always the read that is tagmented, the first base sequenced in the R2s identify
    the integration site.  Thus, each read pair is examined and a site is called if both reads map
    to the same contig.  Next, the counts for each integration site in the genome is aggregated.
    Finally, sites within inter-site-slop distance from each other are aggregated.

    For sites that are aggregated, the position is the site with the highest read count.

    Args:
        in_bam: path to the input BAM either sorted by queryname or grouped by query.
        out_txt: the list of sites found
        r2_max_start_soft_clip: the maximum number of soft-clipped bases at the start of R2 in
                                sequencing order.
        inter_site_slop: aggregate sites that are within this distance of each other.  The counts
                         are summed, and the coordinate is the site with the largest count.  Set to
                         zero to not aggregate counts.
        min_mapq: the minimum mapping quality to consider for a read pair
    """

    site_generator = CrypticSeqSitesGenerator(r2_max_start_soft_clip=r2_max_start_soft_clip)

    all_sites: List[FindSitesMetric] = site_generator.generate(
        in_bam=in_bam,
        inter_site_slop=inter_site_slop,
        min_mapq=min_mapq,
        use_duplicates=False,
        same_reference=True,
    )

    FindSitesMetric.write(out_txt, *all_sites)
