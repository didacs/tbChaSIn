import logging
from pathlib import Path
from typing import List
from typing import Optional
from typing import Set

from fgpyo.sam import Cigar
from fgpyo.sam import CigarOp
from fgpyo.sam import Template
from pysam import AlignedSegment

from pytomebio.core.sites import FindSitesMetric
from pytomebio.core.sites import Side
from pytomebio.core.sites import SiteKey
from pytomebio.core.sites import SitesGenerator


class CrypticSeqSitesGenerator(SitesGenerator):
    _MAPPED_CIGAR_OPS: Set = {CigarOp.M, CigarOp.EQ or CigarOp.X}

    def __init__(
        self,
        donor_name: str = "attD",
        max_template_length: int = 1000,
        min_r1_mapped_bases: int = 25,
    ) -> None:
        self._logger: logging.Logger = logging.getLogger(__name__)
        self._donor_name = donor_name
        self._template_length = max_template_length
        self._min_r1_genome_mapped_bases: int = min_r1_mapped_bases

    @classmethod
    def _mapped_bases(cls, rec: AlignedSegment) -> int:
        cigar: Cigar = Cigar.from_cigartuples(rec.cigartuples)
        mapped_bases: int = sum(
            elem.length for elem in cigar.elements if elem.operator in cls._MAPPED_CIGAR_OPS
        )
        return mapped_bases

    def get_key(self, template: Template) -> Optional[SiteKey]:
        """Returns the site key after filtering.

        Ensures that the observed template length and that R1 has enough mapped bases to the genome
        """
        if abs(template.r1.template_length) > self._template_length:
            return None

        position: int
        if not template.r2.is_reverse:
            position = template.r2.reference_start
        else:
            position = template.r2.reference_end

        # Only keep if R1 has enough mapped bases to the genome
        if self._mapped_bases(rec=template.r1) < self._min_r1_genome_mapped_bases:
            return None

        # read1 always reads on the opposite strand relative to the plasmid integration
        positive_strand = template.r1.is_reverse

        return SiteKey(
            reference_name=template.r1.reference_name,
            position=position,
            attachment_site=self._donor_name,
            positive_strand=positive_strand,
        )

    def get_side(self, template: Template) -> Optional[Side]:
        return Side.Left


def find_sites(
    *,
    in_bam: Path,
    out_txt: Path,
    inter_site_slop: int = 5,
    min_mapq: int = 20,
    donor_name: str = "attD",
    max_template_length: int = 1000,
) -> None:
    """Finds integration sites by examining read pair alignments.

    Args:
        in_bam: path to the input BAM either sorted by queryname or grouped by query.
        out_txt: the list of sites found
        inter_site_slop: aggregate sites that are within this distance of each other.  The counts
                         are summed, and the coordinate is the site with the largest count.  Set to
                         zero to not aggregate counts.
        min_mapq: the minimum mapping quality to consider for a read pair
        donor_name: the name of the donor (attD) sequence in the genome.
        max_template_length: the maximum template length to allow
    """

    site_generator = CrypticSeqSitesGenerator(
        donor_name=donor_name, max_template_length=max_template_length
    )

    all_sites: List[FindSitesMetric] = site_generator.generate(
        in_bam=in_bam,
        inter_site_slop=inter_site_slop,
        min_mapq=min_mapq,
        use_duplicates=True,
        same_reference=True,
    )

    FindSitesMetric.write(out_txt, *all_sites)
