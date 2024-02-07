import logging
from pathlib import Path
from typing import Iterator
from typing import List
from typing import Optional
from typing import Set

from fgpyo import sam
from fgpyo.sam import Cigar
from fgpyo.sam import CigarOp
from fgpyo.sam import Template
from pysam import AlignedSegment
from pysam import AlignmentFile


class Filterer:
    _MAPPED_CIGAR_OPS: Set = {CigarOp.M, CigarOp.EQ or CigarOp.X}

    def __init__(
        self,
        donor_name: str,
        max_r1_donor_mapped_bases: int = 55,
    ) -> None:
        self.num_templates_genome_disagree: int = 0
        self._logger: logging.Logger = logging.getLogger(__name__)
        self._donor_name: str = donor_name
        self._max_r1_donor_mapped_bases: int = max_r1_donor_mapped_bases

    @classmethod
    def _mapped_bases(cls, rec: AlignedSegment) -> int:
        cigar: Cigar = Cigar.from_cigartuples(rec.cigartuples)
        mapped_bases: int = sum(
            elem.length for elem in cigar.elements if elem.operator in cls._MAPPED_CIGAR_OPS
        )
        return mapped_bases

    def _get_donor(self, template: Template, from_r1: bool) -> Optional[AlignedSegment]:
        recs: List[AlignedSegment]
        if from_r1:
            recs = [template.r1] + template.r1_supplementals
        else:
            recs = [template.r2] + template.r2_supplementals
        rec_to_return: Optional[AlignedSegment] = next(
            (rec for rec in recs if rec.reference_name == self._donor_name), None
        )
        return rec_to_return

    def keep_read(self, template: Template) -> bool:
        """Returns true if the template should be kept

        - R1 and R2 must each have at least one alignment to the genome (non-donor)
        - If R1 has an alignment to the donor, it cannot have too many mapped bases to the donor,
          otherwise it is assumed to be a linear plasmid template
        - R2 must have an alignment to the donor on the forward strand
        """

        # Both R1 and R2 must have genome alignments
        r1_genomes = [
            r1
            for r1 in [template.r1] + template.r1_supplementals
            if r1.reference_name != self._donor_name
        ]
        if len(r1_genomes) == 0:
            return False
        r2_genomes = [
            r2
            for r2 in [template.r2] + template.r2_supplementals
            if r2.reference_name != self._donor_name
        ]
        if len(r2_genomes) == 0:
            return False

        # Do not keep R1 if it has too many mapped bases to the donor
        r1_donor: Optional[AlignedSegment] = self._get_donor(template=template, from_r1=True)
        if (
            r1_donor is not None
            and self._mapped_bases(rec=r1_donor) >= self._max_r1_donor_mapped_bases
        ):
            return False

        # R2 should map to the donor on the forward strand
        r2_donor: Optional[AlignedSegment] = self._get_donor(template=template, from_r1=False)
        if r2_donor is None or r2_donor.is_reverse:
            return False

        return True


def filter_reads(
    *, in_bam: Path, keep_bam: Path, reject_bam: Path, donor_name: str = "attD"
) -> None:
    """Filters templates for reads with evidence of cryptic integration.

    Assumes reads have been mapped to the genome containing the donor sequence (attD)

    Finds:
    - R1 has an alignment to the genome
    - R1 has too many bases mapped to the donor (if present)
    - R2 has an alignment to the genome
    - R2 should map to the donor on the forward strand

    Args:
        in_bam: the input mapped BAM
        keep_bam: the output BAM with reads to keep
        reject_bam: the output BAM with reads that were rejected
        donor_name: the name of the donor (attD) sequence in the genome.
    """
    logger = logging.getLogger(__name__)

    filterer = Filterer(donor_name=donor_name)

    with AlignmentFile(str(in_bam), check_sq=False) as reader, sam.writer(
        keep_bam, header=reader.header
    ) as keep_writer, sam.writer(reject_bam, header=reader.header) as reject_writer:
        record_number: int = 1

        template_iter: Iterator[Template] = Template.iterator(reader)

        for template in template_iter:
            assert template.r1 is not None, f"Missing r1 for {template.name}"
            assert template.r2 is not None, f"Missing r2 for {template.name}"
            template.validate()

            if filterer.keep_read(template=template):
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
