"""Main entry point for all PYTOMEBIO tools."""

import logging
import sys
from typing import Callable
from typing import Dict
from typing import List
from typing import Optional

import defopt

from pytomebio.tools.change_seq.fastq_to_bam import AttachmentSite
from pytomebio.tools.change_seq.fastq_to_bam import fastq_to_bam as change_seq_fastq_to_bam
from pytomebio.tools.change_seq.find_sites import find_sites as change_seq_find_sites
from pytomebio.tools.change_seq.trim_for_tn5me import trim_for_tn5me as change_seq_trim_for_tn5me
from pytomebio.tools.common.collate_sites import collate_sites as common_collate_sites
from pytomebio.tools.cryptic_seq.create_config_from_metasheet import (
    create_config_from_metasheet as cryptic_seq_create_config_from_metasheet,
)
from pytomebio.tools.cryptic_seq.create_metasheet_from_benchling import (
    create_metasheet_from_benchling as cryptic_seq_create_metasheet_from_benchling,
)
from pytomebio.tools.cryptic_seq.find_sites import find_sites as cryptic_seq_find_sites
from pytomebio.tools.cryptic_seq.resolve_references import (
    resolve_references as cryptic_seq_resolve_references,
)
from pytomebio.tools.cryptic_seq.trim_for_tn5me import trim_for_tn5me as cryptic_seq_trim_for_tn5me
from pytomebio.tools.cryptic_seq.trim_leading_attachment_site import (
    trim_leading_attachment_site as cryptic_trim_leading_attachment_site,
)
from pytomebio.tools.durant.filter_reads import filter_reads as durant_filter_reads
from pytomebio.tools.durant.find_sites import find_sites as durant_find_sites
from pytomebio.tools.durant.trim_leading_r2 import trim_leading_r2 as durant_trim_leading_r2

TOOLS: Dict[str, List[Callable]] = {
    "common": [common_collate_sites],
    "change-seq": [
        change_seq_fastq_to_bam,
        change_seq_find_sites,
        change_seq_trim_for_tn5me,
    ],
    "cryptic-seq": [
        cryptic_seq_create_metasheet_from_benchling,
        cryptic_seq_create_config_from_metasheet,
        cryptic_seq_find_sites,
        cryptic_seq_resolve_references,
        cryptic_seq_trim_for_tn5me,
        cryptic_trim_leading_attachment_site,
    ],
    "durant": [durant_filter_reads, durant_find_sites, durant_trim_leading_r2],
}


def main(argv: Optional[List[str]] = None) -> None:
    if argv is None:
        argv = sys.argv[1:]
    logger = logging.getLogger(__name__)
    if len(argv) != 0 and all(arg not in argv for arg in ["-h", "--help"]):
        logger.info("Running command: tomebio-tools " + " ".join(argv))
    try:
        defopt.run(
            funcs=TOOLS, argv=argv, parsers={AttachmentSite: AttachmentSite.parse_attachment_site}
        )
        logger.info("Completed successfully.")
    except Exception as e:
        logger.info("Failed on command: " + " ".join(argv))
        raise e
