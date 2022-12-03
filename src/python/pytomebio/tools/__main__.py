"""Main entry point for all PYTOMEBIO tools."""

import defopt
import logging
import sys
from typing import Callable
from typing import Dict
from typing import List
from typing import Optional

from pytomebio.tools.change_seq.collate_sites import collate_sites as change_seq_collate_sites
from pytomebio.tools.change_seq.fastq_to_bam import fastq_to_bam as change_seq_fastq_to_bam
from pytomebio.tools.change_seq.find_sites import find_sites as change_seq_find_sites
from pytomebio.tools.change_seq.trim_for_tn5me import trim_for_tn5me as change_seq_trim_for_tn5me

TOOLS: Dict[str, List[Callable]] = {
    "change-seq": [
        change_seq_collate_sites,
        change_seq_fastq_to_bam,
        change_seq_find_sites,
        change_seq_trim_for_tn5me,
    ]
}


def main(argv: Optional[List[str]] = None) -> None:
    if argv is None:
        argv = sys.argv[1:]
    logger = logging.getLogger(__name__)
    if len(argv) != 0 and all(arg not in argv for arg in ["-h", "--help"]):
        logger.info("Running command: tomebio-tools " + " ".join(argv))
    try:
        defopt.run(funcs=TOOLS, argv=argv)
        logger.info("Completed successfully.")
    except Exception as e:
        logger.info("Failed on command: " + " ".join(argv))
        raise e
