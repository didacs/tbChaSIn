"""Tests for the change-seq pipeline"""

from pathlib import Path
from py._path.local import LocalPath as TmpDir
from typing import Any
from typing import Dict

from pytomebio.pipeline.tests.util import run_snakemake


def _touch_path(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w"):
        pass
    return path


def test_change_seq(tmpdir: TmpDir) -> None:
    """Basic unit test that runs the snakefile in dry-run mode to ensure it
    parses correctly.
    """

    # Set up reference data
    ref_fasta: Path = _touch_path(Path(tmpdir) / "ref" / "ref.fasta")
    for ext in ["amb", "ann", "bwt", "pac", "sa"]:  # BWA files
        _touch_path(Path(f"{ref_fasta}.{ext}"))

    # Set up FASTQs
    samples = ["foo", "bar", "two"]
    fq_dir: Path = Path(tmpdir) / "fastqs"
    for sample in samples:
        _touch_path(fq_dir / f"{sample}_R1_001.fastq.gz")
        _touch_path(fq_dir / f"{sample}_R2_001.fastq.gz")

    rules: Dict[str, int] = {
        "align": 3,
        "all": 1,
        "collate_sites": 1,
        "fastq_to_bam": 3,
        "fastqc": 6,
        "fgsv_aggregatesvpileup": 3,
        "fgsv_svpileup": 3,
        "find_sites": 3,
        "mark_duplicates": 3,
        "multiqc": 1,
        "picard_collect_alignment_summary_metrics": 3,
        "picard_collect_multiple_metrics": 3,
        "trim_for_tn5me": 3,
    }

    config: Dict[str, Any] = {
        "settings": [
            {
                "name": "first",
                "fq_dir": fq_dir,
                "ref_fasta": ref_fasta,
                "samples": ["foo"],
            },
            {
                "name": "second",
                "fq_dir": fq_dir,
                "ref_fasta": ref_fasta,
                "samples": ["bar", "two"],
            },
        ],
    }

    run_snakemake(
        pipeline="change-seq",
        workdir=tmpdir,
        rules=rules,
        config=config,
    )
