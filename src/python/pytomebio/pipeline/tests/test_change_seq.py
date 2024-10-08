"""Tests for the change-seq pipeline"""

from pathlib import Path
from typing import Any
from typing import Dict

import yaml

from pytomebio.pipeline.tests import touch_path
from pytomebio.pipeline.tests.util import run_snakemake


def test_change_seq(tmp_path: Path) -> None:
    """Basic unit test that runs the snakefile in dry-run mode to ensure it
    parses correctly.
    """

    # Set up reference data
    ref_fasta: Path = touch_path(tmp_path / "ref" / "ref.fasta")
    for ext in ["amb", "ann", "bwt", "pac", "sa"]:  # BWA files
        touch_path(Path(f"{ref_fasta}.{ext}"))

    # Set up FASTQs
    samples = ["foo", "bar", "two"]
    fq_dir: Path = tmp_path / "fastqs"
    for sample in samples:
        touch_path(fq_dir / f"{sample}_R1_001.fastq.gz")
        touch_path(fq_dir / f"{sample}_R2_001.fastq.gz")

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
                "fq_dir": f"{fq_dir}",
                "ref_fasta": f"{ref_fasta}",
                "attachment_sites": [
                    "attB:CACCACGCGTGGCCGGCTTGTCGACGACGGCG:GT:CTCCGTCGTCAGGATCATCCGGGGATCCCGGG"
                ],
                "samples": [
                    {
                        "name": "foo",
                        "replicate": 1,
                        "fq1": f"{fq_dir}/foo_R1_001.fastq.gz",
                        "fq2": f"{fq_dir}/foo_R2_001.fastq.gz",
                    }
                ],
            },
            {
                "name": "second",
                "fq_dir": f"{fq_dir}",
                "ref_fasta": f"{ref_fasta}",
                "attachment_sites": [
                    "attP:GCCGCTAGCGGTGGTTTGTCTGGTCAACCACCGCG:GT:"
                    "GACCGGTAGCTGGGTTTGTACCGTACACCACTGAG"
                ],
                "samples": [
                    {
                        "name": "bar",
                        "replicate": 1,
                        "fq1": f"{fq_dir}/bar_R1_001.fastq.gz",
                        "fq2": f"{fq_dir}/bar_R2_001.fastq.gz",
                    },
                    {
                        "name": "two",
                        "replicate": 2,
                        "fq1": f"{fq_dir}/two_R1_001.fastq.gz",
                        "fq2": f"{fq_dir}/two_R2_001.fastq.gz",
                    },
                ],
            },
        ],
    }

    config_yml: Path = tmp_path / "config.yml"
    with config_yml.open("w") as writer:
        yaml.dump(config, writer, default_flow_style=False)

    run_snakemake(
        pipeline="change-seq",
        workdir=tmp_path,
        rules=rules,
        config={"config_yml": config_yml},
    )
