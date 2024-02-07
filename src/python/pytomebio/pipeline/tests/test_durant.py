"""Tests for the cryptic-seq pipeline"""

from pathlib import Path
from typing import Any
from typing import Dict

import yaml

from pytomebio.pipeline.tests import touch_path
from pytomebio.pipeline.tests.util import run_snakemake


def test_cryptic_seq(tmp_path: Path) -> None:
    """Basic unit test that runs the snakefile in dry-run mode to ensure it
    parses correctly.
    """

    # Set up reference data
    ref_fasta: Path = touch_path(tmp_path / "ref" / "full.fasta")
    for ext in ["amb", "ann", "bwt", "pac", "sa"]:  # BWA files
        touch_path(Path(f"{ref_fasta}.{ext}"))
    genome_fasta: Path = touch_path(tmp_path / "ref" / "genome.fasta")
    for ext in ["amb", "ann", "bwt", "pac", "sa"]:  # BWA files
        touch_path(Path(f"{genome_fasta}.{ext}"))

    # Set up FASTQs
    samples = ["foo", "bar", "two"]
    fq_dir: Path = tmp_path / "fastqs"
    for sample in samples:
        touch_path(fq_dir / f"{sample}_R1_001.fastq.gz")
        touch_path(fq_dir / f"{sample}_R2_001.fastq.gz")

    rules: Dict[str, int] = {
        "align_full": 3,
        "align_genome": 3,
        "all": 1,
        "fastp": 3,
        "filter_reads": 3,
        "find_sites": 3,
        "flagstat": 6,
        "samtools_import": 3,
        "trim_leading_r2": 3,
    }

    config: Dict[str, Any] = {
        "settings": [
            {
                "name": "first",
                "fq_dir": f"{fq_dir}",
                "ref_fasta": f"{ref_fasta}",
                "genome_fasta": f"{genome_fasta}",
                "min_aln_score": 20,
                "inter_site_slop": 10,
                "samples": [
                    {
                        "name": "foo",
                        "replicate": 1,
                        "fq1": f"{fq_dir}/foo_R1_001.fastq.gz",
                        "fq2": f"{fq_dir}/foo_R2_001.fastq.gz",
                        "stagger": "A",
                        "donor_inner_primer": "CAGCGAGTCAGTGAGCGAGG",
                        "umi_length": 0,
                        "r1_adapter": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
                        "r2_adapter": "CTGTCTCTTATACACATCTGACGCTGCCGACGA",
                    }
                ],
            },
            {
                "name": "second",
                "fq_dir": f"{fq_dir}",
                "ref_fasta": f"{ref_fasta}",
                "genome_fasta": f"{genome_fasta}",
                "min_aln_score": 15,
                "inter_site_slop": 20,
                "samples": [
                    {
                        "name": "bar",
                        "replicate": 1,
                        "fq1": f"{fq_dir}/bar_R1_001.fastq.gz",
                        "fq2": f"{fq_dir}/bar_R2_001.fastq.gz",
                        "stagger": "",
                        "donor_inner_primer": "GATTACA",
                        "umi_length": 0,
                        "r1_adapter": "ATTTATATAT",
                        "r2_adapter": "TATATATATAT",
                    },
                    {
                        "name": "two",
                        "replicate": 2,
                        "fq1": f"{fq_dir}/two_R1_001.fastq.gz",
                        "fq2": f"{fq_dir}/two_R2_001.fastq.gz",
                        "stagger": "T",
                        "donor_inner_primer": "ACGT",
                        "umi_length": 0,
                        "r1_adapter": "CGCGCGCGC",
                        "r2_adapter": "CGCGCGCGC",
                    },
                ],
            },
        ],
    }

    config_yml: Path = tmp_path / "config.yml"
    with config_yml.open("w") as writer:
        yaml.dump(config, writer)
    with config_yml.open("r") as reader:
        for line in reader:
            print(line.rstrip())

    run_snakemake(
        pipeline="durant",
        workdir=tmp_path,
        rules=rules,
        config={"config_yml": config_yml},
    )
