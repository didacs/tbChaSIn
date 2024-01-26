"""Tests for the cryptic-seq pipeline"""

import yaml
from pathlib import Path
from py._path.local import LocalPath as TmpDir
from typing import Any
from typing import Dict

from pytomebio.pipeline.tests import touch_path
from pytomebio.pipeline.tests.util import run_snakemake


def test_cryptic_seq(tmpdir: TmpDir) -> None:
    """Basic unit test that runs the snakefile in dry-run mode to ensure it
    parses correctly.
    """

    # Set up reference data
    ref_fasta: Path = touch_path(Path(tmpdir) / "ref" / "ref.fasta")
    for ext in ["amb", "ann", "bwt", "pac", "sa"]:  # BWA files
        touch_path(Path(f"{ref_fasta}.{ext}"))

    # Set up FASTQs
    samples = ["foo", "bar", "two"]
    fq_dir: Path = Path(tmpdir) / "fastqs"
    for sample in samples:
        touch_path(fq_dir / f"{sample}_R1_001.fastq.gz")
        touch_path(fq_dir / f"{sample}_R2_001.fastq.gz")

    rules: Dict[str, int] = {
        "align": 3,
        "all": 1,
        "change_seq_trim_for_tn5me": 3,
        "cryptic_seq_trim_for_tn5me": 3,
        "collate_sites": 1,
        "fgbio_clip_bam": 3,
        "fgbio_fastq_to_bam": 3,
        "fastqc": 6,
        # "fgsv_aggregatesvpileup": 3,
        # "fgsv_svpileup": 3,
        "find_sites": 3,
        "mark_duplicates": 3,
        "multiqc": 1,
        "picard_collect_alignment_summary_metrics": 3,
        "picard_collect_multiple_metrics": 3,
        "trim_leading_attachment_site": 3,
    }

    config: Dict[str, Any] = {
        "settings": [
            {
                "name": "first",
                "fq_dir": f"{fq_dir}",
                "ref_fasta": f"{ref_fasta}",
                "attachment_sites": [
                    "attP:GTGGTTTGTCTGGTCAACCACCGCG:GT:CTCAGTGGTGTACGGTACAAACCCA"
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
                "attachment_sites": ["attB:GGCCGGCTTGTCGACGACGGCG:GT:CTCCGTCGTCAGGATCATCCGG"],
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

    config_yml: Path = Path(tmpdir) / "config.yml"
    with config_yml.open("w") as writer:
        yaml.dump(config, writer)
    with config_yml.open("r") as reader:
        for line in reader:
            print(line.rstrip())

    run_snakemake(
        pipeline="cryptic-seq",
        workdir=tmpdir,
        rules=rules,
        config={"config_yml": config_yml},
    )
