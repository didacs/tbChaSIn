################################################################################
# Pipeline for CRYPTIC Seq Analysis
################################################################################

from pathlib import Path
from typing import Dict
from typing import List
from typing import Optional

from pytomebio.pipeline import snakemake_utils
from pytomebio.pipeline.samples import from_path
from pytomebio.pipeline.samples import Sample
import os
import copy
from attrs import define


################################################################################
# Utility methods and variables
################################################################################

yml: Path = Path(config["config_yml"])
sample_dict: Dict[str, Sample] = from_path(yml=yml)
samples: List[Sample] = list(sample_dict.values())

################################################################################
# Terminal files
################################################################################

# File extensions to generate
extensions = [
    "alignment_summary_metrics.txt",
    "mark_duplicates.txt",
    "quality_yield_metrics.txt",  # for picard CollectMultipleMetrics
    "sites.txt",
]


all_terminal_files: List[Path] = [
    Path(f"{sample.group}/{sample.name}/{sample.name}.{ext}")
    for sample in samples
    for ext in extensions
]

all_terminal_files.append(Path("sequencing_quality_report.html"))
all_terminal_files.append(Path("sites.per_sample.txt"))


################################################################################
# Snakemake rules
################################################################################


# Block of code that gets called if the snakemake pipeline exits with an error.
onerror:
    snakemake_utils.on_error(snakefile=Path(__file__), config=config, log=Path(log))


rule all:
    input:
        all_terminal_files,


rule fgbio_fastq_to_bam:
    """Convert FASTQ to BAM and extract UMI information.

    Runs:
    - fgbio's FastqToBam
    """
    input:
        fq1=lambda wildcards: sample_dict[wildcards.sample].fq1,
        fq2=lambda wildcards: sample_dict[wildcards.sample].fq2,
    params:
        read_structure="11M+T +T",  # TODO: make user configurable
    output:
        bam="{group}/{sample}/{sample}.raw.bam",
    log:
        "logs/{group}/{sample}.fgbio_fastq_to_bam.log",
    resources:
        mem_gb=4,
        jvm_gb=4,
    shell:
        """
        fgbio \
          -Dsamjdk.use_async_io_read_samtools=true \
          -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx{resources.jvm_gb}g \
          FastqToBam \
          --input {input.fq1} {input.fq2} \
          --output {output.bam} \
          --read-structure {params.read_structure} \
          --sample {wildcards.sample} \
          --library {wildcards.sample} \
        &> {log}
        """


rule cryptic_seq_trim_for_tn5me:
    """Trims the start of R1s for the Tn5 mosaic end.

    Runs:
    - tomebio-tools cryptic-seq trim-for-tn5me
    """
    input:
        bam="{group}/{sample}/{sample}.raw.bam",
    output:
        keep="{group}/{sample}/{sample}.cryptic_seq.trim_for_tn5me.keep.bam",
        reject="{group}/{sample}/{sample}.cryptic_seq.trim_for_tn5me.reject.bam",
        metric_tsv="{group}/{sample}/{sample}.cryptic_seq.trim_for_tn5me.metrics.tsv",
    log:
        "logs/{group}/{sample}.cryptic_seq.trim_for_tn5me.log",
    resources:
        mem_gb=2,
    shell:
        """
        tomebio-tools cryptic-seq trim-for-tn5me \
            --in-bam {input.bam} \
            --keep-bam {output.keep} \
            --reject-bam {output.reject} \
            --out-metrics {output.metric_tsv} \
        &> {log}
        """


rule trim_leading_attachment_site:
    """Trims the start of R2s for the leading attachment site.

    Runs:
    - tomebio-tools cryptic-seq trim-leading-attachment-site
    """
    input:
        bam=(
            "{group}/{sample}/{sample}.cryptic_seq.trim_for_tn5me.keep.bam"
            if config.get("trim_Tn5", True)
            else "{group}/{sample}/{sample}.raw.bam"
        ),
    params:
        attachment_sites=lambda wildcards: sample_dict[wildcards.sample].attachment_sites,
    output:
        keep="{group}/{sample}/{sample}.trim_leading_attachment_site.keep.bam",
        reject="{group}/{sample}/{sample}.trim_leading_attachment_site.reject.bam",
        metric_tsv="{group}/{sample}/{sample}.trim_leading_attachment_site.metrics.tsv",
    log:
        "logs/{group}/{sample}.trim_leading_attachment_site.log",
    resources:
        mem_gb=2,
    shell:
        """
        tomebio-tools cryptic-seq trim-leading-attachment-site \
            --in-bam {input.bam} \
            --keep-bam {output.keep} \
            --reject-bam {output.reject} \
            --out-metrics {output.metric_tsv} \
            --attachment-site {params.attachment_sites} \
            --max-mismatches 4 \
        &> {log}
        """


# TODO: min score is set low
rule change_seq_trim_for_tn5me:
    """Trims the end of R2s for the Tn5 mosaic end.

    The change-seq tool trims both ends, but this should be fine.

    Runs:
    - tomebio-tools change-seq trim-for-tn5me
    """
    input:
        bam="{group}/{sample}/{sample}.trim_leading_attachment_site.keep.bam",
    output:
        bam="{group}/{sample}/{sample}.change_seq.trim_for_tn5me.bam",
        metric_tsv="{group}/{sample}/{sample}.change_seq.trim_for_tn5me.metrics.tsv",
    log:
        "logs/{group}/{sample}.change_seq.trim_for_tn5me.log",
    resources:
        mem_gb=2,
    shell:
        """
        tomebio-tools change-seq trim-for-tn5me \
            --in-bam {input.bam} \
            --out-bam {output.bam} \
            --out-metrics {output.metric_tsv} \
            --min-score 12 \
            --no-is-circularized \
        &> {log}
        """


rule align:
    """Aligns the reads in the unmapped BAM to the reference genome and coordinate sorts

    Important: runs with -Y so hard-clipping is not performed on supplementary alignments,
    but instead soft-clipping.  This is important to retain all bases for downstream analysis.

    Runs:
    - samtools fastq
    - bwa mem
    - fgbio ZipperBams
    - samtools sort
    """
    input:
        bam="{group}/{sample}/{sample}.change_seq.trim_for_tn5me.bam",
        ref_fasta=lambda wildcards: sample_dict[wildcards.sample].ref_fasta,
        bwa_files=lambda wildcards: [
            Path(f"{sample_dict[wildcards.sample].ref_fasta}.{ext}")
            for ext in ["amb", "ann", "bwt", "pac", "sa"]
        ],
    output:
        bam="{group}/{sample}/{sample}.mapped.bam",
        csi="{group}/{sample}/{sample}.mapped.bam.csi",
    log:
        "logs/{group}/{sample}.align.log",
    threads: 4
    resources:
        mem_gb=8,  # 8 for bwa, 4 for fgbio, and 16 for samtools sort (8 threads times 2G per thread)
        fgbio_mem_gb=4,
        sort_mem_gb=2,
        sort_threads=1,
    shell:
        """
        (samtools fastq -N {input.bam} \
          | bwa mem -Y -K 320000000 -t {threads} -p {input.ref_fasta} - \
          | fgbio -Xmx{resources.fgbio_mem_gb}g --compression 0 ZipperBams --unmapped {input.bam} --ref {input.ref_fasta} \
          | samtools sort -m {resources.sort_mem_gb}G -@ {resources.sort_threads} -o {output.bam} --write-index
        ) &> {log}
        """


rule fgbio_clip_bam:
    """Clip reads in FR pairs that sequence past the far end of their mate

    Runs:
    - fgbio ClipBam
    """
    input:
        bam="{group}/{sample}/{sample}.mapped.bam",
        ref_fasta=lambda wildcards: sample_dict[wildcards.sample].ref_fasta,
    output:
        bam="{group}/{sample}/{sample}.clipped.bam",
        txt="{group}/{sample}/{sample}.clipped.txt",
    log:
        "logs/{group}/{sample}.fgbio_clip_bam.log",
    resources:
        mem_gb=12,
        jvm_gb=8,
    shell:
        """
        samtools sort -n -u {input.bam} \
        | fgbio \
          -Dsamjdk.use_async_io_read_samtools=true \
          -Duse_async_io_write_samtools=true \
          -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx{resources.jvm_gb}g \
          ClipBam \
          --input /dev/stdin \
          --ref {input.ref_fasta} \
          --output {output.bam} \
          --metrics {output.txt} \
          --clipping-mode=Soft \
          --clip-bases-past-mate=true \
          --sort-order=Coordinate \
        &> {log}
        """


rule mark_duplicates:
    """Marks PCR duplicates

    Runs:
    - picard MarkDuplicates
    """
    input:
        bam="{group}/{sample}/{sample}.clipped.bam",
    output:
        bam="{group}/{sample}/{sample}.deduped.bam",
        txt="{group}/{sample}/{sample}.mark_duplicates.txt",
        tmp_dir=temp(directory("{group}/{sample}/tmp_mark_duplicates_dir")),
    log:
        "logs/{group}/{sample}.mark_duplicates.log",
    resources:
        mem_gb=12,
        jvm_gb=8,
    shell:
        """
        picard \
          -Dsamjdk.use_async_io_read_samtools=true \
          -Duse_async_io_write_samtools=true \
          -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx{resources.jvm_gb}g \
          MarkDuplicates \
          --INPUT {input.bam} \
          --OUTPUT {output.bam} \
          --METRICS_FILE {output.txt} \
          --CREATE_INDEX \
          --TMP_DIR {output.tmp_dir} \
          --BARCODE_TAG RX \
        &> {log}
        """


rule picard_collect_alignment_summary_metrics:
    """Collect alignment summary metrics.

    Note: this tools is not run as part of picard CollectMultipleMetrics as the plotting of the
    former tool fails on empty BAM files, and there's no way to skip generating plots in the
    latter tool.  Thus we run this tool on its own without plotting output (HISTOGRAM argument).

    Runs:
    - picard CollectAlignmentSummaryMetrics
    """
    input:
        bam="{group}/{sample}/{sample}.deduped.bam",
        ref_fasta=lambda wildcards: sample_dict[wildcards.sample].ref_fasta,
    output:
        txt="{group}/{sample}/{sample}.alignment_summary_metrics.txt",
    log:
        "logs/{group}/{sample}.picard_collect_alignment_summary_metrics.log",
    resources:
        mem_gb=12,
        jvm_gb=8,
    shell:
        """
        picard \
          -Dsamjdk.use_async_io_read_samtools=true \
          -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx{resources.jvm_gb}g \
          CollectAlignmentSummaryMetrics \
          --INPUT {input.bam} \
          --OUTPUT {output.txt} \
          --REFERENCE_SEQUENCE {input.ref_fasta} \
         &> {log}
        """


rule picard_collect_multiple_metrics:
    """Collect basic sequencing quality metrics.

    Runs:
    - picard CollectMultipleMetrics
    """
    input:
        bam="{group}/{sample}/{sample}.deduped.bam",
        ref_fasta=lambda wildcards: sample_dict[wildcards.sample].ref_fasta,
    params:
        prefix="{group}/{sample}/{sample}",
    output:
        "{group}/{sample}/{sample}.quality_yield_metrics.txt",
    log:
        "logs/{group}/{sample}.picard_collect_multiple_metrics.log",
    resources:
        mem_gb=12,
        jvm_gb=8,
    shell:
        """
        picard \
          -Dsamjdk.use_async_io_read_samtools=true \
          -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx{resources.jvm_gb}g \
          CollectMultipleMetrics \
          --INPUT {input.bam} \
          --OUTPUT {params.prefix} \
          --REFERENCE_SEQUENCE {input.ref_fasta} \
          --PROGRAM null \
          --PROGRAM QualityScoreDistribution \
          --PROGRAM MeanQualityByCycle \
          --PROGRAM CollectBaseDistributionByCycle \
          --PROGRAM CollectGcBiasMetrics \
          --PROGRAM CollectSequencingArtifactMetrics \
          --PROGRAM CollectInsertSizeMetrics \
          --PROGRAM CollectQualityYieldMetrics \
          --FILE_EXTENSION .txt \
         &> {log}
        """


rule fastqc:
    """Collect pre-alignment sequencing quality control

    Note: run as a single thread since FASTQC only parallelizes
    across files, so adding threading doesn't help here.

    Runs:
    - fastqc
    """
    input:
        fq=lambda wildcards: (
            sample_dict[wildcards.sample].fq1
            if wildcards.read_num == "1"
            else sample_dict[wildcards.sample].fq2
        ),
    params:
        outdir="{group}/{sample}",
    output:
        html="{group}/{sample}/{sample}.R{read_num}_fastqc.html",
        zip="{group}/{sample}/{sample}.R{read_num}_fastqc.zip",
    log:
        "logs/{group}/{sample}.fastqc.R{read_num}.log",
    shell:
        """
        (mkdir -p {params.outdir} &&
          (cat {input.fq} \
            | zcat -fc \
            | fastqc \
              --outdir {params.outdir} \
              stdin:{wildcards.sample}.R{wildcards.read_num} \
          )
        ) &> {log}
        """


# TODO: Customize the MultiQC config?
rule multiqc:
    """Aggregates metrics into a HTML (MultiQC) report

    Runs:
    - MultiQC
    """
    input:
        fastqc_html=[
            f"{sample.group}/{sample.name}/{sample.name}.R{r}_fastqc.html"
            for sample in samples
            for r in [1, 2]
        ],
        fastqc_zip=[
            f"{sample.group}/{sample.name}/{sample.name}.R{r}_fastqc.zip"
            for sample in samples
            for r in [1, 2]
        ],
        dupe=[
            f"{sample.group}/{sample.name}/{sample.name}.mark_duplicates.txt"
            for sample in samples
        ],
        asm=[
            f"{sample.group}/{sample.name}/{sample.name}.alignment_summary_metrics.txt"
            for sample in samples
        ],
        qym=[
            f"{sample.group}/{sample.name}/{sample.name}.quality_yield_metrics.txt"
            for sample in samples
        ],
    params:
        directory=f"{os.getcwd()}",
    output:
        multiqc_report="sequencing_quality_report.html",
    log:
        f"logs/multiqc.log",
    shell:
        "multiqc {params.directory} --force --no-ansi --filename {output.multiqc_report} &> {log}"


rule find_sites:
    """Finds putative integration sties.

    Runs:
    - tomebio-tools cryptic-seq find-sites
    """
    input:
        bam="{group}/{sample}/{sample}.deduped.bam",
    output:
        tsv="{group}/{sample}/{sample}.sites.txt",
    log:
        "logs/{group}/{sample}.find_sites.log",
    resources:
        mem_gb=2,
    shell:
        """
        samtools sort \
         -n \
         {input.bam} \
         | tomebio-tools cryptic-seq find-sites \
         --in-bam - \
         --out-txt {output.tsv} \
        &> {log}
        """


rule collate_sites:
    """Collates tables for identified integration sites across samples and replicates.

    Runs:
    - tomebio-tools common collate-sites
    """
    input:
        yml=yml,
        txt=[f"{sample.group}/{sample.name}/{sample.name}.sites.txt" for sample in samples],
    params:
        out_pre="sites",
    output:
        tsv="sites.per_sample.txt",
    log:
        "logs/collate_sites.log",
    resources:
        mem_gb=2,
    shell:
        """
        tomebio-tools common collate-sites \
         --in-yml {input.yml} \
         --out-prefix {params.out_pre} \
        &> {log}
        """
