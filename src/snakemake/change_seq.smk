################################################################################
# Pipeline for CHANGE Seq Analysis
################################################################################

from pathlib import Path
from typing import Dict
from typing import List
from typing import Optional

from pytomebio.pipeline import snakemake_utils
from pytomebio.pipeline.samples import from_config
from pytomebio.pipeline.samples import Sample
import os
import copy
from attrs import define


################################################################################
# Utility methods and variables
################################################################################

sample_dict: Dict[str, Sample] = from_config(config=config)
samples: List[Sample] = list(sample_dict.values())

################################################################################
# Terminal files
################################################################################

# File extensions to generate
extensions = [
    "alignment_summary_metrics.txt",
    "mark_duplicates.txt",
    "quality_yield_metrics.txt",  # for picard CollectMultipleMetrics
    "svpileup.aggregated.txt",
    "sites.txt"
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

onerror:
    """Block of code that gets called if the snakemake pipeline exits with an error."""
    snakemake_utils.on_error(snakefile=Path(__file__), config=config, log=Path(log))


rule all:
    input:
        all_terminal_files,


rule fastq_to_bam:
    """Converts the input FASTQs to an unmapped BAM
    
    Runs:
    - tomebio-tools change-seq fastq-to-bam
    """
    input:
        fq1 = lambda wildcards: sample_dict[wildcards.sample].fq1,
        fq2 = lambda wildcards: sample_dict[wildcards.sample].fq2,
    output:
        keep_bam   = "{group}/{sample}/{sample}.fastq_to_bam.keep.bam",
        reject_bam = "{group}/{sample}/{sample}.fastq_to_bam.reject.bam",
        metric_tsv = "{group}/{sample}/{sample}.fastq_to_bam.metrics.tsv"
    log: "logs/{group}/{sample}.fastq_to_bam.log"
    resources:
        mem_gb = 2
    shell:
        """
        tomebio-tools change-seq fastq-to-bam \
            --r1 {input.fq1} \
            --r2 {input.fq2} \
            --keep-bam {output.keep_bam} \
            --reject-bam {output.reject_bam} \
            --metric-tsv {output.metric_tsv} \
            --read-group {wildcards.sample} \
            --threads 4 \
        &> {log}
        """


rule trim_for_tn5me:
    """Trims the ends of reads for the Tn5 mosaic end and palindromic sequence

    Runs:
    - tomebio-tools change-seq trim-for-tn5me
    """
    input:
        bam="{group}/{sample}/{sample}.fastq_to_bam.keep.bam",
    output:
        bam="{group}/{sample}/{sample}.trim_for_tn5me.bam",
        metric_tsv="{group}/{sample}/{sample}.trim_for_tn5me.metrics.tsv"
    log: "logs/{group}/{sample}.trim_for_tn5me.log"
    resources:
        mem_gb=2
    shell:
        """
        tomebio-tools change-seq trim-for-tn5me \
            --in-bam {input.bam} \
            --out-bam {output.bam} \
            --out-metrics {output.metric_tsv} \
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
        bam = "{group}/{sample}/{sample}.trim_for_tn5me.bam",
        ref_fasta = lambda wildcards: sample_dict[wildcards.sample].ref_fasta,
        bwa_files = lambda wildcards: [Path(f"{sample_dict[wildcards.sample].ref_fasta}.{ext}") for ext in ["amb", "ann", "bwt", "pac", "sa"]]
    output:
        bam = "{group}/{sample}/{sample}.mapped.bam",
        csi = "{group}/{sample}/{sample}.mapped.bam.csi"
    log: "logs/{group}/{sample}.align.log"
    threads: 18
    resources:
        mem_gb = 32,  # 8 for bwa, 4 for fgbio, and 16 for samtools sort (8 threads times 2G per thread)
        fgbio_mem_gb = 4,
        sort_mem_gb = 2,
        sort_threads = 8
    shell:
        """
        (samtools fastq -N {input.bam} \
          | bwa mem -Y -K 320000000 -t {threads} -p {input.ref_fasta} - \
          | fgbio -Xmx{resources.fgbio_mem_gb}g --compression 0 ZipperBams --unmapped {input.bam} --ref {input.ref_fasta} \
          | samtools sort -m {resources.sort_mem_gb}G -@ {resources.sort_threads} -o {output.bam} --write-index
        ) &> {log}
        """


# TODO: use --BARCODE_TAG RX if UMIs are present
rule mark_duplicates:
    """Marks PCR duplicates
    
    Runs:
    - picard MarkDuplicates
    """
    input:
        bam = "{group}/{sample}/{sample}.mapped.bam"
    output:
        bam = "{group}/{sample}/{sample}.deduped.bam",
        txt = "{group}/{sample}/{sample}.mark_duplicates.txt",
        tmp_dir = temp(directory("{group}/{sample}/tmp_mark_duplicates_dir"))
    log: "logs/{group}/{sample}.mark_duplicates.log"
    resources:
        mem_gb=12,
        jvm_gb=8
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
        &> {log}
        """


rule picard_collect_alignment_summary_metrics:
    """Collect alignment summary  metrics.
    
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
        txt = "{group}/{sample}/{sample}.alignment_summary_metrics.txt"
    log:
        "logs/{group}/{sample}.picard_collect_alignment_summary_metrics.log"
    resources:
        mem_gb=12,
        jvm_gb=8
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
        ref_fasta = lambda wildcards: sample_dict[wildcards.sample].ref_fasta,
    output:
        "{group}/{sample}/{sample}.quality_yield_metrics.txt"
    params:
        prefix="{group}/{sample}/{sample}"
    log:
        "logs/{group}/{sample}.picard_collect_multiple_metrics.log"
    resources:
        mem_gb=12,
        jvm_gb=8
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
        fq = lambda wildcards: sample_dict[wildcards.sample].fq_dir / "{sample}_R{read_num}_001.fastq.gz"
    output:
        html = "{group}/{sample}/{sample}.{read_num}_fastqc.html",
        zip = "{group}/{sample}/{sample}.{read_num}_fastqc.zip"
    params:
        outdir = "{group}/{sample}"
    log:
        "logs/{group}/{sample}.fastqc.R{read_num}.log"
    shell:
        """
        (mkdir -p {params.outdir} &&
          (cat {input.fq} \
            | zcat -fc \
            | fastqc \
              --outdir {params.outdir} \
              stdin:{wildcards.sample}.{wildcards.read_num} \
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
        fastqc_html = [f"{sample.group}/{sample.name}/{sample.name}.{r}_fastqc.html" for sample in samples for r in [1, 2]],
        fastqc_zip = [f"{sample.group}/{sample.name}/{sample.name}.{r}_fastqc.zip" for sample in samples for r in [1, 2]],
        dupe = [f"{sample.group}/{sample.name}/{sample.name}.mark_duplicates.txt" for sample in samples],
        asm = [f"{sample.group}/{sample.name}/{sample.name}.alignment_summary_metrics.txt" for sample in samples],
        qym = [f"{sample.group}/{sample.name}/{sample.name}.quality_yield_metrics.txt" for sample in samples],
    output:
        multiqc_report = "sequencing_quality_report.html",
    params:
        directory=f"{os.getcwd()}"
    log:
        f"logs/multiqc.log"
    shell:
        "multiqc {params.directory} --force --no-ansi --filename {output.multiqc_report} &> {log}"


rule fgsv_svpileup:
    """Collates a pileup of putative structural variant supporting reads.

    Runs:
    - fgsv SvPileup
    """
    input:
        bam = "{group}/{sample}/{sample}.deduped.bam"
    output:
        txt = "{group}/{sample}/{sample}.svpileup.txt",
        bam = "{group}/{sample}/{sample}.svpileup.bam"
    params:
        prefix = "{group}/{sample}/{sample}.svpileup",
    log: "logs/{group}/{sample}.fgsv_svpileup.log"
    resources:
        mem_gb=5,
        jvm_gb=4
    shell:
        """
        fgsv \
          -Dsamjdk.use_async_io_read_samtools=true \
          -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx{resources.jvm_gb}g \
          SvPileup \
          --input={input.bam} \
          --output={params.prefix} \
        &> {log}
        """


rule fgsv_aggregatesvpileup:
    """Merges nearby pileups of reads supporting putative breakpoints.

    Runs:
    - fgsv AggregateSvPileup
    """
    input:
        txt = "{group}/{sample}/{sample}.svpileup.txt"
    output:
        txt = "{group}/{sample}/{sample}.svpileup.aggregated.txt",
    log: "logs/{group}/{sample}.fgsv_aggregatesvpileup.log"
    resources:
        mem_gb=5,
        jvm_gb=4
    shell:
        """
        fgsv \
          -Dsamjdk.use_async_io_read_samtools=true \
          -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx{resources.jvm_gb}g \
          AggregateSvPileup \
          --input={input.txt} \
          --output={output.txt} \
        &> {log}
        """


rule find_sites:
    """Finds putative integration sties.

    Runs:
    - tomebio-tools change-seq find-sites
    """
    input:
        bam="{group}/{sample}/{sample}.deduped.bam",
    output:
        tsv="{group}/{sample}/{sample}.sites.txt"
    log: "logs/{group}/{sample}.find_sites.log"
    resources:
        mem_gb=2
    shell:
        """
        samtools sort \
         -n \
         {input.bam} \
         | tomebio-tools change-seq find-sites \
         --in-bam - \
         --out-txt {output.tsv} \
        &> {log}
        """


rule collate_sites:
    """Collates tables for identified integration sites across samples and replicates.

    Runs:
    - tomebio-tools change-seq collate-sites
    """
    input:
        txt=[f"{sample.group}/{sample.name}/{sample.name}.sites.txt" for sample in samples]
    output:
        tsv="sites.per_sample.txt"
    params:
        out_pre="sites"
    log: "logs/collate_sites.log"
    resources:
        mem_gb=2
    shell:
        """
        tomebio-tools change-seq collate-sites \
         --in-txt {input.txt} \
         --out-prefix {params.out_pre} \
        &> {log}
        """