import Meta
import References

include {
    validateParameters;
    paramsHelp;
    paramsSummaryLog;
} from 'plugin/nf-schema'

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
    log.info paramsHelp("nextflow run cryptic-seq --metasheet metasheet.txt")
    exit 0
}

// Validate input parameters
validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

// Collect pre-alignment sequencing quality control
// Note: run as a single thread since FASTQC only parallelizes across files,
// so adding threading doesn't help here.
process fastqc {
    label 'fastqc'
    tag "${meta.uniqueId}"
    ext memory: '4.GB'

    input:
        tuple val(meta), val(read_num), path(fq)
    
    output:
        path html, emit: report_html
        path zip, emit: report_zip
    
    script:
        def prefix = "${meta.uniqueId}.R${read_num}"
        html = "${prefix}_fastqc.html"
        zip = "${prefix}_fastqc.zip"
        """
        zcat -f ${fq} | fastqc "stdin:${prefix}"
        """
}

// Convert FASTQ to unaligned BAM and extract UMI information.
process fastq_to_ubam {
    label 'fgbio'
    tag "${meta.uniqueId}"

    input:
        tuple val(meta), path(fq1), path(fq2)

    output:
        tuple val(meta), path(ubam), emit: ubam
    
    script:
        ubam = "${meta.uniqueId}.raw.bam"
        """
        fgbio \
            -Dsamjdk.use_async_io_read_samtools=true \
            -XX:GCTimeLimit=50 \
            -XX:GCHeapFreeLimit=10 \
            -Xmx${task.memory.giga}g \
            FastqToBam \
            --input "${fq1}" "${fq2}" \
            --output "${ubam}" \
            --read-structure "${params.read_structure_r1}" "${params.read_structure_r2}" \
            --sample "${meta.uniqueId}" \
            --library "${meta.uniqueId}"
        """
}

// Trims the start of R1s for the Tn5 mosaic end.
process trim_r1_tn5me {
    label 'tools'
    tag "${meta.uniqueId}"

    input:
        tuple val(meta), path(bam)
    
    output:
        tuple val(meta), path(keep_bam), emit: keep_bam
        tuple val(meta), path(reject_bam), emit: reject_bam
        path metric_tsv, emit: metric_tsv
    
    script:
        keep_bam = "${meta.uniqueId}.cryptic_seq.trim_for_tn5me.keep.bam"
        reject_bam = "${meta.uniqueId}.cryptic_seq.trim_for_tn5me.reject.bam"
        metric_tsv = "${meta.uniqueId}.cryptic_seq.trim_for_tn5me.metrics.tsv"
        """
        tomebio-tools cryptic-seq trim-r1-tn5me \
            --in-bam "${bam}" \
            --keep-bam "${keep_bam}" \
            --reject-bam "${reject_bam}" \
            --out-metrics "${metric_tsv}" \
            --max-mismatches ${params.trim_Tn5_max_mismatches}
        """
}

// Trims the start of R2s for the leading attachment site.
process trim_leading_attachment_site {
    label 'tools'
    tag "${meta.uniqueId}"

    input:
        tuple val(meta), path(bam)
    
    output:
        tuple val(meta), path(keep_bam), emit: keep_bam
        tuple val(meta), path(reject_bam), emit: reject_bam
        path metric_tsv, emit: metric_tsv
    
    script:
        keep_bam = "${meta.uniqueId}.trim_leading_attachment_site.keep.bam"
        reject_bam = "${meta.uniqueId}.trim_leading_attachment_site.reject.bam"
        metric_tsv = "${meta.uniqueId}.trim_leading_attachment_site.metrics.tsv"
        """
        tomebio-tools cryptic-seq trim-leading-attachment-site \
            --in-bam "${bam}" \
            --keep-bam "${keep_bam}" \
            --reject-bam "${reject_bam}" \
            --out-metrics "${metric_tsv}" \
            --attachment-site ${meta.attachmentSites.join(' ')} \
            --max-mismatches ${params.trim_att_max_mismatches}
        """
}

// Trims the end of R2s for the Tn5 mosaic end.
process trim_r2_tn5me {
    label 'tools'
    tag "${meta.uniqueId}"

    input:
        tuple val(meta), path(bam)
    
    output:
        tuple val(meta), path(keep_bam), emit: keep_bam
        tuple val(meta), path(reject_bam), emit: reject_bam
        path metric_tsv, emit: metric_tsv
    
    script:
        keep_bam = "${meta.uniqueId}.cryptic_seq.trim_r2_for_tn5me.keep.bam"
        reject_bam = "${meta.uniqueId}.cryptic_seq.trim_r2_for_tn5me.reject.bam"
        metric_tsv = "${meta.uniqueId}.cryptic_seq.trim_r2_for_tn5me.metrics.tsv"
        def min_score_opt = ""
        if (params.trim_Tn5_r2_min_score != null) {
            min_score_opt = "--min-score ${params.trim_Tn5_r2_min_score}"
        }
        def min_match_length_opt = ""
        if (params.trim_Tn5_r2_min_match_length) {
            min_match_length_opt = "--min-match-length ${params.trim_Tn5_r2_min_match_length}"
        }
        def min_keep_length_opt = ""
        if (params.trim_Tn5_r2_min_keep_length) {
            min_keep_length_opt = "--min-keep-length ${params.trim_Tn5_r2_min_keep_length}"
        }
        """
        tomebio-tools cryptic-seq trim-r2 \
            --in-bam "${bam}" \
            --keep-bam "${keep_bam}" \
            --reject-bam "${reject_bam}" \
            --out-metrics "${metric_tsv}" \
            ${min_score_opt} \
            ${min_match_length_opt} \
            ${min_keep_length_opt}
        """
}

// Trims the end of R2s for the Tn5 mosaic end.
process trim_r2_adapter {
    label 'tools'
    tag "${meta.uniqueId}"

    input:
        tuple val(meta), path(bam)
    
    output:
        tuple val(meta), path(keep_bam), emit: keep_bam
        tuple val(meta), path(reject_bam), emit: reject_bam
        path metric_tsv, emit: metric_tsv
    
    script:
        keep_bam = "${meta.uniqueId}.cryptic_seq.trim_r2_for_adapter.keep.bam"
        reject_bam = "${meta.uniqueId}.cryptic_seq.trim_r2_for_adapter.reject.bam"
        metric_tsv = "${meta.uniqueId}.cryptic_seq.trim_r2_for_adapter.metrics.tsv"
        def min_score_opt = ""
        if (params.trim_adapter_r2_min_score != null) {
            min_score_opt = "--min-score ${params.trim_adapter_r2_min_score}"
        }
        def min_match_length_opt = ""
        if (params.trim_adapter_r2_min_match_length) {
            min_match_length_opt = "--min-match-length ${params.trim_adapater_r2_min_match_length}"
        }
        def min_keep_length_opt = ""
        if (params.trim_adapter_r2_min_keep_length) {
            min_keep_length_opt = "--min-keep-length ${params.trim_adapter_r2_min_keep_length}"
        }
        """
        tomebio-tools cryptic-seq trim-r2 \
            --in-bam "${bam}" \
            --keep-bam "${keep_bam}" \
            --reject-bam "${reject_bam}" \
            --out-metrics "${metric_tsv}" \
            --sequence ${params.trim_adapter_r2} \
            ${min_score_opt} \
            ${min_match_length_opt} \
            ${min_keep_length_opt}
        """
}

// Concatenates multiple BAMs into a single BAM.
process concatenate_bams {
    label 'fgbio' // TODO: create a separate hts label/image
    tag "${meta.uniqueId}"
    ext cpus: 4

    input:
        tuple val(meta), path(bams)
    
    output:
        tuple val(meta), path(concatenated_bam), emit: concatenated_bam
    
    shell:
        concatenated_bam = "${meta.uniqueId}.concatenated.bam"
        def bam_paths = bams.collect { it.toString() }.join(' ')
        def compression_threads = task.cpus - 1
        """
        readlink ${bam_paths} > bams
        samtools cat \
            -b bams \
            -o "${concatenated_bam}" \
            -@ ${compression_threads}
        """
}

// Aligns the reads in the unmapped BAM to the reference genome and sorts 
// by coordinate. Important: runs with -Y so hard-clipping is not performed 
// on supplementary alignments, but instead soft-clipping. This is important
// to retain all bases for downstream analysis.
process align {
    label 'fgbio'
    tag "${meta.uniqueId}"
    ext cpus: 8, memory: '16.GB'

    input:
        // `reference_dir` is the directory that contains the reference FASTA 
        // (bgzipped and indexed) and the BWA index. The directory name is expected
        // to be the same as the reference prefix, e.g. `/path/to/GRCh38/GRCh38.fasta.gz`.
        tuple val(meta), path(bam), path(reference_dir), val(reference_fasta_name)
    
    output:
        tuple val(meta), path(mapped_bam), path(mapped_index), emit: mapped_bam

    shell:
        mapped_bam = "${meta.uniqueId}.mapped.bam"
        mapped_index = "${mapped_bam}.bai"
        def reference_fasta = reference_dir / "${reference_fasta_name}"
        // TODO: make these configurable
        def sort_threads = 2
        def total_mem_gb = task.memory.giga
        def sort_mem_gb = [Math.round(total_mem_gb * 0.1), 2].max() / sort_threads
        def fgbio_mem_gb = [Math.round(total_mem_gb * 0.1), 1].max()
        """
        samtools fastq -N "${bam}" \
        | bwa mem \
            -Y \
            -K 320000000 \
            -t ${task.cpus} \
            -p "${reference_fasta}" \
            - \
        | fgbio \
            -Xmx${fgbio_mem_gb}g \
            --compression 0 \
            ZipperBams \
            --unmapped "${bam}" \
            --ref "${reference_fasta}" \
        | samtools sort \
            -m ${sort_mem_gb}G \
            -@ ${sort_threads} \
            -u \
            -O bam \
        | samtools view \
            -@ ${task.cpus} \
            -o "${mapped_bam}##idx##${mapped_index}" \
            --write-index \
            -
        """
}

// Clip reads in FR pairs that sequence past the far end of their mate.
process clip_bam {
    label 'fgbio'
    tag "${meta.uniqueId}"
    ext cpus: 4, memory: "5.GB"

    input:
        tuple val(meta), path(bam), path(index), path(reference_fasta), path(reference_index)
    
    output:
        // TODO: verify that ClipBam writes an index
        tuple val(meta), path(clipped_bam), path(clipped_index), emit: clipped_bam
        path metric_tsv, emit: metric_tsv
    
    script:
        clipped_bam = "${meta.uniqueId}.clipped.bam"
        clipped_index = "${meta.uniqueId}.clipped.bai"
        metric_tsv = "${meta.uniqueId}.clipped.txt"
        // The sorting step will complete before ClipBam runs, so we can
        // use all CPUs and memory for sorting
        def sort_mem_gb = [Math.round((task.memory.giga - 1) / task.cpus), 1].max()
        """
        set -eo pipefail
        samtools sort \
            -n \
            -m ${sort_mem_gb}G \
            -@ ${task.cpus} \
            -u \
            -O bam \
            "${bam}" \
        | fgbio \
            -Dsamjdk.use_async_io_read_samtools=true \
            -Duse_async_io_write_samtools=true \
            -XX:GCTimeLimit=50 \
            -XX:GCHeapFreeLimit=10 \
            -Xmx${task.memory.giga}g \
            ClipBam \
            --input /dev/stdin \
            --ref "${reference_fasta}" \
            --output "${clipped_bam}" \
            --metrics "${metric_tsv}" \
            --clipping-mode=Soft \
            --clip-bases-past-mate=true \
            --sort-order=Coordinate
        """
}

process copy_umi_from_read_name {
    label 'fgbio'
    tag "${meta.uniqueId}"

    input:
        tuple val(meta), path(bam), path(index)
    
    output:
        tuple val(meta), path(bam_with_rx), path(index_with_rx), emit: bam_with_rx

    script:
        bam_with_rx = "${meta.uniqueId}.umi_from_read_name.bam"
        index_with_rx = "${meta.uniqueId}.umi_from_read_name.bai"
        def total_mem_gb = task.memory.giga
        // TODO: make configurable
        def fgbio_mem_gb = total_mem_gb >= 10 ? 
            Math.round(total_mem_gb * 0.8) : [total_mem_gb - 1, 1].max()
        """
        fgbio \
            -Dsamjdk.use_async_io_read_samtools=true \
            -XX:GCTimeLimit=50 \
            -XX:GCHeapFreeLimit=10 \
            -Xmx${fgbio_mem_gb}g \
            CopyUmiFromReadName \
            --input "${bam}" \
            --output "${bam_with_rx}"
        """
}

// Marks PCR duplicates.
process mark_duplicates {
    label 'picard'
    tag "${meta.uniqueId}"
    ext memory: '16.GB'

    input:
        tuple val(meta), path(bam), path(index)
    
    output:
        // TODO: the Snakemake rule outputs the temp dir - is that necessary?
        tuple val(meta), path(deduped_bam), path(deduped_index), emit: deduped_bam
        path metric_txt, emit: metric_txt
    
    script:
        deduped_bam = "${meta.uniqueId}.deduped.bam"
        deduped_index = "${meta.uniqueId}.deduped.bai"
        metric_txt = "${meta.uniqueId}.mark_duplicates.txt"
        def total_mem_gb = task.memory.giga
        // TODO: make configurable
        def picard_mem_gb = total_mem_gb >= 10 ?
            Math.round(total_mem_gb * 0.8) : [total_mem_gb - 1, 1].max()
        // TODO: Adding this to the command line causes an error on MacOS
        // -Djdk.lang.Process.launchMechanism=vfork
        """
        picard \
            -Dsamjdk.use_async_io_read_samtools=true \
            -Duse_async_io_write_samtools=true \
            -XX:GCTimeLimit=50 \
            -XX:GCHeapFreeLimit=10 \
            -Xmx${picard_mem_gb}g \
            MarkDuplicates \
            --INPUT "${bam}" \
            --OUTPUT "${deduped_bam}" \
            --METRICS_FILE "${metric_txt}" \
            --CREATE_INDEX \
            --BARCODE_TAG RX
        """
}

// Collect alignment summary metrics.
// Note: this tools is not run as part of picard CollectMultipleMetrics
// as the plotting of the former tool fails on empty BAM files, and
// there's no way to skip generating plots in the latter tool. Thus we
// run this tool on its own without plotting output (HISTOGRAM argument).
process collect_alignment_summary_metrics {
    label 'picard'
    tag "${meta.uniqueId}"
    ext memory: '16.GB'

    input:
        tuple val(meta), path(bam), path(index), path(reference_fasta), path(reference_index)
    
    output:
        path metric_txt, emit: metric_txt 
    
    script:
        metric_txt = "${meta.uniqueId}.alignment_summary_metrics.txt"
        def total_mem_gb = task.memory.giga
        // TODO: make configurable
        def picard_mem_gb = total_mem_gb >= 10 ?
            Math.round(total_mem_gb * 0.8) : [total_mem_gb - 1, 1].max()
        // TODO: Adding this to the command line causes an error on MacOS
        // -Djdk.lang.Process.launchMechanism=vfork
        """
        picard \
            -Dsamjdk.use_async_io_read_samtools=true \
            -XX:GCTimeLimit=50 \
            -XX:GCHeapFreeLimit=10 \
            -Xmx${picard_mem_gb}g \
            CollectAlignmentSummaryMetrics \
            --INPUT "${bam}" \
            --OUTPUT "${metric_txt}" \
            --REFERENCE_SEQUENCE "${reference_fasta}"
        """
}

// Collect basic sequencing quality metrics.
process collect_multiple_metrics {
    label 'picard'
    tag "${meta.uniqueId}"
    ext memory: '16.GB'

    input:
        tuple val(meta), path(bam), path(index), path(reference_fasta), path(reference_index)
    
    output:
        path metric_txt, emit: metric_txt 
    
    script:
        metric_txt = "${meta.uniqueId}.quality_yield_metrics.txt"
        def total_mem_gb = task.memory.giga
        // TODO: make configurable
        def picard_mem_gb = total_mem_gb >= 10 ?
            Math.round(total_mem_gb * 0.8) : [total_mem_gb - 1, 1].max()
        // TODO: adding this to the command line causes an error on MacOS
        // -Djdk.lang.Process.launchMechanism=vfork
        """
        picard \
            -Dsamjdk.use_async_io_read_samtools=true \
            -XX:GCTimeLimit=50 \
            -XX:GCHeapFreeLimit=10 \
            -Xmx${picard_mem_gb}g \
            CollectMultipleMetrics \
            --INPUT "${bam}" \
            --OUTPUT "${meta.uniqueId}" \
            --REFERENCE_SEQUENCE "${reference_fasta}" \
            --PROGRAM null \
            --PROGRAM QualityScoreDistribution \
            --PROGRAM MeanQualityByCycle \
            --PROGRAM CollectBaseDistributionByCycle \
            --PROGRAM CollectGcBiasMetrics \
            --PROGRAM CollectSequencingArtifactMetrics \
            --PROGRAM CollectInsertSizeMetrics \
            --PROGRAM CollectQualityYieldMetrics \
            --FILE_EXTENSION .txt
        """
}

// Aggregates metrics into a HTML (MultiQC) report.
process multiqc {
    label 'multiqc'
    label 'global'
    ext memory: '4.GB'

    input:
        // Channel with a single list of all input files
        path files

    output:
        path report, emit: report
    
    script:
        report = "sequencing_quality_report.html"
        """
        multiqc . --force --no-ansi --filename "${report}"
        """
}

// Finds putative integration sites.
process find_sites {
    label 'tools'
    tag "${meta.uniqueId}"
    ext cpus: 4

    input:
        tuple val(meta), path(bam)
    
    output:
        tuple val(meta), path(sites_txt), emit: sites_txt
    
    script:
        sites_txt = "${meta.uniqueId}.sites.txt"
        // The sorting step will complete before find-sites runs, so we can
        // use all CPUs and memory for sorting
        def sort_mem_gb = [Math.round((task.memory.giga - 1) / task.cpus), 1].max()
        """
        samtools sort \
            -n \
            -m ${sort_mem_gb}G \
            -@ ${task.cpus} \
            -u \
            -O bam \
            "${bam}" \
        | tomebio-tools cryptic-seq find-sites \
            --in-bam - \
            --out-txt "${sites_txt}"
        """
}

// Collates tables for identified integration sites across samples and replicates.
process collate_sites {
    label 'tools'
    label 'global'

    input:
        // list of all sample metadata objects
        val meta
        // list of all sites files
        path sites_txt
    
    output:
        path collated_sites_txt, emit: collated_sites_txt
    
    script:
        def meta_json = Meta.toCollateJson(meta)
        def prefix = "sites"
        collated_sites_txt = "${prefix}.per_sample.txt"
        """
        printf '${meta_json}' > sample_metadata.json
        tomebio-tools common collate-sites \
            --in-json sample_metadata.json \
            --file-grouping flat \
            --out-prefix ${prefix}
        """
}

// Annotates cryptic sites and compiles output tables into a multi-sheet excel file
process annotate_sites {
    label 'tools'
    label 'global'
    ext memory: '8.GB'

    input:
        // list of all sample metadata objects
        val meta
        // list of all sites files
        path sites_txt
        // list of all unique reference fasta
        path reference_fasta, stageAs: 'references/*'
        // list of all unique GTFs
        path gene_set_gtfs, stageAs: 'references/*'
        // the name of the fasta to use
        tuple path(match_reference_dir), val(match_fasta_name)
    
    output:
        path annotated_sites_xlsx, emit: annotated_sites_xlsx

    script:
        def meta_json = Meta.toAnnotateJson(meta)
        def match_reference = "${match_reference_dir}/${match_fasta_name}"
        annotated_sites_xlsx = "sites.annotated.xlsx" 
        """
        printf '${meta_json}' > sample_metadata.json
        tomebio-tools cryptic-seq annotate-sites \
            --in-json sample_metadata.json \
            --ref-dir references \
            --match-index-base ${match_reference} \
            --file-grouping flat \
            --output ${annotated_sites_xlsx}
        """
}

// Native Nextflow workflow.
workflow cryptic_seq {
    prepare_metasheet()
    
    attachment_sites = null
    if (params.attachment_sites != null) {
        attachment_sites = file(params.attachment_sites)
    }
    default_attachment_site = null
    if (params.default_attachment_site != null) {
        default_attachment_site = params.default_attachment_site.split(";")
    }

    // Parse the references. We also need to pass in the name of the'
    // annotation reference, so we need to first concatenate the channels.
    references = prepare_metasheet.out.references
        .concat(prepare_metasheet.out.annotation_reference)
        .toList()
        .map { references, annotation_reference ->
            References.load(references, annotation_reference)
        }
    
    // Parse the metasheet into one `Meta` object per sample. We also
    // need to pass in the file with the references mapping, so we need
    // to first concatenate the channels.
    sample_meta = prepare_metasheet.out.metasheet
        .concat(references)
        .toList()
        .flatMap { metasheet, references ->
            Meta.fromTxt(
                metasheet,
                references,
                attachment_sites,
                default_attachment_site,
                file(params.fastq_dir),
                params.fastq_name_pattern,
                params.verify_files
            )
        }

    // call fastqc with each fastq file separately
    single_input = sample_meta.flatMap { meta ->
        [
            tuple(meta, 1, meta.fq1),
            tuple(meta, 2, meta.fq2)
        ]
    }
    fastqc(single_input)

    split_fastq = params.fastq_chunk_size != null

    // call the main pipeline with pairs of fastqs
    if (split_fastq) {
        // split each file into chunks of the specified size
        // right now the index is only used to make the intermediate file names unique,
        // but it could also be used to sort the chunks prior to merging so that the
        // merged reads are in the same order as the input
        idx = 1
        paired_input = sample_meta
            .map { meta -> tuple(meta, meta.fq1, meta.fq2) }
            .splitFastq(by: params.fastq_chunk_size, file: true, pe: true)
            .map { meta, fq1, fq2 -> tuple(meta.withIndex(idx++), fq1, fq2) }
    } else {
        paired_input = sample_meta
            .map { meta -> tuple(meta, meta.fq1, meta.fq2) }
    }
    
    fastq_to_ubam(paired_input)

    if (params.trim_Tn5) {
        trim_r1_tn5me(fastq_to_ubam.out.ubam)
        trim_ubam_r1 = trim_r1_tn5me.out.keep_bam
    } else {
        trim_ubam_r1 = fastq_to_ubam.out.ubam
    }

    trim_leading_attachment_site(trim_ubam_r1)

    if (params.trim_Tn5_r2) {
        trim_r2_tn5me(trim_leading_attachment_site.out.keep_bam)
        trim_ubam_r2 = trim_r2_tn5me.out.keep_bam
    } else if (params.trim_adapter_r2 != null) {
        trim_r2_adapter(trim_leading_attachment_site.out.keep_bam)
        trim_ubam_r2 = trim_r2_adapter.out.keep_bam
    } else {
        trim_ubam_r2 = trim_leading_attachment_site.out.keep_bam
    }

    if (split_fastq) {
        // remove the index from the meta so all the metas with same id
        // will group together
        bam_group = trim_ubam_r2
            .map { meta, bam -> tuple(meta.withoutIndex(), bam) }
            .groupTuple()
        concatenate_bams(bam_group)
        final_ubam = concatenate_bams.out.concatenated_bam
    } else {
        final_ubam = trim_ubam_r2
    }
    
    align_input = final_ubam.map { meta, bam ->
        tuple(meta, bam, meta.referenceDir, meta.referenceFastaName)
    }
    align(align_input)

    clip_input = align.out.mapped_bam.map { meta, bam, index ->
        tuple(meta, bam, index, meta.referenceFasta, meta.referenceFastaIndex)
    }
    clip_bam(clip_input)

    if (params.umi_from_read_name) {
        copy_umi_from_read_name(clip_bam.out.clipped_bam)
        mark_duplicates_input = copy_umi_from_read_name.out.bam_with_rx
    } else {
        mark_duplicates_input = clip_bam.out.clipped_bam
    }

    mark_duplicates(mark_duplicates_input)

    metrics_input = mark_duplicates.out.deduped_bam.map { meta, bam, index ->
        tuple(meta, bam, index, meta.referenceFasta, meta.referenceFastaIndex)
    }
    collect_alignment_summary_metrics(metrics_input)
    collect_multiple_metrics(metrics_input)
    
    multiqc_input = channel.empty().mix(
        fastqc.out.report_html,
        fastqc.out.report_zip,
        mark_duplicates.out.metric_txt,
        collect_alignment_summary_metrics.out.metric_txt,
        collect_multiple_metrics.out.metric_txt
    ).collect()
    multiqc(multiqc_input)

    find_sites_input = mark_duplicates.out.deduped_bam.map { meta, bam, index -> tuple(meta, bam) }
    find_sites(find_sites_input)

    aggregate_input = find_sites.out.sites_txt
        .multiMap { meta, sites_txt -> 
            meta: meta
            sites: sites_txt
        }
    all_meta = aggregate_input.meta.collect()
    all_sites = aggregate_input.sites.collect()

    collate_sites(all_meta, all_sites)

    annotate_fastas = aggregate_input.meta
        .map { meta -> meta.referenceFasta }
        .unique()
        .collect()
    annotate_gtfs = aggregate_input.meta
        .map { meta -> meta.geneSetGtf }
        .unique()
        .collect()
    match_reference = references
        .map { tuple(it.annotationDir, it.annotationFastaName) }
    annotate_sites(
        all_meta,
        all_sites,
        annotate_fastas,
        annotate_gtfs,
        match_reference
    )
}

process snakemake_pipeline {
    ext cpus: 8, memory: '16.GB'
    
    input:
        // the output prefix
        val prefix
        // metasheet - ref_dir column must have the reference file basename
        path metasheet
        // root dir where FASTQ files are located - this causes the entire
        // directory to be staged, recursively, so each experiment should have
        // its own directory structure for FASTQ files.
        path fastq_dir
        // list of reference folders, staged in a subdirectory
        path references, stageAs: 'references/*'
        // name of the annotation reference
        val annotation_reference
        // optional mapping of Benchling ID to attachment site
        path attachment_sites
        // optional global config - settings here override the default values
        // for snakemake global parameters, but are overriden by any parameters
        // set in the nextflow config or command line
        path global_config

    output:
        path 'sequencing_quality_report.html', emit: report
        path '*.txt', emit: summary_files
        path data_dir, emit: data_dir
        path 'logs', emit: log_dir
        path 'sequencing_quality_report_data', emit: report_data_dir
    
    script:
        sample_config = "${prefix}.yml"
        data_dir = "${prefix}_output"

        def fastq_pattern_opt = ""
        if (params.fastq_name_pattern != null) {
            // reformat from groovy-style to python-style
            def fastq_name_pattern = params.fastq_name_pattern.replace('${', '{')
            fastq_pattern_opt = "--fastq-name-pattern \"${fastq_name_pattern}\""
        }
        def attachment_sites_opt = ""
        if (attachment_sites != []) {
            attachment_sites_opt = "--attachment-sites \"${attachment_sites}\""
        }
        def default_attachment_site_opt = ""
        if (params.default_attachment_site != null) {
            default_attachment_site_opt = "--default-attachment-site \"${params.default_attachment_site}\""
        }
        
        def global_config_opt = ""
        if (global_config != []) {
            global_config_opt = "--config \"${global_config}\""
        }
        def annotation_ref_opt = ""
        if (annotation_reference != []) {
            fasta_path = "references/${annotation_reference}/${annotation_reference}.fasta.gz"
            annotation_ref_opt = "--annotation-fasta \"${fasta_path}\""
        }
        def trim_Tn5_opt = ""
        if (params.trim_Tn5 == true) {
            trim_Tn5_opt = "--trim-tn5"
        } else if (params.trim_Tn5 == false) {
            trim_Tn5_opt = "--no-trim-tn5"
        }
        def trim_Tn5_max_mismatches_opt = ""
        if (params.trim_Tn5_max_mismatches != null) {
            trim_Tn5_max_mismatches_opt = "--trim-Tn5-max-mismatches ${params.trim_Tn5_max_mismatches}"
        }
        def trim_att_max_mismatches_opt = ""
        if (params.trim_att_max_mismatches != null) {
            trim_att_max_mismatches_opt = "--trim-att-max-mismatches ${params.trim_att_max_mismatches}"
        }
        def read_structure_r1_opt = ""
        if (params.read_structure_r1 != null) {
            read_structure_r1_opt = "--read-structure-r1 ${params.read_structure_r1}"
        }
        def read_structure_r2_opt = ""
        if (params.read_structure_r2 != null) {
            read_structure_r2_opt = "--read-structure-r2 ${params.read_structure_r2}"
        }
        
        """
        # create the config file from the metasheet
        tomebio-tools cryptic-seq create-config-from-metasheet \
            --metasheet "${metasheet}" \
            --ref-dir references \
            --fastq-dir "${fastq_dir}" \
            ${fastq_pattern_opt} \
            ${attachment_sites_opt} \
            ${default_attachment_site_opt} \
            --output-file "${sample_config}" \
            --groups-file "${prefix}.groups.txt"
        
        # run the pipeline
        cryptic-seq run \
            ${global_config_opt} \
            ${annotation_ref_opt} \
            ${trim_Tn5_opt} \
            ${trim_Tn5_max_mismatches_opt} \
            ${trim_att_max_mismatches_opt} \
            ${read_structure_r1_opt} \
            ${read_structure_r2_opt} \
            --config-yml "${sample_config}" \
            --cores ${task.cpus} \
            --directory .
        
        # move group data directories into a subdirectory so we don't have to
        # know how they are named to output them
        mkdir ${data_dir}
        while read group; do
            mv "./\$group" "./${data_dir}/"
        done < "${prefix}.groups.txt"
        """
}

// Workflow that just calls the snakemake pipeline via a wrapper.
workflow snakemake_wrapper {
    // prepare the metasheet and references
    prepare_metasheet()

    prefix = params.prefix ?: prepare_metasheet.out.prefix
    metasheet = prepare_metasheet.out.metasheet
    fastq_dir = file(params.fastq_dir)
    // parse the references JSON file and get a list of reference paths
    // that need to be staged in
    references = prepare_metasheet.out.references
        .splitJson()
        .map { file(it['value']) }
        .toList()
    annotation_reference = prepare_metasheet.out.annotation_reference
    attachment_sites = []
    if (params.attachment_sites != null) {
        attachment_sites = file(params.attachment_sites)
    }
    global_config = []
    if (params.global_config != null) {
        global_config = file(params.global_config)
    }
    
    // run snakemake pipeline
    snakemake_pipeline(
        prefix,
        metasheet,
        fastq_dir,
        references,
        annotation_reference,
        attachment_sites,
        global_config
    )
}

// Given a CTB ID, extracts the associated metadata from Benchling and creates a metasheet.
process metasheet_from_benchling {
    label 'tools'
    label 'global'
    secret 'WAREHOUSE_USERNAME'
    secret 'WAREHOUSE_PASSWORD'

    input:
        path genomes_json
    
    output:
        path metasheet
    
    script:
        def prefix = params.prefix ?: "metasheet"
        metasheet = "${prefix}.txt"
        def genomes_json_opt = ""
        if (genomes_json!= []) {
            genomes_json_opt = "--genomes ${genomes_json}"
        }
        def default_genome_opt = ""
        if (params.default_genome != null) {
            "--default-genome ${params.default_genome}"
        }
        """
        tomebio-tools cryptic-seq create-metasheet-from-benchling \
            --ctb-id ${params.ctb_id} \
            ${genomes_json_opt} \
            ${default_genome_opt} \
            --output-file ${metasheet}
        """
}

// Given a metasheet, an optional root folder for references, and an optional 
// mapping between keys used in the metasheet ref_dir column and relative or
// absolute paths, returns a list of paths to the unique references and an 
// updated metasheet with the reference file name in the ref_dir column.
process resolve_references {
    label 'tools'
    label 'global'

    input:
        path metasheet
        path references_json
    
    output:
        path updated_metasheet, emit: metasheet
        path updated_references, emit: references
        path annotation_reference, emit: annotation
    
    script:
        def prefix = params.prefix ?: metasheet.simpleName
        updated_metasheet = "${prefix}.updated.txt"
        updated_references = "${prefix}.references.json"
        annotation_reference = "${prefix}.annotation.txt"
        def references_dir_opt = ""
        if (params.references_dir != null) {
            references_dir_opt = "--root ${params.references_dir}"
        }
        def references_json_opt = ""
        if (references_json!= []) {
            references_json_opt = "--mapping ${references_json}"
        }
        def default_reference_opt = ""
        if (params.default_reference != null) {
            default_reference_opt = "--default ${params.default_reference}"
        }
        def annotation_reference_opt = ""
        if (params.annotation_reference != null) {
            annotation_reference_opt = "--annotation ${params.annotation_reference}"
        }
        """
        tomebio-tools cryptic-seq resolve-references \
            --metasheet ${metasheet} \
            ${references_dir_opt} \
            ${references_json_opt} \
            ${default_reference_opt} \
            ${annotation_reference_opt} \
            --output-metasheet ${updated_metasheet} \
            --output-references ${updated_references} \
            --output-annotation ${annotation_reference}
        """
}

workflow prepare_metasheet {
    main:
        if (params.metasheet != null) {
            metasheet = file(params.metasheet)
        } else if (params.ctb_id != null) {
            // TODO: this code path should move inside of the Lambda function
            // that executes workflows
            genomes_json = []
            if (params.genomes_json != null) {
                genomes_json = file(params.genomes_json)
            }
            metasheet = metasheet_from_benchling(genomes_json)
        } else {
            exit(1, "Either 'metasheet' or 'ctb_id' must be specified.")
        }

        references_json = []
        if (params.references_json != null) {
            references_json = file(params.references_json)
        }
        // update metasheet with reference name, and extract unique references
        // to stage in the snakemake process
        resolve_references(metasheet, references_json)
    
    emit:
        prefix = metasheet.simpleName
        metasheet = resolve_references.out.metasheet
        references = resolve_references.out.references
        // read the file with the name of the annotation reference and return it as value
        annotation_reference = resolve_references.out.annotation
            .splitText() { it.trim() }
            .first()
}

workflow {
    // Right now we have two separate workflows - one that runs the
    // snakemake pipeline via a wrapper, and another that runs the
    // native Nextflow pipeline. For now we'll maintain both and default
    // to the native version.
    if (params.snakemake) {
        snakemake_wrapper()
    } else {
        cryptic_seq()
    }
}
