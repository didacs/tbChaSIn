// Give a CTB ID, extracts the associated metadata from Benchling and creates a metasheet.
process metasheet_from_benchling {
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
    input:
        path metasheet
        path references_json
    
    output:
        path updated_metasheet, emit: metasheet
        path reference_paths, emit: references
    
    script:
        def prefix = params.prefix ?: metasheet.simpleName
        updated_metasheet = "${prefix}.updated.txt"
        reference_paths = "${prefix}.references.txt"
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
        """
        tomebio-tools cryptic-seq resolve-references \
            --metasheet ${metasheet} \
            ${references_dir_opt} \
            ${references_json_opt} \
            ${default_reference_opt} \
            --output-metasheet ${updated_metasheet} \
            --output-paths ${reference_paths}
        """
}

process cryptic_seq {
    publishDir "${params.output_dir ?: workflow.launchDir}/${prefix}"

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
            fastq_pattern_opt = "--fastq-name-pattern \"${params.fastq_name_pattern}\""
        }
        def attachment_sites_opt = ""
        if (attachment_sites != []) {
            attachment_sites_opt = "--attachment-sites \"${attachment_sites}\""
        }
        
        def global_config_opt = ""
        if (global_config != []) {
            global_config_opt = "--config \"${global_config}\""
        }
        def trim_tn5_opt = ""
        if (params.trim_tn5 == true) {
            trim_tn5_opt = "--trim-tn5"
        } else if (params.trim_tn5 == false) {
            trim_tn5_opt = "--no-trim-tn5"
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
            --output-file "${sample_config}" \
            --groups-file "${prefix}.groups.txt"
        
        # run the pipeline
        cryptic-seq run \
            ${global_config_opt} \
            ${trim_tn5_opt} \
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

workflow {
    if (params.metasheet != null) {
        metasheet = file(params.metasheet)
    } else if (params.ctb_id != null) {
        // TODO: this code path should move inside of the Lambda function
        // that executes workflows
        genomes_json = []
        if (params.genomes_json != null) {
            genomes_json = file(params.genomes_json)
        }
        metasheet = metasheet_from_benchling(params.ctb_id, genomes_json)
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
    // create a channel with the paths to all the references
    references = resolve_references.out.references
        .splitText { it.strip() }
        .map { path -> file(path) }
        .toList()
    
    prefix = params.prefix ?: metasheet.simpleName
    fastq_dir = file(params.fastq_dir)
    global_config = []
    if (params.global_config != null) {
        global_config = file(params.global_config)
    }
    attachment_sites = []
    if (params.attachment_sites != null) {
        attachment_sites = file(params.attachment_sites)
    }
    
    // run snakemake pipeline
    cryptic_seq(
        prefix,
        resolve_references.out.metasheet,
        fastq_dir,
        references,
        attachment_sites,
        global_config
    )
}