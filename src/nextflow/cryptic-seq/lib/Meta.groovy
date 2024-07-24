import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import groovy.text.SimpleTemplateEngine
import groovy.text.Template
import groovy.transform.CompileDynamic
import groovy.transform.Immutable
import java.nio.file.Files
import java.nio.file.Path
import nextflow.Nextflow
import nextflow.io.ValueObject
import nextflow.util.KryoHelper

@CompileDynamic
@ValueObject
@Immutable(copyWith = true, knownImmutables = ['referenceDir', 'fq1', 'fq2'])
class Meta extends Nextflow{
    static {
        KryoHelper.register(Meta)
    }

    // required metasheet column names
    static final String ID = 'sample_name'
    static final String REPLICATE = 'replicate'
    static final String GROUP = 'group'
    static final String REFERENCE = 'reference'
    static final String ATTACHMENT_SITES = 'attachment_sites'
    static final String ATTACHMENT_SITE_ID = 'benchling_pl_id'
    static final String FQ1 = 'fq1'
    static final String FQ2 = 'fq2'
    static final List<String> REQUIRED_FIELDS = [
        ID,
        REPLICATE,
        GROUP,
        REFERENCE
    ]
    static final String DEFAULT_FASTQ_PATTERN = '${sample_name}/${sample_name}_R${read}.fq.gz'

    /**
     * Parses a metasheet at the specified path and creates a `Meta` object for each sample.
     * @param path The path to the metasheet.
     * @param references The references used by samples in the metasheet.
     * @param attachmentSitesJson The path to the JSON file mapping benchling IDs to attachment sites.
     * @param defaultAttachmentSites The list of default attachment sites.
     * @param fastqDir The path to the base directory containing FASTQ files for this run.
     * @param fastqPattern The pattern for creating the FASTQ file name from the columns in the
     * metasheet. The pattern can have `${column}` placeholders that will be replaced with the
     * corresponding column value. The `${read}` placeholder will be replaced with the FASTQ file
     * pair ('1' or '2'). Defaults to '${sample_name}/${sample_name}_R${read}.fq.gz'.
     * @return A list of `Meta` objects.
     */
    static List<Meta> fromTxt(
        Path metasheet,
        References references,
        Path attachmentSitesJson,
        List<String> defaultAttachmentSites,
        Path fastqDir,
        String fastqPattern,
        Boolean verify
    ) {
        List<List<String>> rows = metasheet.readLines()
            .collect { line -> line.split('\t', -1).toList() }
            .findAll { row ->
                // remove comments and empty lines
                !(
                    row.empty ||
                    row[0].startsWith('#') ||
                    row.size() == 1 && row[0].trim().empty
                )
            }
        
        if (rows.size() < 2) {
            throw new RuntimeException("Metasheet is empty")
        }

        List<String> header = rows.head()*.toLowerCase()
        List<String> missing = REQUIRED_FIELDS.findAll { !header.contains(it) }
        if (!missing.empty) {
            throw new RuntimeException("Missing required columns: ${missing}")
        }

        JsonSlurper jsonSlurper = new JsonSlurper()

        Map<String, List<String>> attachmentSiteById = [:]
        if (attachmentSitesJson != null) {
            attachmentSiteById = jsonSlurper.parse(attachmentSitesJson)
        }

        SimpleTemplateEngine engine = new SimpleTemplateEngine()
        Template fastqTemplate = engine.createTemplate(fastqPattern ?: DEFAULT_FASTQ_PATTERN)

        rows.tail().collect { row ->
            Map<String, String> values = [header, row].transpose().collectEntries()
            
            String id = values[ID]
            Integer replicate = values[REPLICATE] as Integer
            String group = values[GROUP]

            // resolve reference dir
            String referenceName = values[REFERENCE]
            Path referenceDir = references.dir(referenceName)
            
            // resolve attachment sites
            List<String> attachmentSites = []
            if (values.containsKey(ATTACHMENT_SITES)) {
                attachmentSites = values[ATTACHMENT_SITES].split(';')
            } else if (values.containsKey(ATTACHMENT_SITE_ID)) {
                def missingKeys = []
                (attachmentSites, missingKeys) = values[ATTACHMENT_SITE_ID]
                    .split(/[;,]/)  // Split by either ',' or ';'
                    .collect { key -> key.trim() }  // Trim whitespace from each key
                    .inject([[], []]) { accu, key ->
                        if (key in attachmentSiteById) {
                            accu[0] << attachmentSiteById[key]
                        } else {
                            accu[1] << key
                        }
                        accu
                    }
                if (!missingKeys.empty) {
                    throw new RuntimeException("No attachment site found for: ${missingKeys}")
                }
            }
            if (!attachmentSites) {
                if (defaultAttachmentSites) {
                    attachmentSites = defaultAttachmentSites
                } else {
                    throw new RuntimeException("No attachment sites specified for sample ${id}")
                }
            }

            String fq1Name = values[FQ1] ?: fastqTemplate.make(values + [read: 1])
            Path fq1 = fastqDir.resolve(fq1Name)

            String fq2Name = values[FQ2] ?: fastqTemplate.make(values + [read: 2])
            Path fq2 = fastqDir.resolve(fq2Name)
            
            if (verify) {
                for (path in [referenceDir, fq1, fq2]) {
                    if (!Files.exists(path)) {
                        throw new RuntimeException("Path does not exist: ${path}")
                    }
                }
            }

            new Meta(
                id, replicate, group, attachmentSites, referenceDir, referenceName, fq1, fq2, null
            )
        }
    }

    /**
     * Creates a JSON file at `path` with one entry for each item in `metas`. Each item is an
     * object with the fields required by `collate_sites`: id, replicate, group.
     */
    static String toCollateJson(List<Meta> metas) {
        JsonOutput.toJson(
            metas.collect { meta ->
                [
                    name: meta.id,
                    replicate: meta.replicate,
                    group: meta.group
                ]
            }
        )
    }

    /**
     * Creates a JSON file at `path` with one entry for each item in `metas`. Each item is an
     * object with the fields required by `annotate_sites`: id, replicate, group, and ref_fasta.
     */
    static String toAnnotateJson(List<Meta> metas) {
        JsonOutput.toJson(
            metas.collect { meta ->
                [
                    name: meta.id,
                    replicate: meta.replicate,
                    group: meta.group,
                    ref_fasta: meta.referenceFastaName,
                    gene_set_gtf: meta.geneSetGtfName
                ]
            }
        )
    }

    /** The sample ID */
    String id
    /** Which replicate is the sample */
    Integer replicate
    /** The group to which the sample belongs */
    String group
    /** The list of attachment sites ("<name>:<left-seq>:<overhang>:<right-seq>") */
    List<String> attachmentSites
    /** The path to the directory that contains the reference FASTA and associated index files. */
    Path referenceDir
    /** The name of the reference. All reference files must have this prefix. */
    String referenceName
    /** Path to the read1 fastq file */
    Path fq1
    /** Path to the read2 fastq file */
    Path fq2
    /** Index that is set when the source FASTQs are split to make the ID unique. */
    Integer index

    /** 
     * Returns a copy of this Meta with `index` set to the specified value. 
     * Throws an exception if `index` is not `null`.
     */
    Meta withIndex(Integer index) {
        if (this.index != null) {
            throw new IllegalStateException("index has already been set")
        }
        return new Meta(
            id, replicate, group, attachmentSites, referenceDir, referenceName, fq1, fq2, index
        )
    }

    Meta withoutIndex() {
        if (index != null) {
            return new Meta(
                id, replicate, group, attachmentSites, referenceDir, referenceName, fq1, fq2, null
            )
        } else {
            return this
        }
    }

    String getUniqueId() {
        if (index != null) {
            return "${id}.${index}"
        } else {
            return id
        }
    }

    String getReferenceFastaName() {
        return "${referenceName}.fasta.gz"
    }

    Path getReferenceFasta() {
        return referenceDir.resolve(referenceFastaName)
    }

    String getGeneSetGtfName() {
        return "${referenceName}.gtf"
    }

    Path getGeneSetGtf() {
        return referenceDir.resolve(geneSetGtfName)
    }

    Path getReferenceFastaIndex() {
        return referenceDir.resolve("${referenceFastaName}.fai")
    }
}