# tbChaSIn

Tools and pipelines for various off-target detection assays:

- CHANGE-Seq
- Cryptic-Seq
- Integration site mapping assay from Durant et al. 2022 

<!---toc start-->
- [tbChaSIn](#tbchasin)
  - [Local Setup](#local-setup)
    - [Running Tests](#running-tests)
  - [Run the Pipeline](#run-the-pipeline)
    - [Reference Preparation](#reference-preparation)
      - [CRISPR Reference Preparation](#crispr-reference-preparation)
      - [Integrase Reference Preparation](#integrase-reference-preparation)
    - [Execution](#execution)
      - [Config](#config)
        - [CHANGE-Seq](#change-seq)
        - [Cryptic-seq](#cryptic-seq)
        - [Durant et al.](#durant-et-al)
      - [Nextflow](#nextflow)

<!---toc end-->

## Local Setup

- [Install poetry][poetry-link]

It is important that poetry is *not* installed in the same environment as your package dependencies. We recommend using the "official" installer:

```console
curl -sSL https://install.python-poetry.org | python3 -
```

- [Install mamba][mamba-link]

We recommend using the [miniforge installer][miniforge-link]. Download the installer for your operating system and run it. For example:

```console
chmod +x Miniforge3-MacOSX-arm64.sh
./Miniforge3-MacOSX-arm64.sh
```

- Get a local copy of the `tbChaSIn` repo

```console
git clone git@github.com:tomebio/tbChaSIn.git
cd tbChaSIn
```

- Create the `tbChaSIn` conda environment

```console
mamba env create -f environment.yml
```

- Activate the `tbChaSIn` conda environment

```console
mamba activate tbChaSIn
```

- Install `pytomebio` (developer mode)

```console
poetry install
```

- Ensure `realpath` is available

On OSX:

```bash
brew install coreutils
```

On Ubuntu 16.04 or higher:

```bash
sudo apt-get install coreutils
```

To be able to use any of the tools that interact with the Benchling warehouse, you need to configure SSL/TLS as described [here][benchling-link].

### Running Tests

To run tests, execute:

```console
bash ci/precommit.sh [-f]
```

This will run:

1. Unit test for Python code and Snakemake plumbing (with `pytest`)
2. Linting of Python code (with `ruff`)
3. Code style checking of Python code (with `ruff`)
4. Type checking of Python code (with `mypy`)
5. Code style checking of Shell code (with `shellcheck`)
6. Code style checking of Snakemake code (with `snakefmt`)

If the optional `-f` flag is specified, then the Python and Snakemake files will be automatically formatted prior to applying the checks.

## Run the Pipeline

### Reference Preparation

#### CRISPR Reference Preparation

TBD

#### Integrase Reference Preparation

For Integrase samples, the GRCh38 (p14) genome is modified to append:

1. The full length attB sequence
2. The full length attP sequence
3. The full length attB-containing plasmid sequence


The following requires [`seqkit`](https://bioconda.github.io/recipes/seqkit/README.html) and the latest development version of fgbio (2.0.3 with git-hash `da9ecbcc` or higher):

```console
# Get the FASTA and assembly report.  The latter is needed to deterministicall sort and name
# the contigs downloaded in the FASTA
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt

# The fasta file needs to be unzipped (it can be re-compressed with bgzip if you like)
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz

# The fasta file needs to be indexed
samtools faidx GCF_000001405.40_GRCh38.p14_genomic.fna

# Create a sequence dictionary (.dict) to sort and name the contigs
fgbio CollectAlternateContigNames \
    -i GCF_000001405.40_GRCh38.p14_assembly_report.txt \
    -o GCF_000001405.40_GRCh38.p14_assembly_report.dict \
    -p UcscName \
    -a SequenceName AssignedMolecule GenBankAccession RefSeqAccession \
    -s AssembledMolecule UnlocalizedScaffold UnplacedScaffold AltScaffold

# Update the contig names and sort the contigs in the FASTA
fgbio UpdateFastaContigNames \
   -i GCF_000001405.40_GRCh38.p14_genomic.fna \
   -d GCF_000001405.40_GRCh38.p14_assembly_report.dict \
   -o GRCh38.p14.fasta \
   --sort-by-dict \
   --skip-missing

# Build the final fasta and index it
cat GRCh38.p14.fasta attB.fasta attP.fasta PL312.fasta | seqkit seq -w 60 - > GRCh38.p14.full.fasta
samtools faidx GRCh38.p14.full.fasta
samtools dict GRCh38.p14.full.fasta > GRCh38.p14.full.dict
bwa index GRCh38.p14.full.fasta
```

### Execution

Execute the following command to run the pipeline. `pipeline` is the name of one of the pipelines in the `src/snakemake` folder. Note that the `-d` argument is only required if you are executing the script from somewhere other than the root folder of the `tbChaSIn` project.

```console
bash src/scripts/run_snakemake.sh \
    [-d /path/to/snakemake/dir] \
    -p <pipeline> \
    -c /path/to/config.yml \
    [-g /path/to/global_config.yml] \
    -o /path/to/output \
    -t /path/to/large/temp/directory
```

See [Reference Preparation](#reference-preparation) for how to prepare the reference genome.

#### Config

There are two configuration files: run and global.

The run config is required and provides metadata for the samples to be processed.
The run config is organized at three levels: run, group, and sample level.
The run level contains configuration that applies to all groups and samples, for example tool-level parameters.
The group level contains configuration that applies to related samples, for example the guide for specific CRISPR samples, or tool-specific parameters recommended for integrase samples.
The sample level contains configuration that applies to each sample, for example the name, replicate number, and paths to the input FASTQs.

| Config Key         | Description                                                                                 | Level  | Required | Default |
| ------------------ | ------------------------------------------------------------------------------------------- | ------ | -------- | ------- |
| `name`             | The name of the group                                                                       | Group  | Yes      | NA      |
| `ref_fasta`        | The absolute path to the reference FASTA, with accompanying BWA index files and FASTA index | Group  | Yes      | NA      |
| `attachment_sites` | The list of attachment sites (`<name>:<left-seq>:<overhang>:<right-seq>`)                   | Group  | Yes      | NA      |
| `name`             | The name of the sample                                                                      | Sample | Yes      | NA      |
| `replicate`        | The replicate number (e.g. 1, 2, 3)                                                         | Sample | Yes      | NA      |
| `fq1`              | The absolute path to the FASTQ for read 1 (R1)                                              | Sample | Yes      | NA      |
| `fq2`              | The absolute path to the FASTQ for read 2 (R2)                                              | Sample | Yes      | NA      |

The global config is optional and overrides default values for parameters that apply to all samples.
The global parameters are specific to each workflow.

##### CHANGE-Seq

The following steps are performed in the CHANGE-Seq pipeline:

1. Convert FASTQ to BAM, trim and annotate the leading attachment sites.  Keep only read pairs where the left/right side 
   of the same attachment site are found.  
2. Search for the Tn5 mosaic end in the reads trimming it and all subsequent bases (due to short inserts).  All reads kept.
3. Align the reads.
4. Mark duplicates (no UMIs).
5. Collect various sequencing quality control metrics.
6. Collate putative integration sites per sample.
7. Collate putative integration sites across samples.

An example `config.yml` for CHANGE-Seq is shown below:

```yaml
settings:
  - name: bxb1
    ref_fasta: /path/to/GRCh38.p14.full/GRCh38.p14.full.fasta
    attachment_sites:
      - attB:CACCACGCGTGGCCGGCTTGTCGACGACGGCG:GT:CTCCGTCGTCAGGATCATCCGGGGATCCCGGG
      - attP:GCCGCTAGCGGTGGTTTGTCTGGTCAACCACCGCG:GT:CTCAGTGGTGTACGGTACAAACCCAGCTACCGGTC
    samples:
      - name: BxbI_alone_rep1_S10
        replicate: 1
        fq1: /path/to/cryptic-seq/fastqs/BxbI_alone_rep1_S10/BxbI_alone_rep1_S10_L001_R1_001.fastq.gz
        fq2: /path/to/cryptic-seq/fastqs/BxbI_alone_rep1_S10/BxbI_alone_rep1_S10_L001_R2_001.fastq.gz
      - name: BxbI_alone_rep2_S11
        replicate: 2
        fq1: /path/to/cryptic-seq/fastqs/BxbI_alone_rep2_S11/BxbI_alone_rep2_S11_L001_R1_001.fastq.gz
        fq2: /path/to/cryptic-seq/fastqs/BxbI_alone_rep2_S11/BxbI_alone_rep2_S11_L001_R2_001.fastq.gz
      - name: BxbI_alone_rep3_S12
        replicate: 3
        fq1: /path/to/cryptic-seq/fastqs/BxbI_alone_rep3_S12/BxbI_alone_rep3_S12_L001_R1_001.fastq.gz
        fq2: /path/to/cryptic-seq/fastqs/BxbI_alone_rep3_S12/BxbI_alone_rep3_S12_L001_R2_001.fastq.gz
      - name: BxbI_attB_rep1_S16
        replicate: 1
        fq1: /path/to/cryptic-seq/fastqs/BxbI_attB_rep1_S16/BxbI_attB_rep1_S16_L001_R1_001.fastq.gz
        fq2: /path/to/cryptic-seq/fastqs/BxbI_attB_rep1_S16/BxbI_attB_rep1_S16_L001_R2_001.fastq.gz
      - name: BxbI_attB_rep2_S17
        replicate: 2
        fq1: /path/to/cryptic-seq/fastqs/BxbI_attB_rep2_S17/BxbI_attB_rep2_S17_L001_R1_001.fastq.gz
        fq2: /path/to/cryptic-seq/fastqs/BxbI_attB_rep2_S17/BxbI_attB_rep2_S17_L001_R2_001.fastq.gz
      - name: BxbI_attB_rep3_S18
        replicate: 3
        fq1: /path/to/cryptic-seq/fastqs/BxbI_attB_rep3_S18/BxbI_attB_rep3_S18_L001_R1_001.fastq.gz
        fq2: /path/to/cryptic-seq/fastqs/BxbI_attB_rep3_S18/BxbI_attB_rep3_S18_L001_R2_001.fastq.gz
      - name: BxbI_attP_rep1_S13
        replicate: 1
        fq1: /path/to/cryptic-seq/fastqs/BxbI_attP_rep1_S13/BxbI_attP_rep1_S13_L001_R1_001.fastq.gz
        fq2: /path/to/cryptic-seq/fastqs/BxbI_attP_rep1_S13/BxbI_attP_rep1_S13_L001_R2_001.fastq.gz
      - name: BxbI_attP_rep2_S14
        replicate: 2
        fq1: /path/to/cryptic-seq/fastqs/BxbI_attP_rep2_S14/BxbI_attP_rep2_S14_L001_R1_001.fastq.gz
        fq2: /path/to/cryptic-seq/fastqs/BxbI_attP_rep2_S14/BxbI_attP_rep2_S14_L001_R2_001.fastq.gz
      - name: BxbI_attP_rep3_S15
        replicate: 3
        fq1: /path/to/cryptic-seq/fastqs/BxbI_attP_rep3_S15/BxbI_attP_rep3_S15_L001_R1_001.fastq.gz
        fq2: /path/to/cryptic-seq/fastqs/BxbI_attP_rep3_S15/BxbI_attP_rep3_S15_L001_R2_001.fastq.gz
      - name: HEK293c12_BxbI_attP_rep1_S25
        replicate: 1
        fq1: /path/to/cryptic-seq/fastqs/HEK293c12_BxbI_attP_rep1_S25/HEK293c12_BxbI_attP_rep1_S25_L001_R1_001.fastq.gz
        fq2: /path/to/cryptic-seq/fastqs/HEK293c12_BxbI_attP_rep1_S25/HEK293c12_BxbI_attP_rep1_S25_L001_R2_001.fastq.gz
      - name: HEK293c12_BxbI_attP_rep2_S26
        replicate: 2
        fq1: /path/to/cryptic-seq/fastqs/HEK293c12_BxbI_attP_rep2_S26/HEK293c12_BxbI_attP_rep2_S26_L001_R1_001.fastq.gz
        fq2: /path/to/cryptic-seq/fastqs/HEK293c12_BxbI_attP_rep2_S26/HEK293c12_BxbI_attP_rep2_S26_L001_R2_001.fastq.gz
      - name: HEK293c12_BxbI_attP_rep3_S27 
        replicate: 3
        fq1: /path/to/cryptic-seq/fastqs/HEK293c12_BxbI_attP_rep3_S27/HEK293c12_BxbI_attP_rep3_S27_L001_R1_001.fastq.gz
        fq2: /path/to/cryptic-seq/fastqs/HEK293c12_BxbI_attP_rep3_S27/HEK293c12_BxbI_attP_rep3_S27_L001_R2_001.fastq.gz
```

##### Cryptic-seq

The following steps are performed in the Cryptic-Seq pipeline:

1. Convert FASTQ to BAM, trim and annotate the leading UMIs
2. Trims the start of R1s for the Tn5 mosaic end, keeping only read pairs where the Tn5 mosaic end was found in R1.
3. Trims the start of R2s for the leading attachment site, keeping only read pairs where the attachment site was found in R2.
4. Search for the Tn5 mosaic end in the R2, trimming it and all subsequent bases (due to short inserts).  All reads kept.
5. Align the reads.
6. Clip reads in FR pairs that sequence past the far end of their mate.
7. Mark duplicates using the UMIs.
8. Collect various sequencing quality control metrics.
9. Collate putative integration sites per sample.
10. Collate putative integration sites across samples.

An example `config.yml` for CRYPTIC-Seq is shown below:

```yaml
settings:
  - name: Bxb1attP
    ref_fasta: /path/to/GRCh38.p14.full/GRCh38.p14.full.fasta
    attachment_sites:
      - attP:GTGGTTTGTCTGGTCAACCACCGCG:GT:CTCAGTGGTGTACGGTACAAACCCA
    samples:
      - name: JM-CS004a-N7-Bxb1attP-rep1-Sid3_S3
        replicate: 1
        fq1: /path/to/cryptic-seq/fastqs/JM-CS004a-N7-Bxb1attP-rep1-Sid3_S3_L001_R1_001.fastq.gz
        fq2: /path/to/cryptic-seq/fastqs/JM-CS004a-N7-Bxb1attP-rep1-Sid3_S3_L001_R2_001.fastq.gz
      - name: JM-CS004a-N7-Bxb1attP-rep2-Sid8_S8
        replicate: 2
        fq1: /path/to/cryptic-seq/fastqs/JM-CS004a-N7-Bxb1attP-rep2-Sid8_S8_L001_R1_001.fastq.gz
        fq2: /path/to/cryptic-seq/fastqs/JM-CS004a-N7-Bxb1attP-rep2-Sid8_S8_L001_R2_001.fastq.gz
  - name: Bxb1attB
    ref_fasta: /path/to/GRCh38.p14.full/GRCh38.p14.full.fasta
    attachment_sites:
      - attB:GGCCGGCTTGTCGACGACGGCG:GT:CTCCGTCGTCAGGATCATCCGG
    samples:
      - name: JM-CS004a-N7-Bxb1attB-rep1-Sid4_S4
        replicate: 1
        fq1: /path/to/cryptic-seq/fastqs/JM-CS004a-N7-Bxb1attB-rep1-Sid4_S4_L001_R1_001.fastq.gz
        fq2: /path/to/cryptic-seq/fastqs/JM-CS004a-N7-Bxb1attB-rep1-Sid4_S4_L001_R2_001.fastq.gz
      - name: JM-CS004a-N7-Bxb1attB-rep2-Sid9_S9
        replicate: 2
        fq1: /path/to/cryptic-seq/fastqs/JM-CS004a-N7-Bxb1attB-rep2-Sid9_S9_L001_R1_001.fastq.gz
        fq2: /path/to/cryptic-seq/fastqs/JM-CS004a-N7-Bxb1attB-rep2-Sid9_S9_L001_R2_001.fastq.gz
  - name: attBalone
    ref_fasta: /path/to/GRCh38.p14.full/GRCh38.p14.full.fasta
    attachment_sites:
      - attB:GGCCGGCTTGTCGACGACGGCG:GT:CTCCGTCGTCAGGATCATCCGG
    samples:
      - name: JM-CS004a-N7-attBalone-rep1-Sid2_S2
        replicate: 1
        fq1: /path/to/cryptic-seq/fastqs/JM-CS004a-N7-attBalone-rep1-Sid2_S2_L001_R1_001.fastq.gz
        fq2: /path/to/cryptic-seq/fastqs/JM-CS004a-N7-attBalone-rep1-Sid2_S2_L001_R2_001.fastq.gz
      - name: JM-CS004a-N7-attBalone-rep2-Sid7_S7
        replicate: 2
        fq1: /path/to/cryptic-seq/fastqs/JM-CS004a-N7-attBalone-rep2-Sid7_S7_L001_R1_001.fastq.gz
        fq2: /path/to/cryptic-seq/fastqs/JM-CS004a-N7-attBalone-rep2-Sid7_S7_L001_R2_001.fastq.gz
  - name: attPalone
    ref_fasta: /path/to/GRCh38.p14.full/GRCh38.p14.full.fasta
    attachment_sites:
      - attP:GTGGTTTGTCTGGTCAACCACCGCG:GT:CTCAGTGGTGTACGGTACAAACCCA
    samples:
      - name: JM-CS004a-N7-attPalone-rep1-Sid1_S1
        replicate: 1
        fq1: /path/to/cryptic-seq/fastqs/JM-CS004a-N7-attPalone-rep1-Sid1_S1_L001_R1_001.fastq.gz
        fq2: /path/to/cryptic-seq/fastqs/JM-CS004a-N7-attPalone-rep1-Sid1_S1_L001_R2_001.fastq.gz
      - name: JM-CS004a-N7-attPalone-rep2-Sid6_S6
        replicate: 2
        fq1: /path/to/cryptic-seq/fastqs/JM-CS004a-N7-attPalone-rep2-Sid6_S6_L001_R1_001.fastq.gz
        fq2: /path/to/cryptic-seq/fastqs/JM-CS004a-N7-attPalone-rep2-Sid6_S6_L001_R2_001.fastq.gz
```

The global parameters supported for Cryptic-seq are listed in the following table. A configuration file with the default values is [provided](src/snakemake/cryptic_seq.config.yaml).

| Config Key              | Description                                                                                                 | Default |
|-------------------------|-------------------------------------------------------------------------------------------------------------|----|
| read_structure_r1       | R1 read structure for _fgbio FastqToBam_, made up of '\<number>\<operator>' pairs                           | 11M+T |
| read_structure_r2       | R2 read structure for _fgbio FastqToBam_, made up of '\<number>\<operator>' pairs                           | +T |
| trim_Tn5                | Whether to trim the Tn5 mosiac end sequence from the start of R1                                            | `True` |
| trim_Tn5_max_mismatches | Maximum number of mismatches to allow Tn5 trimming                                                          | 1  |
| trim_att_max_mismatches | Maximum number of mismatches to allow leading attachment site trimming                                      | 4  |

##### Durant et al.

The following steps are performed in the Durant et al. pipeline:

1. uses `fastp` to trim adapters (r1: custom, r2: nextera)
2. uses `tomebio-tools durant trim-leading-r2` to keep only reads where R2 starts with the sample-specific stagger and inner donor primer
3. aligns to the genome, which includes attB and attD
4. keeps only reads that:
    - R1 and R2 must each have at least one alignment to the genome (non-donor)
    - If R1 has an alignment to the donor, it cannot have too many mapped bases to the donor,
      otherwise it is assumed to be a linear plasmid template (default: 55bp)
    - R2 must have an alignment to the donor on the forward strand
5. re-aligns these reads to the genome without attB and attD
6. does not de-duplicate
7. uses the `tomebio-tools durant find-sites` tool to find integration sites.  Keeps templates that:
  - both R1 and R2 map to the genome
  - the template length is 1kb or smaller (implies R1 and R2 map to the same contig!)
  - R1 has enough mapped bases (default: 25bp) 
For a given template (read pair), the following must be true to be considered an integration:

An example `config.yml` for Durant et al. is shown below:

```yaml
settings:
- name: Durant
  ref_fasta: /path/to/GRCh38.p14.full/GRCh38.p14.full.fasta
  genome_fasta: /path/to/GRCh38.p14/GRCh38.p14.fasta
  min_aln_score: 20
  inter_site_slop: 10
  samples:
  - name: SRR21306552
    replicate: 1
    bio_rep: 2
    tech_rep: 2
    stagger: ''
    donor_inner_primer: CAGCGAGTCAGTGAGCGAGG
    umi_length: 0
    r1_adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    r2_adapter: CTGTCTCTTATACACATCTGACGCTGCCGACGA
    fq1: /path/to/SRR21306552_1.fastq.gz
    fq2: /path/to/SRR21306552_2.fastq.gz
  - name: SRR21306553
    replicate: 2
    bio_rep: 2
    tech_rep: 1
    stagger: T
    donor_inner_primer: CAGCGAGTCAGTGAGCGAGG
    umi_length: 0
    r1_adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    r2_adapter: CTGTCTCTTATACACATCTGACGCTGCCGACGA
    fq1: /path/to/SRR21306553_1.fastq.gz
    fq2: /path/to/SRR21306553_2.fastq.gz
  - name: SRR21306554
    replicate: 3
    bio_rep: 1
    tech_rep: 1
    stagger: ATCGAT
    donor_inner_primer: TCGATCGAGGTTGCATTCGG
    umi_length: 0
    r1_adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    r2_adapter: CTGTCTCTTATACACATCTGACGCTGCCGACGA
    fq1: /path/to/SRR21306554_1.fastq.gz
    fq2: /path/to/SRR21306554_2.fastq.gz
  - name: SRR21306560
    replicate: 4
    bio_rep: 5
    tech_rep: 1
    stagger: T
    donor_inner_primer: CAGCGAGTCAGTGAGCGAGG
    umi_length: 0
    r1_adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    r2_adapter: CTGTCTCTTATACACATCTGACGCTGCCGACGA
    fq1: /path/to/SRR21306560_1.fastq.gz
    fq2: /path/to/SRR21306560_2.fastq.gz
  - name: SRR21306561
    replicate: 5
    bio_rep: 2
    tech_rep: 3
    stagger: ''
    donor_inner_primer: CAGCGAGTCAGTGAGCGAGG
    umi_length: 0
    r1_adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    r2_adapter: CTGTCTCTTATACACATCTGACGCTGCCGACGA
    fq1: /path/to/SRR21306561_1.fastq.gz
    fq2: /path/to/SRR21306561_2.fastq.gz
  - name: SRR21306626
    replicate: 6
    bio_rep: 4
    tech_rep: 1
    stagger: CGAT
    donor_inner_primer: TCGATCGAGGTTGCATTCGG
    umi_length: 12
    r1_adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    r2_adapter: CTGTCTCTTATACACATCTGACGCTGCCGACGA
    fq1: /path/to/SRR21306626_1.fastq.gz
    fq2: /path/to/SRR21306626_2.fastq.gz
  - name: SRR21306627
    replicate: 7
    bio_rep: 3
    tech_rep: 1
    stagger: T
    donor_inner_primer: TCGATCGAGGTTGCATTCGG
    umi_length: 12
    r1_adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
    r2_adapter: CTGTCTCTTATACACATCTGACGCTGCCGACGA
    fq1: /path/to/SRR21306627_1.fastq.gz
    fq2: /path/to/SRR21306627_2.fastq.gz
```

*** Important ***: `ref_fasta` is the genome _with_ attD, while `genome_fasta` _does not_ contain attD.

#### Nextflow

These pipelines are in the process of being ported to Nextflow. The first step is to wrap the Snakemake pipeline using [snk][snk-link] inside a Docker container and call it from a Nextflow process. There is an initial step that creates the configuration file for the Snakemake workflow from a metasheet.

To build the Docker container:

```
docker build -t tomebio/cryptic-seq:1.0 .
```

On an ARM-based Mac, some additional options are required:

```
docker build \
  --platform linux/amd64 \
  --load \
  -t tomebio/cryptic-seq:1.0 .
```

To run the workflow:

```
nextflow run \
  src/nextflow/cryptic-seq \
  --metasheet <metasheet.[xlsx|txt]> \
  --genomes_json <genomes.json> \
  --references_json <references.json> \
  --fastq_dir <fastq_dir> \
  [--output_dir <output_dir> ] \
  [--prefix <output_prefix>] \
  [-with-report report.html] \
  # this option only required on ARM Mac
  [-profile rosetta] \
  # this option only required on linux systems where it is required to run docker as root
  [-profile linux]  
```

The `genomes.json` file contains mappings between species name and genome build.

The `references.json` file contains mappings between genome build and path to the folder that contains the reference (FASTA file and BWA index).

The `fastq_dir` is the root directory where FASTQ files live. The FASTQ file names may be specified in the metasheet `fq1` and `fq2` columns as relative paths under the `fastq_dir`, or the `--fastq_name_prefix` option may be specified with a glob expression that can contain placeholders for any of the columns in the metasheet, as well as the special `read` placeholder which has a value of `1` for read 1 files and `2` for read 2 files, e.g. `**/{sample_name}*/*_R{read}_*.gz`.

The `output_dir` is the directory where pipeline outputs will be published. It defaults to the directory where the workflow is launched.

The `output_prefix` is the name of the subdirectory within `output_dir` where pipeline outputs will be published, and is also used to name the run-level outputs. It defaults to the name of the `metasheet` (without extension).

[poetry-link]: https://python-poetry.org/docs/#installation
[mamba-link]: https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html
[miniforge-link]: https://github.com/conda-forge/miniforge
[snk-link]: https://github.com/Wytamma/snk
[benchling-link]: https://docs.benchling.com/docs/getting-started