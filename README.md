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

<!---toc end-->

## Local Setup

- [Install conda][conda-link]

- Install mamba

```console
conda install -c conda-forge mamba
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
python setup.py develop
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

### Running Tests

To run tests, execute:

```console
bash ci/precommit.sh
```

This will run:

1. Unit test for Python code and Snakemake plumbing (with `pytest`)
2. Linting of Python code (with `flake8`)
3. Code style checking of Python code (with `black`)
4. Type checking of Python code (with `mypy`)
5. Code style checking of Shell code (with `shellcheck`)

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

# Create a sequence dictionary (.dict) to sort and name the contigs
java -Xmx2g -jar ~/work/git/fgbio/target/scala-2.13/fgbio-2.0.3-da9ecbcc-SNAPSHOT.jar CollectAlternateContigNames \
    -i GCF_000001405.40_GRCh38.p14_assembly_report.txt \
    -o GCF_000001405.40_GRCh38.p14_assembly_report.dict \
    -p UcscName \
    -a SequenceName AssignedMolecule GenBankAccession RefSeqAccession \
    -s AssembledMolecule UnlocalizedScaffold UnplacedScaffold AltScaffold \
   --sort-by-sequencing-role

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

Execute the following command to run the pipeline

```console
bash src/scripts/run_snakemake.sh \
    -t /path/to/large/temp/directory \
    -s src/snakemake/<pipeline>smk \
    -c /path/to/config.yml \
    -o /path/to/output
```

See [Reference Preparation](#reference-preparation) for how to prepare the reference genome.

#### Config

The config is organized at three levels: run, group, and sample level.
The run level contains configuration that applies to all groups and samples, for example tool-level parameters
The group level contains configuration that applies to related samples, for example the guide for specific CRISPR samples, or tool-specific parameters recommended for integrase samples.
The sample level contains configuration that applies to each sample, for example the name, replicate number, and paths to the input FASTQs.

| Config Key         | Description                                                                                 | Level  | Required | Default |
|--------------------|---------------------------------------------------------------------------------------------|--------|----------|---------|
| `name`             | The name of the group                                                                       | Group  | Yes      | NA      |
| `ref_fasta`        | The absolute path to the reference FASTA, with accompanying BWA index files and FASTA index | Group  | Yes      | NA      |
| `attachment_sites` | The list of attachment sites (`<name>:<left-seq>:<overhang>:<right-seq>`)                   | Group  | Yes      | NA      |
 | `name`             | The name of the sample                                                                      | Sample | Yes      | NA      |
 | `replicate`        | The replicate number (e.g. 1, 2, 3)                                                         | Sample | Yes      | NA      |
 | `fq1`              | The absolute path to the FASTQ for read 1 (R1)                                              | Sample | Yes      | NA      |
 | `fq2`              | The absolute path to the FASTQ for read 2 (R2)                                              | Sample | Yes      | NA      |

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


[conda-link]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/
