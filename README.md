# tbChaSIn

Tools and pipelines for various off-target detection assays:

- CHANGE-Seq
- Cryptic-Seq
- Integration site mapping assay from Durant et al. 2022

<!---toc start-->
- [tbChaSIn](#tbchasin)
  - [Local Setup](#local-setup)
    - [Running Code Checks](#running-code-checks)
  - [Design](#design)
  - [Run the Pipeline](#run-the-pipeline)
    - [Reference Preparation](#reference-preparation)
      - [CRISPR Reference Preparation](#crispr-reference-preparation)
      - [Integrase Reference Preparation](#integrase-reference-preparation)
    - [Execution (Nextflow)](#execution-nextflow)
      - [Building Docker images](#building-docker-images)
      - [Setting secrets](#setting-secrets)
      - [Running the workflow](#running-the-workflow)
    - [Execution (Snakemake)](#execution-snakemake)
      - [Config](#config)
        - [CHANGE-Seq](#change-seq)
        - [Cryptic-seq](#cryptic-seq)
        - [Durant et al.](#durant-et-al)
    - [Development](#development)
      - [Python toolkit](#python-toolkit)
      - [Converting Snakemake to Nextflow](#converting-snakemake-to-nextflow)
        - [Adding a process](#adding-a-process)
      - [Integration tests](#integration-tests)

<!---toc end-->

## Local Setup

- [Install `poetry`][poetry-link]

It is important that `poetry` is *not* installed in the same environment as your package dependencies. We recommend using the "official" installer:

```console
curl -sSL https://install.python-poetry.org | python3 -
```

If you already have poetry installed, make sure it is version `1.7.0` or greater:

```console
poetry --version
```

- [Install `mamba`][mamba-link]

We recommend using the [miniforge installer][miniforge-link]. Download the installer for your operating system and run it. For example:

```console
chmod +x Miniforge3-MacOSX-arm64.sh
./Miniforge3-MacOSX-arm64.sh
```

If you already have `mamba` installed, make sure it is version `1.3.0` or greater:

```console
mamba --version
```

- Get a local copy of the `tbChaSIn` repo

```console
git clone git@github.com:tomebio/tbChaSIn.git
cd tbChaSIn
```

- Create and activate the environment

```console
bash scripts/mamba.sh -A [-m micromamba]
```

You can also activate the environment manually:

```console
mamba activate tbChaSIn
```

- Install `pytomebio` (developer mode)

```console
cd src/python && poetry install
```

- If you plan to run the Snakemake workflows, ensure `realpath` is available

On OSX:

```bash
brew install coreutils
```

On Ubuntu 16.04 or higher:

```bash
sudo apt-get install coreutils
```

- To be able to use any of the tools that interact with the Benchling warehouse, you need to configure SSL/TLS as described [here][benchling-link].

### Running Code Checks

To run code checks, execute:

```console
bash scripts/precommit.sh [-f] [-l]
```

This will run:

1. Unit test for Python code and Snakemake plumbing (with `pytest`)
2. Linting of Python code (with `ruff`)
3. Code style checking of Python code (with `ruff`)
4. Type checking of Python code (with `mypy`)
5. Code style checking of Shell code (with `shellcheck`)
6. Code style checking of Snakemake code (with `snakefmt`)

If the optional `-f` flag is specified, then the Python and Snakemake files will be automatically formatted prior to applying the checks. Similarly, use the `-l` flag to fix any lints that can be fixed automatically.

## Design

The Nextflow workflows are designed to run either locally or using AWS Batch. To facilitate this, we must consider that input files (e.g., FASTQs and references) may be stored remotely (e.g., in AWS S3) and "staged" (i.e., Nextflow manages the automatic downloading of remote inputs to the local system) at runtime.

The minimal information to run the pipeline is:

* The CTB ID of the experiment, which is used to look up the experiment metadata in Benchling
* The root location of the references used in the experiment
* The root location of the FASTQ files for the experiment

The workflow runs the following steps:

1. Retrieve the metadata from Benchling for the given CTB ID, specified by the `--ctb_id` option, and create a "metasheet" (aka a sample sheet).
   1. If you already have a metasheet in the correct format, you can specify it using `--metasheet` rather than using the `--ctb_id` option.
2. For each sample in the metasheet:
   1. If the genome build is not known, then use the mapping between species and genome build (specified by the `--genomes_json` option) to determine it, or fall back to the default genome (specified by the `--default_genome` option).
   2. Map the genome build to a reference using the mapping specified by the `--references_json` option, otherwise assume the reference has the same name as the genome build, or fall back to the default reference (specified by the `--default_reference` option) if the genome build is not known>
   3. If the reference is specified as a relative path, then resolve it to an absolute path using the root folder specified by the `--references_dir` option. Note that `--references_dir` may specify a URI, such as an S3 bucket.
   4. Determine the names of the FASTQ files using the pattern string specified by `--fastq_name_pattern`. This pattern string can contain variables of the form `${var}`, which are replaced with the value from the corresponding column in the metasheet. There is also a special `${read}` variable with the read number for the FASTQ (1 or 2).
   5. If `--fastq_name_pattern` specifies a relative path, then resolve it to an absolute path using the root folder specified by `--fastq_dir`, which may also specify a URI, such as an S3 bucket.
3. Load each sample in the metasheet into a [native Groovy object](src/nextflow/cryptic-seq/lib/Meta.groovy).
4. Run all of the sample-specific steps in the pipeline in parallel over all the samples.
5. Aggregate the per-sample results and run the aggregate steps in the pipeline.

## Run the Pipeline

### Reference Preparation

#### CRISPR Reference Preparation

TBD

#### Integrase Reference Preparation

For Integrase samples, a reference genome (e.g. GRCh38.p14) is modified to append:

1. The full length attB sequence
2. The full length attP sequence
3. The full length attB-containing plasmid sequence

The [build_reference.sh](scripts/build_reference.sh) script will build the reference for you:

```console
bash scripts/build_reference.sh \
  [-i <reference folder URI> | -g <genome_name>] \
  [-n <reference name>] \
  [-D <root dir> | -G <genomes_dir> -R <references_dir>] \
  [-A] \
  attB.fasta attP.fasta plasmid.fasta
```

where:

* `-i` specifies the folder from which the genome files can be downloaded. Right now, this is assumed to be an NCBI URI, and defaults to https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14. The folder must contain a `<name>_genomic.fna.gz` file with the genome sequence, and a `<name>_assembly_report.txt` file with the contig names, where `<name>` is the folder name (e.g. `GCF_000001405.40_GRCh38.p14` in the default URI).
* `-g` specifies the genome name. Use this instead of `-i` if the genome files already exist in the genomes folder. Alternatively, the script will try to infer the NCBI URI from the genome name.
* `-n` specifies the name of the reference to create. This defaults to the genome name.
* `-D` specifies a root directory where both genomes and references will be stored, in the `<root>/genomes/` and `<root>/references` directories, respectively. This defaults to the system temp dir and is not required if you specify both `-G` and `-R`.
* `-G` specifies the directory where genomes are downloaded.
* `-R` specifies the directory where references are created.
* `-A` specifies that alt contigs should not be included in the final reference.

All the files for a reference are created in a directory with the reference name, and with all files prefixed with the reference name. For example, with reference name `GRCh38.p14_PL2312`, the reference directory will contain:

* GRCh38.p14_PL2312.dict
* GRCh38.p14_PL2312.gtf
* GRCh38.p14_PL2312.fasta.gz
* GRCh38.p14_PL2312.fasta.gz.amb
* GRCh38.p14_PL2312.fasta.gz.ann
* GRCh38.p14_PL2312.fasta.gz.bwt
* GRCh38.p14_PL2312.fasta.gz.fai
* GRCh38.p14_PL2312.fasta.gz.gzi
* GRCh38.p14_PL2312.fasta.gz.pac

### Execution (Nextflow)

The Cryptic-seq pipeline has been ported to Nextflow, and this is now the preferred method of execution. For the other pipelines, use [Snakemake](#execution-snakemake). For regression testing, the Nextflow workflow can run a wrapped version of the original Snakemake workflow using the `-profile snakemake` option.

#### Building Docker images

Each process uses a single Docker image, but multiple processes can use the same image. Each environment configuration file in [mamba/](mamba/) corresponds to a Docker image. To build all the Docker images, run:

```console
bash scripts/docker.sh [-m] [-v VERSION] [-x] [-p]
```

The `-m` option builds images that are compatible with an ARM Mac. The `-v` option specifies the version with which to tag the images, and must match the `docker_image_version` parameter in [nextflow.config](nextflow.config).

The `-x` option tags images with "prod" (in production) image names.
The `-p` option pushes tagged images to ECR.

To build a specific image, run the following command, where `<target>` is the name of the environment configuration file without the `.yml` extension.

```console
bash scripts/docker.sh [-m] [-v VERSION] [-x] <target>
```

To be able to run the Snakemake version of the workflow, you'll need to build a single "monolithic" Docker image instead. Run the following from the root directory of the project:

```
docker build \
  -f docker/Dockerfile.snakemake \
  -t tomebio/cryptic-seq:1.0 .
```

On an ARM-based Mac, some additional options are required:

```
docker build \
  -f docker/Dockerfile.snakemake \
  -t tomebio/cryptic-seq:1.0 \
  --platform linux/amd64 \
  --load .
```

#### Setting secrets

If you plan to use CTB IDs to fetch metadata from benchling, then you must configure secrets for the Benchling Warehouse credentials:

```console
nextflow secrets set WAREHOUSE_USERNAME 'username'
nextflow secrets set WAREHOUSE_PASSWORD 'password'
nextflow secrets set WAREHOUSE_HOST 'postgres-warehouse.tome.benchling.com'
nextflow secrets set WAREHOUSE_PORT '5432'
nextflow secrets set WAREHOUSE_DBNAME 'warehouse'
nextflow secrets set WAREHOUSE_SSLMODE 'verify-ca'
```

Note: make sure to enclose the secret in single-quotes.

#### Running the workflow

```
nextflow run \
  src/nextflow/cryptic-seq \
  [--metasheet <metasheet.[xlsx|txt]> | --ctb_id <ctb_id>] \
  --references_json <references.json> \
  --reference_dir <ref_dir> \
  --annotation_reference <name_or_path> \
  --fastq_dir <fastq_dir> \
  [--fastq_name_pattern <pattern>] \
  [--output_dir <output_dir> ] \
  [--prefix <output_prefix>] \
  [-with-report report.html] \
  # this option configures the pipeline for running locally
  -profile local \
  # this option configures the pipeline to run with "prod" or "dev" docker images
  # use "prod" if you are unsure which images to use.
  -profile {dev|prod}
  # this option only required on ARM Mac
  [-profile rosetta] \
  # this option only required on linux systems where it is required to run docker as root
  [-profile linux]
  # this option causes the snakemake workflow to be run
  [-profile snakemake] \
  # this option provides the global config when running the snakemake workflow
  [--global_config <config_yml> ]
```

The most common options are:

* `--metasheet` specifies the sample sheet. Alternatively, `--ctb_id` can be used to specify a CTB ID in benchling from which to create the metasheet. This may require specifying `--genomes_json`, which is a mapping between species and reference name.
* `--references_json` specifies a file that contains mappings between genome build and path to the folder that contains the reference (FASTA file and BWA index).
* `--reference_dir` specifies the root folder that contains references. This is only necessary if `--references_json` contains relative paths.
* `--annotation_reference` the name of, or path to, the reference to use in the `annotate_sites` step when searching for exact site sequence hits in the genome.
* `--fastq_dir` is the root directory where FASTQ files live. The FASTQ file names may be specified in the metasheet `fq1` and `fq2` columns as relative paths under the `fastq_dir`, or the `--fastq_name_prefix` option may be specified with a glob expression that can contain placeholders for any of the columns in the metasheet, as well as the special `read` placeholder which has a value of `1` for read 1 files and `2` for read 2 files, e.g. `**/{sample_name}*/*_R{read}_*.gz`.
* `--output_dir` is the directory where pipeline outputs will be published. It defaults to the directory where the workflow is launched.
* `--output_prefix` is the name of the subdirectory within `output_dir` where pipeline outputs will be published, and is also used to name the run-level outputs. It defaults to the name of the `metasheet` (without extension).

To see the full list of options that can be set, run:

```console
nextflow run src/nextflow/cryptic-seq -h
```

Note that the options that are specific to the nextflow command only start with one dash (e.g., `-with-report`) while options that override workflow parameters start with two dashes (e.g. `--metasheet`).

Instead of (or in addition to) setting options on the command line, you can also put them in a configuration file. Options specified on the command line override those in config files.

`example.config`
```
metasheet = "my_samples.txt"
references_json = "my_references.json"
trim_Tn5 = true
```

```console
nextflow run \
  src/nextflow/cryptic-seq \
  -c example.config \
  --fastq_dir my_fastq_dir
  --output_dir my_output_dir
  -profile local,prod
```

### Execution (Snakemake)

Execute the following command to run the pipeline. `pipeline` is the name of one of the pipelines in the `src/snakemake` folder. Note that the `-d` argument is only required if you are executing the script from somewhere other than the root folder of the `tbChaSIn` project.

```console
bash scripts/run_snakemake.sh \
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
| umi_from_read_name | Set to `True` if the UMI is found in the read name instead of the sequence        | `False`  |

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

---

## Development

#### Python toolkit

All custom scripts used by these workflows are written in python and available as a command line toolkit called `pytomebio`. To add a new tool:

* Create a `<tool>.py` file in the appropriate package within `pytomebio.tools`
* Add the tool to the `TOOLS` dict in `pytomebio.__main__.py`
* Add any new dependencies to `pyproject.toml`

#### Converting Snakemake to Nextflow

We follow a two-step process:

1. Wrap the Snakemake pipeline using [snk][snk-link] inside a Docker container and call it from a Nextflow process.
2. Convert each Snakemake rule to a corresponding Nextflow process and write a "native" Nextflow workflow.

The workflow then has the option of running either the "native" version or the Snakemake version. This is useful to be able to compare the results for regression testing.

##### Adding a process

To add a new step to the workflow, create a new `process` in the appropriate `main.nf` file in one of the subfolders of `src/nextflow`. The process should have the following directives:

* `label` corresponding to the Docker image that should be used
* `tag "${meta.id}"` for processes that run on individual samples
* `label 'global'` for processes that run on the aggregated results
* `ext resource: value` to specify resource requirements that are above the base level of 1 CPU, 2 GB memory, and 100 GB of disk space

Each process needs to have an associated Docker image. You can use one of the existing images by giving your process the associated `label` directive. If you need to create a new image:

* Create a new Mamba configuration file in [mamba/](mamba/)
* Create a new target in [Dockerfile](Dockerfile)
* Add the necessary configuration for the new image in `src/nextflow/<workflow>/nextflow.config`

#### Integration tests

First prepare test genomes with:

```console
bash scripts/build_reference.sh -D ./data -g GCF_000001405.40_GRCh38.p14 -c chr1 -x
bash scripts/build_reference.sh -D ./data -g GCF_000001635.27_GRCm39 -c chr1 -x
```

From `src/nextflow/cryptic-seq` run:

```console
REF_DIR=<ref_dir> PROFILES=local,test[,rosetta|linux] pytest --git-aware tests/integration/
```

where `ref_dir` is the root directory for references.

[poetry-link]: https://python-poetry.org/docs/#installation
[mamba-link]: https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html
[miniforge-link]: https://github.com/conda-forge/miniforge
[snk-link]: https://github.com/Wytamma/snk
[benchling-link]: https://docs.benchling.com/docs/getting-started

---

### Cutting a new release and deploying to Seqera Platform

To cut a new release of the pipeline to use on Seqera Platform, you must: update the pipeline metadata with the appropriate [semantic version](https://semver.org), update the pipeline's ECR images, and update the corresponding launch template(s) on Seqera Platform.

To update the repository with a new semantic version:
1. Bump the [semantic version](https://semver.org) of the pipeline in [src/nextflow/cryptic-seq/pipeline.env](src/nextflow/cryptic-seq/pipeline.env)
The semver should _not_ have a leading "v".
1. Merge the changes to `main`.
1. Tag `main` with the same version and push the tag. I.e   `git tag -a "X.Y.Z"`, then `git push origin tag "X.Y.Z"`.
1. Build and push prod or dev images using `scripts/docker.sh` with the `-p` option (with the `-x` version if pushing to production).

To update launch templates in Seqera Platform workspace(s):
1. Open the options drop-down to the right of "Launch", and click "Edit".
1. Change the revision number to match git repository tag of desired release version.

If there are any Nextflow secrets to update and expose to the pipeline, those secrets must be updated and/or added to the Seqera Platform workspace secrets.

#### Cleaning up AWS secrets

As of Seqera Platform v24.2.0 (Aug 20 2024), there is a bug where AWS secrets are not deleted (as described in the documentation) after pipeline termination (after a successful, cancelled, or failed run). For secret name `SECRETNAME` a run with ID `XYZ`, Seqera Platform writes a corresponding secret name `tower-XYZ/SECRETNAME` to AWS Secrets Manager.

For convenience, we've included a script [`scripts/cleanup_secrets.sh`](scripts/cleanup_secrets.sh) that can be used to delete AWS secrets provisioned for this pipeline
(matching `tower-*/TBCHASIN_*`).

Running the following will prompt you to confirm deletion of each matching secret found in AWS Secrets Manager.
```console
bash scripts/cleanup_secrets.sh [-a aws_profile] [-r region]
```

Be careful to not delete secrets used for pipeline runs that are currently running and in progress.