# tb-bfx

Tools and pipelines for various off-target detection assays:

- CHANGE-Seq
- Cryptic-Seq

<!---toc start-->
   * [Local Setup](#local-setup)
      * [Running Tests](#running-tests)
   * [Run the Pipeline](#run-the-pipeline)
      * [Reference Preparation](#reference-preparation)
      * [Execution](#execution)

<!---toc end-->

## Local Setup

- [Install conda][conda-link]

- Install mamba

```console
conda install -c conda-forge mamba
```

- Get a local copy of the `tb-bfx` repo

```console
git clone git@github.com:tomebio/tb-bfx.git
cd tb-bfx
```

- Create the `tb-bfx` conda environment


```console
mamba env create -f environment.yml
```

- Activate the `tb-bfx` conda environment

```console
mamba activate tb-bfx
```


- Install `pytomebio` (developer mode)


```console
python setup.py develop
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

An example `config.yml` is shown below:


```yaml
TBD
```

The FASTQs for a given group are assumed to be in the given FASTQ directory, and named `<sample-name>_R1_001.fastq.gz` and `<sample-name>_R2_001.fastq.gz`.

See [Reference Preparation](#reference-preparation) for how to prepare the reference genome.

The config is organized at two levels: run and group level.
The run level contains configuration that applies to all groups and samples, for example the path to the `digenomitas` JAR.
The group level contains configuration that applies to related samples, for example the guide for specific CRISPR samples, or tool-specific parameters recommended for integrase samples.

| Config Key                | Description                                                                                        | Level  | Required | Default |
|---------------------------|----------------------------------------------------------------------------------------------------|--------|----------|---------|
| `fq_dir`                  | The directory containing FASTQs with suffixes `_R<#>_001.fastq.gz`                                 | Group  | Yes      | NA      |
| `ref_fasta`               | The path to the reference FASTA, with accompanying BWA index files and FASTA index                 | Group  | Yes      | NA      |
| `name`                    | The name of the sample (FASTQs must be `<name>_R<1 or 2>_001.fastq.gz`                             | Group  | Yes      | NA      |


[conda-link]: https://docs.conda.io/projects/conda/en/latest/user-guide/install/
