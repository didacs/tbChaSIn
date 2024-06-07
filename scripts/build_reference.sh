#!/usr/bin/env bash
# TODO: add support for curl as an alternative to wget

set -exo pipefail

usage() {
  cat <<EOF
Usage: $0 [-i <uri> | -g <genome>] [options] attB.fasta attP.fasta plasmid.fasta

Build a reference genome for a given integrase. The integrase-specific FASTA files are given as
positional arguments. At least three integrase-specific FASTAs are expected unless the '-f' option
is specified.

Optional:
    -g NAME  Name of the genome to use. Defaults to the folder name in the URI. If '-i' is not specified then
             this is assumed to be an NCBI reference and the URI is inferred from the name.
    -i URI   The URI of the folder containing the genomic reference. Defaults to NCBI GRCh38.p14.
    -n NAME  Name of the reference to generate. Defaults to the short genome name (everything after the last '_').
    -c CHRS  Comma-separted list of contigs to include from the reference. Defaults to the full genomic reference.
    -D DIR   Root directory to use for storing genomes and references. Defaults to the system
             temporary directory. Not required if both '-g' and '-r' are given.
    -G DIR   Directory where genomes are to be stored. Defaults to 'ROOT/genomes/'.
    -R DIR   Directory where references are to be stored. Defaults to 'ROOT/references/'.
    -f       Overwrite the reference even if it already exists.
    -x       Force building of the reference even if < 3 integrase-specific FASTAs are given.
    -A       Do not include alt contigs in the final reference.
EOF
  # shellcheck disable=SC2188
  >&2
  exit 1
}

NCBI_BASE_URI="https://ftp.ncbi.nlm.nih.gov/genomes/all"
DEFAULT_GENOME="GCF_000001405.40_GRCh38.p14"
GENOME_REGEX="(GC.)_([0-9]{3})([0-9]{3})([0-9]{3}).*_(.+)"

genome_folder_uri=""
ref_name=""
work_dir="/var/tmp"
genome_root=""
genome=""
chrs=""
ref_root=""
if [ -n "${TMPDIR}" ]; then
  work_dir="${TMPDIR}"
fi
force=false
no_check=false
no_alt=false
while getopts "i:g:c:n:D:G:R:fxA" flag; do
  case "${flag}" in
  i) genome_folder_uri=${OPTARG} ;;
  g) genome=${OPTARG} ;;
  c) chrs=${OPTARG} ;;
  n) ref_name=${OPTARG} ;;
  D) work_dir=${OPTARG} ;;
  G) genome_root=${OPTARG} ;;
  R) ref_root=${OPTARG} ;;
  f) force=true ;;
  x) no_check=true ;;
  A) no_alt=true ;;
  *) usage ;;
  esac
done
shift $((OPTIND - 1))

# treat positional arguments as integrase-specific fasta files
int_fastas=("$@")
if [ ! ${no_check} ] && [ "${#int_fastas[@]}" -lt 3 ]; then
  usage "Expected at least three FASTA files for the integrase"
fi

if [ -z "${genome_folder_uri}" ] && [ -z "${genome}" ]; then
  genome="${DEFAULT_GENOME}"
fi

if [ -z "${genome}" ]; then
  genome="$(basename "${genome_folder_uri}")"
elif [ -z "${genome_folder_uri}" ]; then
  if [[ ${genome} =~ ${GENOME_REGEX} ]]; then
    genome_folder_uri="${NCBI_BASE_URI}/${BASH_REMATCH[1]}/${BASH_REMATCH[2]}/${BASH_REMATCH[3]}/${BASH_REMATCH[4]}/${genome}"
  else
    usage "Genome ${genome} does not match NCBI format"
  fi
fi

if [ -z "${ref_name}" ]; then
  if [[ ${genome} =~ ${GENOME_REGEX} ]]; then
    ref_name="${BASH_REMATCH[5]}"
  else
    ref_name="${genome}"
  fi
fi

if [ -z "${ref_root}" ]; then
  ref_root="${work_dir}/references"
fi

ref_dir="${ref_root}/${ref_name}"
echo "Creating reference in directory ${ref_dir}"
if [ -d "${ref_dir}" ]; then
  if [ ${force} ]; then
    echo "Overwriting reference directory"
    rm -rf "${ref_dir}"
  else
    usage "Reference directory already exists"
  fi
fi
mkdir -p "${ref_dir}"

if [ -z "${genome_root}" ]; then
  genome_root="${work_dir}/genomes"
fi

genome_dir="${genome_root}/${genome}"
echo "Downloading genome files to ${genome_dir}"
if [ -d "${genome_dir}" ]; then
  echo "Genome directory already exists"
else
  mkdir -p "${genome_dir}"
fi

# download a file to $genome_dir
# 1. source folder URI
# 2. target folder
# 3. file name
download_genome_file() {
  folder_uri="${1}"
  outdir="${2}"
  filename="${3}"
  file_uri="${folder_uri}/${filename}"
  outfile="${outdir}/${filename}"
  # for some reason wget exits with code 4 even after successfully downloading
  set +e
  wget -P "${outdir}" -F "${filename}" "${file_uri}"
  rc=$?
  if [ $rc -ne 0 ]; then
    if [ $rc -ne 4 ]; then
      echo "Error downloading genome file ${file_uri}"
      exit $rc
    elif [ ! -f "${outfile}" ]; then
      echo "Error downloading genome file ${file_uri}"
      exit 4
    fi
  fi
  set -e
}

chr_dir="${genome_dir}/chr"
chr_temp_file="${chr_dir}/temp.fna"
if [ -n "${chrs}" ]; then
  chr_folder_uri="${genome_folder_uri}/${genome}_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA"
  echo "Downloading chromosome files to ${chr_dir}"
  if [ -d "${chr_dir}" ]; then
    echo "Chromosome directory already exists"
  else
    mkdir -p "${chr_dir}"
  fi
  IFS=',' read -ra chrs_array <<<"${chrs}"
  chr_file_array=()
  for chr in "${chrs_array[@]}"; do
    chr_filename="${chr}.fna.gz"
    chr_file="${chr_dir}/${chr_filename}"
    if [ ! -f "${chr_file}" ]; then
      download_genome_file "${chr_folder_uri}" "${chr_dir}" "${chr_filename}"
    fi
    chr_file_array+=("${chr_file}")
  done
  genome_file="${chr_temp_file}"
  zcat "${chr_file_array[@]}" >"${genome_file}"
else
  genome_filename="${genome}_genomic.fna.gz"
  genome_unzipped="${genome_dir}/${genome}_genomic.fna"
  genome_file="${genome_dir}/${genome_filename}"
  if [ -f "${genome_unzipped}" ]; then
    echo "Unzipped genome file already exists, skipping download"
  else
    if [ -f "${genome_file}" ]; then
      echo "Genome file already exists, skipping download"
    else
      download_genome_file "${genome_folder_uri}" "${genome_dir}" "${genome_filename}"
    fi
    gunzip -f "${genome_file}"
  fi
  genome_file="${genome_unzipped}"
fi

report_filename="${genome}_assembly_report.txt"
report_file="${genome_dir}/${report_filename}"
if [ -f "${report_file}" ]; then
  echo "Genome report file already exists, skipping download"
else
  download_genome_file "${genome_folder_uri}" "${genome_dir}" "${report_filename}"
fi

gtf_filename="${genome}_genomic.gtf.gz"
gtf_unzipped="${genome_dir}/${genome}_genomic.gtf"
gtf_file="${genome_dir}/${gtf_filename}"
if [ -f "${gtf_unzipped}" ]; then
  echo "Unzipped GTF file already exists, skipping download"
else
  if [ -f "${gtf_file}" ]; then
    echo "Genome GTF file already exists, skipping download"
  else
    download_genome_file "${genome_folder_uri}" "${genome_dir}" "${gtf_filename}"
  fi
  gunzip -f "${gtf_file}"
fi
gtf_file="${gtf_unzipped}"

# The fasta file needs to be indexed
if [ -f "${genome_file}.fai" ] && [ ! "${genome_file}.fai" -ot "${genome_file}" ]; then
  echo "FASTA index already exists"
else
  samtools faidx "${genome_file}"
fi

# Create a sequence dictionary (.dict) to sort and name the contigs
dict_file="${genome_dir}/${genome}_assembly_report.dict"
if [ ! -f "${dict_file}" ]; then
  alt_opt=""
  if [ ! ${no_alt} ]; then
    alt_opt="AltScaffold"
  fi
  fgbio CollectAlternateContigNames \
    -i "${report_file}" \
    -o "${dict_file}" \
    -p UcscName \
    -a SequenceName AssignedMolecule GenBankAccession RefSeqAccession \
    -s AssembledMolecule UnlocalizedScaffold UnplacedScaffold ${alt_opt}
fi

# Update the contig names and sort the contigs in the FASTA
fixed_file="${ref_dir}/${ref_name}_genomic.fasta"
fgbio UpdateFastaContigNames \
  -i "${genome_file}" \
  -d "${dict_file}" \
  -o "${fixed_file}" \
  --sort-by-dict \
  --skip-missing

# Build the final fasta
ref_file="${ref_dir}/${ref_name}.fasta.gz"
cat "${fixed_file}" "${int_fastas[@]}" |
  seqkit seq -w 60 - |
  bgzip >"${ref_file}"

# Index the final fasta
ref_index="${ref_dir}/${ref_name}.fasta.gz.fai"
samtools faidx "${ref_file}"
if [ ! -f "${ref_index}" ]; then
  echo "FASTA index ${ref_index} was not created"
  exit 1
fi

# Create dict for final fasta
ref_dict="${ref_dir}/${ref_name}.dict"
samtools dict "${ref_file}" >"${ref_dict}"

# Build BWA index for final fasta
bwa index "${ref_file}"

# Remove temporary files
rm -f "${fixed_file}"
if [ -f "${chr_temp_file}" ]; then
  rm -f "${chr_temp_file}"
fi

# Update the contig names in the GTF
ref_gtf="${ref_dir}/${ref_name}.gtf"
fgbio UpdateGffContigNames \
  -i "${gtf_file}" \
  -d "${dict_file}" \
  -o /dev/stdout \
  --skip-missing |
  bedtools sort -i - >"${ref_gtf}"
