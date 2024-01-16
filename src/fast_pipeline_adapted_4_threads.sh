#! /bin/bash

# Usage: ./fast_pipeline.sh <projectname> <reference.fasta> <input.fq> <out_dir>

#-------------------------------------------------------------------------------
### check commandline
#-------------------------------------------------------------------------------
if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ] || [ -z "$4" ]; then
  echo "Pipeline"
  echo "--------"
  echo "Usage: ./fast_pipeline_adapted_single_thread.sh <projectname> <reference.fasta> <input.fq> <out_dir>"
  echo ""
  exit
fi

project_name=$1
reference_genome=$2
input_fq=$3
out_dir=$4
current_dir="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

echo "  i : aligning reads to reference genome for $project_name .."
bwa aln -f "${out_dir}"/"${project_name}".aln -t 4 "${reference_genome}" "${input_fq}"

echo "  i : creating SAM file for $project_name .."
bwa samse -f "${out_dir}"/"${project_name}".sam "${reference_genome}" "${out_dir}"/"${project_name}".aln "${input_fq}"
rm "${out_dir}"/"${project_name}".aln

echo "  i : creating results table for $project_name .."
chrom_arm_py_script="${current_dir}"/faster_seq_analysis_chr_arm.py
# check existence of the chromosome arm counting script
if [ ! "${chrom_arm_py_script}" ]; then
  echo "Python script for counting fragments on chromosome arms not found! Expected location: '${chrom_arm_py_script}'"
  exit
fi

python3 "${chrom_arm_py_script}" "${out_dir}"/"${project_name}".sam \
"${out_dir}"/"${project_name}".chrarm.txt "${current_dir}"/chrom_arm_hg19.txt