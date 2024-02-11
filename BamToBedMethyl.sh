#!/bin/sh
#SBATCH --mem=16g
#SBATCH -c 8
#SBATCH --time=2:00:00
#SBATCH -o modkit-logs


error_exit() {
  echo "Error: $1" >&2
  exit 1
}

# Get the input file path from the first argument
input_file_path="$1"

# Check if input file path is provided
if [[ -z "$input_file_path" ]]; then
  error_exit "Input file path not provided. Usage: $0 <input_file_path>"
fi

# Extract file name and base path from input path
input_file_name=$(basename "$input_file_path")
input_base_path=$(dirname "$input_file_path")

# Set output file names based on input file name
output_fastq="$input_base_path/$input_file_name.fastq"
output_sam="$input_base_path/$input_file_name.fq_mapped.sam"
output_bam="$input_base_path/$input_file_name.fq_mapped.sam.bam"
output_sorted_bam="$input_base_path/$input_file_name.fq_mapped.sam_sorted.bam"
output_bed="$input_base_path/$input_file_name.fq_mapped.sam_sorted.bam.CpG.bed"

# Commands with error handling
samtools fastq -T MM,ML - < "$input_file_path" > "$output_fastq" || error_exit "samtools fastq failed"
minimap2 --MD -L -t 9 -ax map-ont -y /sci/data/reference_data/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa "$output_fastq" > "$output_sam" || error_exit "minimap2 failed"
samtools view -bS -o "$output_bam" "$output_sam" || error_exit "samtools view failed"
samtools sort "$output_bam" -o "$output_sorted_bam" || error_exit "samtools sort failed"
samtools index "$output_sorted_bam" || error_exit "samtools index failed"
resources/tools/modkit/dist/modkit pileup "$output_sorted_bam" "$output_bed" --ref /sci/data/reference_data/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa --preset traditional --only-tabs --threads 4 || error_exit "modkit pileup failed"

echo "Pipeline completed successfully!"
