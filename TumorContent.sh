#!/bin/bash

#SBATCH --mem=16g
#SBATCH -n 12
#SBATCH -t 10:00:00
#SBATCH -o pipeline.log


# Check if the required arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <base_dir> <output_dir>"
    exit 1
fi

# Assign command-line arguments to variables
base_dir="$1"
output_dir="$2"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Use find to list FASTQ files recursively
find "$base_dir" -type f -name "*.fastq" | while read -r fastq_file; do
    # Print the name of the current FASTQ file
    #echo "Processing: $fastq_file"

    # Run minimap2 and convert to SAM
    sam_output="$output_dir/$(basename "$fastq_file" .fastq)_mapped.sam"
    minimap2 --MD -L -t 9 -ax map-ont -y /sci/data/reference_data/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa "$fastq_file" > "$sam_output"

    # Convert SAM to BAM
    bam_output="$output_dir/$(basename "$fastq_file" .fastq)_mapped.bam"
    samtools view -bS -o "$bam_output" "$sam_output"
done

# Merge and sort the BAM files directly
output_bam="$output_dir/all_samples_merged_sorted.bam"
find "$output_dir" -type f -name '*.bam' | xargs samtools cat | samtools sort -o "$output_bam"

# Create the .bam.bai file
bai_file="$output_bam.bai"
samtools index "$output_bam" "$bai_file"

# Use the .bam and .bam.bai files for hmmcopy_ReadCount
window=1000000
quality=50
chromosome="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"
output_wig="$bai_file.wig"

/sci/labs/aviadz/nirdagan/hmmcopy_utils/bin/readCounter --window "$window" --quality "$quality" --chromosome "$chromosome" "$output_bam" > "$output_wig"
