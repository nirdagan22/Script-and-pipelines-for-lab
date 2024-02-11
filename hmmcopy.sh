#!/bin/bash
#Script for using Hmmcopy libary


# Directory containing barcode subdirectories
base_dir="/sci/labs/aviadz/nirdagan/projects/minimap2-samples"

# Use find to list sorted BAM files recursively
find "$base_dir" -type f -name "*.bam" | while read -r bam_file; do

    # Create the output BAI file path in the same directory as BAM file
    bai_file="${bam_file%.bam}.bam.bai"

    # Create the .bam.bai file
    samtools index "$bam_file" "$bai_file"

    # Use the .bam and .bam.bai files for hmmcopy_ReadCount
    window=1000000
    quality=50
    chromosome="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22"
    output_wig="${bam_file%.bam}.wig"
    
    /sci/labs/aviadz/nirdagan/hmmcopy_utils/bin/readCounter --window "$window" --quality "$quality" --chromosome "$chromosome" "$bam_file" > "$output_wig"
done
