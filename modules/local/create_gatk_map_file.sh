#!/bin/bash

# Script to create a GATK map file for GenomicsDBImport from a directory of GVCF files
# Usage: ./create_gatk_map_file.sh <gvcf_directory> <interval> [output_file]

# Check arguments
if [ $# -lt 2 ]; then
    echo "Usage: $0 <gvcf_directory> <interval> [output_file]"
    echo "Example: $0 ./results/gatk_haplotype 1 interval_1_map.txt"
    echo "If output_file is not specified, it will default to <interval>_sample_map.txt"
    exit 1
fi

# Directory containing the GVCF files
GVCF_DIR="$1"
INTERVAL="$2"
OUTPUT_FILE="${3:-${INTERVAL}_sample_map.txt}"

# Check if directory exists
if [ ! -d "$GVCF_DIR" ]; then
    echo "Error: Directory $GVCF_DIR does not exist"
    exit 1
fi

# Find all GVCF files for the specified interval
echo "Searching for GVCF files in $GVCF_DIR for interval $INTERVAL..."

# Initialize counter and map file
count=0
> "$OUTPUT_FILE"

# Search for GVCF files with the specific interval pattern
for gvcf_file in "$GVCF_DIR"/*".$INTERVAL.g.vcf.gz"; do
    if [ -f "$gvcf_file" ]; then
        # Extract sample name from filename (remove everything from the interval onwards)
        base_name=$(basename "$gvcf_file")
        # This pattern assumes filenames like "samplename.interval.g.vcf.gz"
        sample_id=${base_name%".$INTERVAL.g.vcf.gz"}
        
        # Write sample ID and GVCF path to map file (tab-delimited)
        echo -e "$sample_id\t$gvcf_file" >> "$OUTPUT_FILE"
        count=$((count + 1))
    fi
done

if [ $count -eq 0 ]; then
    echo "Warning: No GVCF files found matching *.$INTERVAL.g.vcf.gz in $GVCF_DIR"
    exit 1
else
    echo "Created $OUTPUT_FILE with $count samples for interval $INTERVAL"
    echo "Map file location: $(readlink -f $OUTPUT_FILE)"
fi

# Make the script executable
chmod +x "$0"

