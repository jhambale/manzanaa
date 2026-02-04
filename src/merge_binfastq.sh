#!/bin/bash

######################################################

# This script performs the following:
# combines specified fastqs
# outputs are written to directory where fastqs are stored


# Usage: bash merge_binfastq.sh \
# ../../manzanaa_data/ngs_raw/ex_001/demultiplex/ex_merged_R1_001.fastq \
# ../../manzanaa_data/ngs_raw/ex_001/demultiplex/*R1_001.fastq.gz \

######################################################


# Simple FASTQ merger
# Usage: ./simple_merge.sh output.fastq.gz *.fastq.gz ...

if [ $# -lt 2 ]; then
    echo "Usage: $0 output.fastq.gz *.fastq.gz ..."
    exit 1
fi

OUTPUT="$1"
shift  # Remove first argument (output file)

# Clear/create output file
> "$OUTPUT"

# Concatenate all input files
for file in "$@"; do
    if [ -f "$file" ]; then
        echo "Merging: $file"
        zcat "$file" >> "$OUTPUT"
    else
        echo "Warning: $file not found, skipping..."
    fi
done

# gzip > "$OUTPUT"

echo "Merge complete! Output saved to: $OUTPUT"
