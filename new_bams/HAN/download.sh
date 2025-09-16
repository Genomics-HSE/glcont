#!/bin/bash

# File containing list of URLs (one per line)
URL_FILE="download_list.txt"

# Directory to save downloaded files
OUTDIR="mtdna_bams"

# Create output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Loop through each URL in the file
while IFS= read -r url; do
    # Skip empty lines
    [ -z "$url" ] && continue
    
    echo "Downloading: $url"
    
    # Use wget (you can switch to curl if you prefer)
    samtools view -T ../../rCRS.fa -b $url chrM > "$OUTDIR/$(basename "$url" .cram).bam"
    
    # Or with curl:
    # curl -L -O --output-dir "$OUTDIR" "$url"
done < "$URL_FILE"

echo "All downloads completed!"

