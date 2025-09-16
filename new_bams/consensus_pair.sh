#!/usr/bin/env bash
# Usage: ./pairs_to_consensus.sh pairs.csv

set -euo pipefail
PAIRS_CSV="$1"

# Skip header, parse CSV
tail -n +2 "$PAIRS_CSV" | while IFS=, read -r seq1 seq2 dist
do
    # Map "TSI__HG001" → "TSI/consensuses/HG001.fasta"
    cons1=$(echo "$seq1" | sed 's|__|/consensuses/|').fasta
    cons2=$(echo "$seq2" | sed 's|__|/consensuses/|').fasta

    echo "$cons1,$cons2"
done
