#!/usr/bin/env python3
# Usage: python parse_pairs.py pairs.csv

import sys

def parse_name(name: str):
    # Remove ".final" suffix if present
    if name.endswith(".final"):
        name = name[:-6]
    group, sample = name.split("__", 1)
    return f"{group},{sample}"
with open(f"{sys.argv[1][:-3]}samples.txt", "w") as out:
    with open(sys.argv[1], "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("seq1"):  # skip header if any
                continue
            a, b, dist = line.split(",")
            pair = f'{parse_name(a)},{parse_name(b)}'
            
            print(pair, file=out)