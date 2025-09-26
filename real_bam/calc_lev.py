import Levenshtein as lev
import argparse
from pathlib import Path
from typing import Dict, List, Tuple
import concurrent.futures as futures
import sys

def read_first_fasta_sequence(fn, remove_gaps: bool = True, uppercase: bool = True) -> str:
    seq_lines = []
    with open(fn) as fh:
        in_seq = False
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if in_seq and seq_lines:
                    break
                in_seq = True
                continue
            if in_seq:
                seq_lines.append(line.strip())
    seq = "".join(seq_lines)
    if uppercase:
        seq = seq.upper()
    if remove_gaps:
        seq = seq.replace("-", "").replace(".", "")
    return seq

s1 = read_first_fasta_sequence("in1.fa")
s2 = read_first_fasta_sequence("in2.fa")

print(lev.distance(s1,s2))