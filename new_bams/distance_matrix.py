#!/usr/bin/env python3
"""
Compute pairwise Levenshtein distances for consensus FASTA files found under:
  TSI/consensuses, IBS/consensuses, HAN/consensuses

Outputs a CSV matrix with filenames as row/column headers.

Requirements:
  pip install numpy python-Levenshtein
"""

import os
import sys
import gzip
import argparse
from typing import Dict, List, Tuple
import numpy as np
from Levenshtein import distance as lev_distance

FASTA_EXTS = {".fa", ".fasta", ".fna"}

def is_fasta(path: str) -> bool:
    base, ext = os.path.splitext(path)
    if ext == ".gz":
        base2, ext2 = os.path.splitext(base)
        return ext2.lower() in FASTA_EXTS
    return ext.lower() in FASTA_EXTS

def read_fasta_sequence(path: str) -> str:
    """Reads FIRST sequence in a (possibly gzipped) FASTA file and returns it as a single uppercase string."""
    opener = gzip.open if path.endswith(".gz") else open
    seq_lines: List[str] = []
    with opener(path, "rt") as fh:
        in_seq = False
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if in_seq and seq_lines:
                    break  # already captured first record
                in_seq = True
                continue
            if in_seq:
                seq_lines.append(line.strip())
    seq = "".join(seq_lines).upper()
    return seq

def find_fasta_files(root_dirs: List[str]) -> List[str]:
    files: List[str] = []
    for grp in root_dirs:
        consensuses_dir = os.path.join(grp, "consensuses")
        if not os.path.isdir(consensuses_dir):
            continue
        for dirpath, _, filenames in os.walk(consensuses_dir):
            for fn in filenames:
                path = os.path.join(dirpath, fn)
                if is_fasta(path):
                    files.append(os.path.abspath(path))
    files.sort()
    return files

def label_for(path: str) -> str:
    """Create a compact label: GROUP__filename (without extensions), keeping subgroup path if present."""
    # group is the top-level folder (TSI/IBS/HAN)
    parts = os.path.normpath(path).split(os.sep)
    group = parts[-3] if len(parts) >= 3 else parts[0]
    # filename without .gz and fasta extension
    fn = os.path.basename(path)
    if fn.endswith(".gz"):
        fn = fn[:-3]
    for ext in FASTA_EXTS:
        if fn.lower().endswith(ext):
            fn = fn[: -len(ext)]
            break
    return f"{group}__{fn}"

def build_sequences(paths: List[str]) -> Tuple[List[str], Dict[str, str]]:
    labels: List[str] = []
    seqs: Dict[str, str] = {}
    for p in paths:
        try:
            s = read_fasta_sequence(p)
            if not s:
                print(f"[WARN] Empty sequence in {p}, skipping.", file=sys.stderr)
                continue
            lab = label_for(p)
            labels.append(lab)
            seqs[lab] = s
        except Exception as e:
            print(f"[WARN] Failed to read {p}: {e}", file=sys.stderr)
    return labels, seqs

def compute_distance_matrix(labels: List[str], seqs: Dict[str, str]) -> np.ndarray:
    n = len(labels)
    M = np.zeros((n, n), dtype=np.int32)
    for i in range(n):
        si = seqs[labels[i]]
        M[i, i] = 0
        for j in range(i + 1, n):
            sj = seqs[labels[j]]
            d = lev_distance(si, sj)
            M[i, j] = d
            M[j, i] = d
    return M

def save_csv(out_path: str, labels: List[str], M: np.ndarray) -> None:
    with open(out_path, "w", encoding="utf-8") as f:
        # header
        f.write("," + ",".join(labels) + "\n")
        for lab, row in zip(labels, M):
            f.write(lab + "," + ",".join(str(int(x)) for x in row) + "\n")

def main():
    ap = argparse.ArgumentParser(description="Pairwise Levenshtein distances for consensus FASTA files.")
    ap.add_argument("-d", "--dirs", nargs="+", default=["TSI", "IBS", "HAN"],
                    help="Top-level folders to scan (default: TSI IBS HAN). Each must contain a 'consensuses' subfolder.")
    ap.add_argument("-o", "--out", default="pairwise_levenshtein.csv",
                    help="Output CSV file (default: pairwise_levenshtein.csv)")
    args = ap.parse_args()

    fasta_paths = find_fasta_files(args.dirs)
    if not fasta_paths:
        print("[ERROR] No FASTA files found in provided directories.", file=sys.stderr)
        sys.exit(1)

    labels, seqs = build_sequences(fasta_paths)
    if not labels:
        print("[ERROR] No valid sequences read.", file=sys.stderr)
        sys.exit(1)

    print(f"[INFO] Computing distances for {len(labels)} sequences...", file=sys.stderr)
    M = compute_distance_matrix(labels, seqs)
    save_csv(args.out, labels, M)
    print(f"[INFO] Done. Wrote {args.out}", file=sys.stderr)

if __name__ == "__main__":
    main()
