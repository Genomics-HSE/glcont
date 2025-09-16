#!/usr/bin/env python3
"""
Filter pairs of sequences by distance range from a distance-matrix CSV.

Input CSV format:
  ,label1,label2,...
  label1,0,12,7,...
  label2,12,0,9,...
  ...

Usage:
  python filter_pairs_by_distance.py -i pairwise_levenshtein.csv -m 5 -n 20 -o pairs_5_20.csv

Options:
  -i / --input   Path to distance-matrix CSV
  -m / --min     Minimum distance (inclusive)
  -n / --max     Maximum distance (inclusive)
  -o / --out     Output CSV for pairs (default: stdout)
  --include-diagonal  Include i==j pairs (default: False)
  --allow-duplicates  Include both (A,B) and (B,A) (default: False)
"""

import argparse
import csv
import sys
from typing import List, Tuple
import numpy as np

def load_distance_matrix(csv_path: str) -> Tuple[List[str], np.ndarray]:
    labels: List[str] = []
    rows: List[List[int]] = []

    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        if header is None or len(header) < 2:
            raise ValueError("Malformed CSV: missing header with labels.")
        labels = header[1:]  # first cell is empty

        for r in reader:
            if not r:
                continue
            row_label = r[0]
            if row_label not in labels:
                # tolerate, but not required
                pass
            try:
                row_vals = [int(float(x)) for x in r[1:]]
            except ValueError as e:
                raise ValueError(f"Non-numeric distance in row starting with {row_label}: {e}")
            rows.append(row_vals)

    M = np.array(rows, dtype=np.int64)
    if M.shape[0] != len(labels) or M.shape[1] != len(labels):
        raise ValueError(f"Matrix shape {M.shape} does not match number of labels {len(labels)}.")

    return labels, M

def find_pairs(labels: List[str],
               M: np.ndarray,
               dmin: int,
               dmax: int,
               include_diagonal: bool = False,
               allow_duplicates: bool = False) -> List[Tuple[str, str, int]]:
    n = len(labels)
    pairs: List[Tuple[str, str, int]] = []

    if allow_duplicates:
        for i in range(n):
            j_start = 0 if include_diagonal else 0
            for j in range(j_start, n):
                if not include_diagonal and i == j:
                    continue
                d = int(M[i, j])
                if dmin <= d <= dmax:
                    pairs.append((labels[i], labels[j], d))
    else:
        # unique unordered pairs i<j (and optionally i==j)
        if include_diagonal:
            for i in range(n):
                d = int(M[i, i])
                if dmin <= d <= dmax:
                    pairs.append((labels[i], labels[i], d))
        for i in range(n):
            for j in range(i + 1, n):
                d = int(M[i, j])
                if dmin <= d <= dmax:
                    pairs.append((labels[i], labels[j], d))
    return pairs

def write_pairs(pairs: List[Tuple[str, str, int]], out_path: str = "-") -> None:
    out_fh = sys.stdout if out_path in ("-", "", None) else open(out_path, "w", encoding="utf-8", newline="")
    close = out_fh is not sys.stdout
    try:
        w = csv.writer(out_fh)
        w.writerow(["seq1", "seq2", "distance"])
        for a, b, d in pairs:
            w.writerow([a, b, d])
    finally:
        if close:
            out_fh.close()

def main():
    ap = argparse.ArgumentParser(description="Show pairs of FASTAs with distance in [m, n] from a distance-matrix CSV.")
    ap.add_argument("-i", "--input", required=True, help="Path to distance-matrix CSV.")
    ap.add_argument("-m", "--min", type=int, required=True, help="Minimum distance (inclusive).")
    ap.add_argument("-n", "--max", type=int, required=True, help="Maximum distance (inclusive).")
    ap.add_argument("-o", "--out", default="-", help="Output CSV path (default: stdout).")
    ap.add_argument("--include-diagonal", action="store_true", help="Include self-pairs (i==j).")
    ap.add_argument("--allow-duplicates", action="store_true", help="Include both (A,B) and (B,A).")
    args = ap.parse_args()

    if args.min > args.max:
        print("[ERROR] min cannot be greater than max.", file=sys.stderr)
        sys.exit(2)

    labels, M = load_distance_matrix(args.input)
    pairs = find_pairs(labels, M, args.min, args.max,
                       include_diagonal=args.include_diagonal,
                       allow_duplicates=args.allow_duplicates)
    write_pairs(pairs, args.out)

if __name__ == "__main__":
    main()
