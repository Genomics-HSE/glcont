#!/usr/bin/env python3
"""
Build a pairwise Levenshtein distance matrix across consensus FASTA files,
using the python-Levenshtein library for fast distance computation.

Default layout assumed:
./TSI/consensuses/*.fa|*.fasta|*.fna
./IBS/consensuses/*.fa|*.fasta|*.fna
./HAN/consensuses/*.fa|*.fasta|*.fna

Outputs:
- distances.csv             : raw Levenshtein distances (integers)
- normalized_distances.csv  : distance / max(len(seq_i), len(seq_j)) in [0,1]

Usage (from the directory that contains TSI, IBS, HAN):
    python3 consensus_levenshtein_lib.py

Custom root or folders:
    python3 consensus_levenshtein_lib.py --root /path/to/root --folders TSI IBS HAN

Install dependency (if missing):
    pip install python-Levenshtein pandas
"""

import argparse
from pathlib import Path
from typing import Dict, List, Tuple
import concurrent.futures as futures
import sys

try:
    import Levenshtein as lev  # python-Levenshtein
except Exception as e:
    sys.stderr.write(
        "ERROR: Could not import 'Levenshtein'. Please install it first:\n"
        "    pip install python-Levenshtein\n"
        f"Original exception: {e}\n"
    )
    raise SystemExit(1)

try:
    import pandas as pd
except ImportError:
    pd = None


def read_first_fasta_sequence(path: Path, remove_gaps: bool = True, uppercase: bool = True) -> str:
    seq_lines = []
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
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


def collect_sequences(root: Path, folders: List[str], consensus_dir: str,
                      exts: Tuple[str, ...]) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    for pop in folders:
        base = root / pop / consensus_dir
        if not base.exists():
            continue
        for ext in exts:
            for fp in base.glob(f"*{ext}"):
                label = f"{pop}:{fp.stem}"
                seq = read_first_fasta_sequence(fp, remove_gaps=True, uppercase=True)
                if not seq:
                    continue
                seqs[label] = seq
    return seqs


def pair_indices(n: int):
    for i in range(n):
        for j in range(i, n):
            yield i, j


def compute_matrix(labels: List[str], seqs: List[str], threads: int = 1):
    n = len(labels)
    mat = np.zeros((n, n))

    def job(i_j):
        i, j = i_j
        if i == j:
            return i, j, 0
        d = lev.distance(seqs[i], seqs[j])
        return i, j, d

    if threads > 1:
        with futures.ThreadPoolExecutor(max_workers=threads) as ex:
            for i, j, d in ex.map(job, pair_indices(n)):
                mat[i][j] = d
                mat[j][i] = d
    else:
        for i, j in pair_indices(n):
            d = 0 if i == j else lev.distance(seqs[i], seqs[j])
            mat[i][j] = d
            mat[j][i] = d

    return mat


def normalize_matrix(mat: List[List[int]], seqs: List[str]):
    n = len(seqs)
    norm = [[0.0] * n for _ in range(n)]
    lens = [len(s) for s in seqs]
    for i in range(n):
        for j in range(n):
            denom = max(lens[i], lens[j]) if max(lens[i], lens[j]) > 0 else 1
            norm[i][j] = mat[i][j] / denom
    return norm


def main():
    ap = argparse.ArgumentParser(description="Pairwise Levenshtein distances (python-Levenshtein) for consensus FASTA files.")
    ap.add_argument("--root", type=str, default=".", help="Root directory containing population folders (default: .)")
    ap.add_argument("--folders", nargs="+", default=["TSI", "IBS", "HAN"],
                    help="Population folders to scan (default: TSI IBS HAN)")
    ap.add_argument("--consensus-dir", type=str, default="consensuses",
                    help="Subdirectory name holding consensus FASTA files (default: consensuses)")
    ap.add_argument("--exts", nargs="+", default=[".fa", ".fasta", ".fna"],
                    help="FASTA file extensions to include (default: .fa .fasta .fna)")
    ap.add_argument("--threads", type=int, default=1, help="Number of worker threads (default: 1)")
    ap.add_argument("--out-prefix", type=str, default="",
                    help="Prefix for output files (default: none)")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    seqs_dict = collect_sequences(root, args.folders, args.consensus_dir, tuple(args.exts))

    if not seqs_dict:
        raise SystemExit(f"No FASTA files found under {root} in {args.folders}/*/{args.consensus_dir}/ with exts {args.exts}")

    labels = sorted(seqs_dict.keys())
    seqs = [seqs_dict[k] for k in labels]

    mat = compute_matrix(labels, seqs, threads=args.threads)
    norm = normalize_matrix(mat, seqs)

    out_pref = args.out_prefix
    if pd is None:
        def write_csv(path: Path, header: List[str], rows: List[List]):
            with path.open("w", encoding="utf-8") as fw:
                fw.write(",".join(["id"] + header) + "\n")
                for name, row in zip(header, rows):
                    fw.write(",".join([name] + [str(x) for x in row]) + "\n")

        (root / f"{out_pref}distances.csv").write_text("", encoding="utf-8")
        write_csv(root / f"{out_pref}distances.csv", labels, mat)
        write_csv(root / f"{out_pref}normalized_distances.csv", labels, norm)
        print(f"Wrote: {root / f'{out_pref}distances.csv'}")
        print(f"Wrote: {root / f'{out_pref}normalized_distances.csv'}")
    else:
        import pandas as pd  # type: ignore
        df = pd.DataFrame(mat, index=labels, columns=labels)
        dfn = pd.DataFrame(norm, index=labels, columns=labels)
        df.to_csv(root / f"{out_pref}distances.csv")
        dfn.to_csv(root / f"{out_pref}normalized_distances.csv")
        print(f"Wrote: {root / f'{out_pref}distances.csv'}")
        print(f"Wrote: {root / f'{out_pref}normalized_distances.csv'}")


if __name__ == "__main__":
    main()
