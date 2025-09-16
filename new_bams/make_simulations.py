#!/usr/bin/env python3
import argparse
import math
import os
import random
import sys
import tempfile
import pysam

def iter_primary_mapped(bam):
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        yield r

def mean_aligned_length(bam_path, max_reads_for_mean=None):
    bam = pysam.AlignmentFile(bam_path, "rb")
    total = 0
    cnt = 0
    for i, r in enumerate(iter_primary_mapped(bam), start=1):
        # query_alignment_length — число выровненных оснований рида
        al = r.query_alignment_length
        if al and al > 0:
            total += al
            cnt += 1
        if max_reads_for_mean and cnt >= max_reads_for_mean:
            break
    bam.close()
    if cnt == 0:
        raise RuntimeError(f"No primary mapped reads found in {bam_path}.")
    return total / cnt, cnt

def reservoir_sample_bam(bam_path, k, seed=None):
    """
    Случайно выбирает k ридов (primary mapped) из BAM,
    возвращает список строк SAM (r.to_string()).
    Использует reservoir sampling, чтобы не держать всё в памяти.
    """
    if k <= 0:
        return []

    if seed is not None:
        random.seed(seed)

    bam = pysam.AlignmentFile(bam_path, "rb")
    reservoir = []
    i = 0
    for r in iter_primary_mapped(bam):
        i += 1
        s = r.to_string()  # сам-строка рида
        if i <= k:
            reservoir.append(s)
        else:
            j = random.randint(1, i)
            if j <= k:
                reservoir[j-1] = s
    bam.close()

    if i < k:
        # если в файле меньше ридов, чем нужно — вернём столько, сколько есть
        sys.stderr.write(f"[warn] {bam_path}: requested {k}, available {i}. Using {i}.\n")
    return reservoir[:min(k, i)]

def write_sam_records_to_bam(records, header_bam_path, out_bam_path):
    """
    Пишет список SAM-строк в BAM с заголовком, взятым из header_bam_path.
    """
    with pysam.AlignmentFile(header_bam_path, "rb") as hbam:
        header = hbam.header

        with pysam.AlignmentFile(out_bam_path, "wb", header=header) as outbam:
            for s in records:
                aln = pysam.AlignedSegment.fromstring(s, header)
                outbam.write(aln)

def main():
    ap = argparse.ArgumentParser(description="Mix two BAMs to target depth h with proportion p of reads from BAM1.")
    ap.add_argument("--bam1", required=True, help="Path to first BAM")
    ap.add_argument("--bam2", required=True, help="Path to second BAM")
    ap.add_argument("--out", required=True, help="Path to output (sorted) BAM")
    ap.add_argument("--h", type=float, required=True, help="Target mean depth (coverage)")
    ap.add_argument("--p", type=float, required=True, help="Proportion of reads from BAM1 (0..1)")
    ap.add_argument("--genome-length", type=int, default=16569, help="Genome length L (default 16569)")
    ap.add_argument("--seed", type=int, default=None, help="Random seed (optional)")
    ap.add_argument("--max-reads-for-mean", type=int, default=None,
                    help="Limit reads for mean length estimation (speed-up); default: use all")
    args = ap.parse_args()

    if not (0.0 <= args.p <= 1.0):
        args.p /= 100.0  # если ввели в процентах, переведём в долю
    if args.h <= 0:
        ap.error("h must be > 0.")
    if args.genome_length <= 0:
        ap.error("genome-length must be > 0.")

    if args.seed is not None:
        random.seed(args.seed)

    # 1) Оценка средних выровненных длин
    l1, n1_est = mean_aligned_length(args.bam1, args.max_reads_for_mean)
    l2, n2_est = mean_aligned_length(args.bam2, args.max_reads_for_mean)

    sys.stderr.write(f"[info] mean aligned length bam1 = {l1:.2f} (n={n1_est} used)\n")
    sys.stderr.write(f"[info] mean aligned length bam2 = {l2:.2f} (n={n2_est} used)\n")

    # 2) Сколько ридов нужно, чтобы получить B = h*L покрытых оснований
    B = args.h * args.genome_length
    denom = args.p * l1 + (1.0 - args.p) * l2
    if denom <= 0:
        raise RuntimeError("Non-positive expected aligned length; check input BAMs.")

    N = math.ceil(B / denom)
    n1_target = int(round(args.p * N))
    n2_target = max(0, N - n1_target)

    sys.stderr.write(f"[info] target total reads N ≈ {N} (bam1: {n1_target}, bam2: {n2_target})\n")

    # 3) Сэмплинг ридов
    rec1 = reservoir_sample_bam(args.bam1, n1_target, seed=None if args.seed is None else args.seed + 1)
    rec2 = reservoir_sample_bam(args.bam2, n2_target, seed=None if args.seed is None else args.seed + 2)

    sys.stderr.write(f"[info] sampled bam1: {len(rec1)} reads, bam2: {len(rec2)} reads\n")

    # 4) Запись во временный BAM, сортировка и индексация
    with tempfile.TemporaryDirectory() as td:
        unsorted_bam = os.path.join(td, "mixed.unsorted.bam")
        write_sam_records_to_bam(rec1 + rec2, args.bam1, unsorted_bam)

        # сортируем и индексируем
        pysam.sort("-o", args.out, unsorted_bam)
        pysam.index(args.out)

    # 5) Оценка фактической средней глубины (грубая)
    # Посчитаем суммарную выровненную длину выбранных ридов и разделим на L
    sum_aligned = 0
    for s in (rec1 + rec2):
        # Быстро вытащим поле CIGAR из SAM-строки и посчитаем M/I/D/=/X как выровненную длину по query.
        # Надёжнее — снова распарсить через pysam, но это уже сделано при записи.
        # Поэтому подсчитаем через повторное чтение итогового BAM (надёжно).
        pass

    # Надёжная проверка через чтение отсортированного BAM:
    sum_aligned = 0
    with pysam.AlignmentFile(args.out, "rb") as obam:
        for r in iter_primary_mapped(obam):
            if r.query_alignment_length and r.query_alignment_length > 0:
                sum_aligned += r.query_alignment_length
    est_depth = sum_aligned / args.genome_length
    sys.stderr.write(f"[info] crude realized depth ≈ {est_depth:.2f}x\n")

if __name__ == "__main__":
    main()
