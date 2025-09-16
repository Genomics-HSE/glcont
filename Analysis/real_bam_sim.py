#!/usr/bin/env python3
import argparse, random, math, sys
from collections import defaultdict
import pysam

def load_genome_size_from_fai(fai_path):
    size = 0
    with open(fai_path) as f:
        for line in f:
            if not line.strip(): continue
            toks = line.split('\t')
            size += int(toks[1])
    return size

def estimate_read_length(bam, nmax=20000):
    n = 0
    total = 0
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped: continue
        if r.is_secondary or r.is_supplementary: continue
        alen = r.query_alignment_length
        if alen and alen > 0:
            total += alen
            n += 1
            if n >= nmax: break
    if n == 0:
        raise RuntimeError("Не удалось оценить длину чтения (нет выровненных чтений). Укажите --read-length вручную.")
    return total / n

def count_primary_mapped_reads(bam):
    cnt = 0
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped: continue
        if r.is_secondary or r.is_supplementary: continue
        # дубликаты обычно лучше исключить при оценке покрытия
        if r.is_duplicate: continue
        cnt += 1
    return cnt

def build_output_header(bam1, bam2):
    header = bam1.header.to_dict()

    # Убедимся, что есть SQ из bam1 (предпочтительно согласованные с рефом).
    # Добавим RG для двух источников.
    rg_list = header.get('RG', [])
    rg_ids = set([rg.get('ID') for rg in rg_list if 'ID' in rg])

    def add_rg(rgid, sample):
        if rgid in rg_ids:
            return
        rg = {'ID': rgid, 'SM': sample}
        rg_list.append(rg)
        rg_ids.add(rgid)

    add_rg('SRC1', 'mix_src1')
    add_rg('SRC2', 'mix_src2')
    header['RG'] = rg_list

    return pysam.AlignmentHeader.from_dict(header)

def write_with_rg(out, read, rgid):
    # Проставим RG, сохранив существующие теги
    try:
        read.set_tag('RG', rgid, value_type='Z')
    except Exception:
        # если запись ридов защищена, создадим копию
        read = read.to_string()
        read = pysam.AlignedSegment.fromstring(read, out.header)
        read.set_tag('RG', rgid, value_type='Z')
    out.write(read)

def stream_subsample_and_write(bam, out, target_reads, keep_prob, rgid, paired_policy_cache):
    """
    Возвращает фактически записанное число чтений (алигнментов).
    keep_prob — вероятность принять *кортеж* с данным QNAME (одинакова для всех фрагментов с этим именем).
    paired_policy_cache — общий кэш решений по QNAME (чтобы решения были согласованы).
    """
    written = 0
    for r in bam.fetch(until_eof=True):
        if r.is_secondary or r.is_supplementary:
            # вторичные/добавочные можно пропускать, чтобы не перекашивать покрытие
            continue

        qn = r.query_name
        keep = paired_policy_cache.get(qn)
        if keep is None:
            # единое решение для всех фрагментов с этим именем
            keep = (random.random() < keep_prob)
            paired_policy_cache[qn] = keep

        if keep and written < target_reads:
            write_with_rg(out, r, rgid)
            written += 1

        # По желанию можно очищать кэш, но не знаем, встретится ли парный фрагмент позже.
        if written >= target_reads:
            # Можно не останавливаться, чтобы "добрать пары", но мы таргетируем ровно число чтений.
            # Для строгого парного отбора можно собирать по имени и писать обе части вместе.
            break
    return written

def main():
    ap = argparse.ArgumentParser(
        description="Смешивает два BAM в заданной пропорции и с целевой средней глубиной."
    )
    ap.add_argument("--bam1", required=True, help="Первый входной BAM (может быть несортированным).")
    ap.add_argument("--bam2", required=True, help="Второй входной BAM.")
    ap.add_argument("--out-bam", required=True, help="Выходной несортированный BAM.")
    ap.add_argument("--proportion", type=float, required=True,
                    help="Доля первого BAM в финальном наборе чтений (0..1). Например, 0.3 означает 30%% из bam1 и 70%% из bam2.")
    ap.add_argument("--target-depth", type=float, required=True,
                    help="Целевая средняя глубина по геному (например, 30 для ~30x).")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--genome-size", type=int, help="Размер референса (сумма длин хромосом).")
    g.add_argument("--fai", help="FASTA index (.fai) для вычисления размера генома.")
    ap.add_argument("--read-length", type=int, default=0,
                    help="Средняя длина чтения (если 0 — оценим автоматически). Для парных — длина одного рида.")
    ap.add_argument("--seed", type=int, default=1, help="Сид генератора случайных чисел.")
    args = ap.parse_args()

    if not (0.0 <= args.proportion <= 1.0):
        sys.exit("proportion должен быть в [0,1].")

    random.seed(args.seed)

    # Открываем входы
    bam1 = pysam.AlignmentFile(args.bam1, "rb", check_sq=False)
    bam2 = pysam.AlignmentFile(args.bam2, "rb", check_sq=False)

    # Размер генома
    if args.genome_size:
        genome_size = args.genome_size
    else:
        genome_size = load_genome_size_from_fai(args.fai)

    # Оценим среднюю длину чтения при необходимости
    read_len = args.read_length if args.read_length > 0 else int(round(estimate_read_length(pysam.AlignmentFile(args.bam1, "rb", check_sq=False))))

    # Оценим число доступных "годных" чтений в каждом BAM
    # (первичные, выровненные, не-дубликаты)
    N1_total = count_primary_mapped_reads(pysam.AlignmentFile(args.bam1, "rb", check_sq=False))
    N2_total = count_primary_mapped_reads(pysam.AlignmentFile(args.bam2, "rb", check_sq=False))

    # Сколько чтений нужно для целевой глубины?
    # depth ~ (N_reads * read_len) / genome_size  =>  N_reads = depth * genome_size / read_len
    N_target = int(math.ceil(args.target_depth * genome_size / read_len))

    # Распределим по пропорции
    N1_target = int(round(args.proportion * N_target))
    N2_target = N_target - N1_target

    if N1_total == 0 or N2_total == 0:
        sys.exit("В одном из BAM нет пригодных чтений (проверьте мэппинг/флаги).")

    # Вероятности отбора из каждого BAM
    s1 = min(1.0, N1_target / float(N1_total))
    s2 = min(1.0, N2_target / float(N2_total))

    # Заголовок и выход
    out_header = build_output_header(bam1, bam2)
    out = pysam.AlignmentFile(args.out_bam, "wb", header=out_header)

    # Общий кэш решений по QNAME, чтобы парные/мульти-фрагментные чтения принимались/отбрасывались вместе
    policy_cache_1 = {}
    policy_cache_2 = {}

    # Потоковое выборочное копирование
    written1 = stream_subsample_and_write(pysam.AlignmentFile(args.bam1, "rb", check_sq=False),
                                          out, N1_target, s1, "SRC1", policy_cache_1)
    written2 = stream_subsample_and_write(pysam.AlignmentFile(args.bam2, "rb", check_sq=False),
                                          out, N2_target, s2, "SRC2", policy_cache_2)

    out.close()

    # Отчёт в stderr
    sys.stderr.write(
        f"[mix_bams] genome_size={genome_size}, read_len≈{read_len}, target_depth={args.target_depth}x\n"
        f"[mix_bams] target_reads={N_target} => bam1:{N1_target} bam2:{N2_target}\n"
        f"[mix_bams] available_reads bam1:{N1_total} bam2:{N2_total}, keep_prob bam1:{s1:.4f} bam2:{s2:.4f}\n"
        f"[mix_bams] written bam1:{written1} bam2:{written2} total:{written1+written2}\n"
        f"[mix_bams] Готово. Не забудьте отсортировать и проиндексировать выходной BAM.\n"
    )

if __name__ == "__main__":
    main()