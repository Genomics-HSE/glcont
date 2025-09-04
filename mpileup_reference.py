#!/usr/bin/env python3
import argparse
import gzip
import sys
from collections import defaultdict

ASCII_OFFSET = 33

def open_maybe_gzip(path, mode="rt"):
    if path == "-" or path is None:
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)

def read_fasta(path):
    fasta = {}
    with open_maybe_gzip(path, "rt") as f:
        name = None
        parts = []
        for line in f:
            if not line:
                continue
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    fasta[name] = "".join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if name is not None:
            fasta[name] = "".join(parts).upper()
    return fasta

def phred_weight(qchar, minQ):
    Q = ord(qchar) - ASCII_OFFSET
    # if Q < minQ:
        # return 0.0
    # expected probability of correctness as weight
    # return 1
    return 1.0 - 10.0 ** (-Q / 10.0)

def parse_mpileup_bases(bases_str, quals_str, ref_base, minQ):
    """
    Parse mpileup 'read bases' and 'base qualities' column for a single position.

    Returns:
      base_weights: dict base->float (A/C/G/T/N) including ref via '.' or ','
      ins_weights: dict ins_seq->float (sequence inserted AFTER this position)
      del_weights: dict del_seq->float (sequence deleted starting AFTER this position)
      total_weight: sum of all base weights that actually voted (for fraction thresholds)
    """
    base_weights = defaultdict(float)
    ins_weights = defaultdict(float)
    del_weights = defaultdict(float)

    i = 0
    qi = 0
    last_w = 0.0  # weight of last *anchoring* base (used for indels)
    ref_base = ref_base.upper()

    while i < len(bases_str):
        ch = bases_str[i]

        # Start of read: '^' then a single mapping-quality char (skip both)
        if ch == '^':
            i += 2
            continue

        # End of read marker
        if ch == '$':
            i += 1
            continue

        # Indels (use last_w as the support for the event)
        if ch in "+-":
            sign = ch
            i += 1
            j = i
            # length may be multiple digits
            while j < len(bases_str) and bases_str[j].isdigit():
                j += 1
            if j == i:
                # malformed; skip
                i += 1
                continue
            try:
                n = int(bases_str[i:j])
            except ValueError:
                i = j
                continue
            i = j
            seq = bases_str[i:i + n]
            i += n
            if n > 0 and seq:
                seqU = seq.upper()
                if sign == '+':
                    ins_weights[seqU] += last_w
                else:
                    del_weights[seqU] += last_w
            continue

        # Placeholder for a deleted base at this position
        if ch == '*':
            if qi < len(quals_str):
                _ = phred_weight(quals_str[qi], minQ)  # consume quality but do not vote
                qi += 1
            i += 1
            last_w = 0.0
            continue

        # Match to reference on fwd/rev strand
        if ch == '.' or ch == ',':
            if qi < len(quals_str):
                w = phred_weight(quals_str[qi], minQ)
                qi += 1
                if w > 0.0:
                    base_weights[ref_base] += w
                    last_w = w
                else:
                    last_w = 0.0
            i += 1
            continue

        # Explicit base
        if ch in 'ACGTNacgtn':
            if qi < len(quals_str):
                w = phred_weight(quals_str[qi], minQ)
                qi += 1
                if w > 0.0:
                    base_weights[ch.upper()] += w
                    last_w = w
                else:
                    last_w = 0.0
            i += 1
            continue

        # Unknown symbol; skip
        i += 1

    total_weight = sum(base_weights.values())
    return base_weights, ins_weights, del_weights, total_weight

def choose_consensus(ref_base, base_weights):
    """Pick the base with the largest weight; break ties in favor of ref."""
    ref_base = ref_base.upper()
    best_base = ref_base
    best_w = base_weights.get(ref_base, 0.0)
    sum_weights = sum(base_weights.values())
    for b, w in base_weights.items():
        if w > best_w or (w == best_w and b == ref_base):
            best_base = b
            best_w = w
    if best_w >= 0.5 * sum_weights:
        # If the best base is supported by at least half of the total weight, return it
        return best_base, best_w
    else:
        return 'N', 100
    return best_base, best_w

def best_indel(indel_dict):
    if not indel_dict:
        return None, 0.0
    seq, w = max(indel_dict.items(), key=lambda kv: kv[1])
    return seq, w

def consensus_from_mpileup(ref, mpileup_path, minQ=13, min_indel_frac=0.6, verbose=False):
    """
    Build consensus sequence(s) for contigs present in 'ref',
    using an mpileup file spanning zero or more of those contigs.
    """
    # Output per contig
    out = {ctg: [] for ctg in ref}
    # Position pointer and skip windows per contig
    pos_ptr = {ctg: 0 for ctg in ref}  # 0-based in reference
    skip_until = {ctg: -1 for ctg in ref}  # if >= current index, skip due to prior deletion

    def flush_to(ctg, target_pos):
        """Append untouched reference bases up to (but not including) target_pos,
        respecting any deletion-skips."""
        seq = ref[ctg]
        p = pos_ptr[ctg]
        while p < target_pos and p < len(seq):
            if p <= skip_until[ctg]:
                # This position is deleted in consensus → skip it
                p += 1
                continue
            out[ctg].append(seq[p])
            p += 1
        pos_ptr[ctg] = p

    with open_maybe_gzip(mpileup_path, "rt") as f:
        for line in f:
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                # mpileup minimal: CHR POS REF DEPTH BASES QUALS
                continue
            ctg, pos_s, ref_base, depth_s, read_bases, quals = cols[:6]
            if ctg not in ref:
                # Unrecognized contig in mpileup; skip
                continue
            # mpileup is 1-based
            pos1 = int(pos_s)
            idx0 = pos1 - 1
            # Move pointer to this position
            flush_to(ctg, idx0)

            # If this position is within a prior deletion, skip writing any base here
            if idx0 < len(ref[ctg]) and idx0 <= skip_until[ctg]:
                # Still consume pileup to maintain skip window decisions, but don't output
                # We won't extend skip here; deletion anchor is at its start only.
                pos_ptr[ctg] = idx0 + 1
                continue

            # Parse pileup bases
            base_w, ins_w, del_w, total_w = parse_mpileup_bases(read_bases, quals, ref_base, minQ)

            # Decide substitution (or keep ref)
            chosen_base, base_support = choose_consensus(ref_base, base_w)

            # Decide indels
            best_ins_seq, best_ins_w = best_indel(ins_w)
            best_del_seq, best_del_w = best_indel(del_w)

            ins_frac = (best_ins_w / total_w) if total_w > 0 else 0.0
            del_frac = (best_del_w / total_w) if total_w > 0 else 0.0

            # Write anchor base (substitution may alter it)
            if idx0 < len(ref[ctg]):
                out[ctg].append(chosen_base)
            pos_ptr[ctg] = idx0 + 1

            # Apply insertion if strongly supported
            if best_ins_seq and ins_frac >= min_indel_frac and best_ins_w >= base_support:
                out[ctg].append(best_ins_seq)

            # Apply deletion if strongly supported
            if best_del_seq and del_frac >= min_indel_frac and best_del_w >= base_support:
                # Deletion starts AFTER this position; remove the next len(del_seq) reference bases
                del_len = len(best_del_seq)
                skip_until[ctg] = max(skip_until[ctg], idx0 + del_len)

    # Flush any remaining tail for each contig
    for ctg in ref:
        flush_to(ctg, len(ref[ctg]))

    # Join
    return {ctg: "".join(seq_parts) for ctg, seq_parts in out.items()}

def write_fasta(consensus_dict, out_path, width=80):
    out = sys.stdout if out_path in (None, "-", "") else open(out_path, "wt")
    try:
        for ctg, seq in consensus_dict.items():
            print(f">{ctg}", file=out)
            for i in range(0, len(seq), width):
                print(seq[i:i+width], file=out)
    finally:
        if out is not sys.stdout:
            out.close()

def main():
    ap = argparse.ArgumentParser(
        description="Naive consensus caller from reference FASTA and samtools mpileup.")
    ap.add_argument("--ref", required=True, help="Reference FASTA (optionally .gz)")
    ap.add_argument("--mpileup", required=True, help="samtools mpileup file (optionally .gz)")
    ap.add_argument("-o", "--out", default="-", help="Output consensus FASTA (default: stdout)")
    ap.add_argument("--minQ", type=int, default=13,
                    help="Minimum base quality to contribute (default: 13)")
    ap.add_argument("--min-indel-frac", type=float, default=0.6,
                    help="Min fraction of total base weight to apply an indel (default: 0.6)")
    args = ap.parse_args()

    ref = read_fasta(args.ref)
    if not ref:
        sys.stderr.write("ERROR: no sequences found in reference FASTA\n")
        sys.exit(1)

    consensus = consensus_from_mpileup(
        ref,
        args.mpileup,
        minQ=args.minQ,
        min_indel_frac=args.min_indel_frac,
    )
    write_fasta(consensus, args.out)

if __name__ == "__main__":
    main()
