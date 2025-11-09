#!/usr/bin/env python3
"""
count.py

Exact-first counting with optional fuzzy fallback (k-mer prefilter + edlib).

Usage:
    python3 count.py -r reads.fasta -p peptides.fasta -o counts.csv [--tolerate-errors] [--threads N]

Notes:
 - By default only exact matching is performed.
 - Use --tolerate-errors to enable fuzzy-mode for reads that do not match exactly.
"""
import argparse
import os
import csv
import gzip
import tempfile
import sys
from typing import List, Dict, Tuple, Iterable
from collections import defaultdict, Counter
from tqdm import tqdm

# edlib used for fast edit-distance local alignment in fuzzy pass
try:
    import edlib
except Exception as e:
    edlib = None

# -----------------------
# I/O helpers
# -----------------------
def open_maybe_gz(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def read_fasta_stream(path: str) -> Iterable[Tuple[str, str]]:
    fh = open_maybe_gz(path)
    header = None
    seq_parts = []
    for line in fh:
        line = line.rstrip("\n\r")
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq_parts).upper()
            header = line[1:].split()[0]
            seq_parts = []
        else:
            seq_parts.append(line.strip())
    if header is not None:
        yield header, "".join(seq_parts).upper()
    fh.close()


def count_headers(path: str) -> int:
    cnt = 0
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    cnt += 1
    else:
        with open(path, "r") as fh:
            for line in fh:
                if line.startswith(">"):
                    cnt += 1
    return cnt


def basename_without_fasta(path: str) -> str:
    bn = os.path.basename(path)
    if bn.endswith(".gz"):
        bn = bn[:-3]
    for ext in (".fasta", ".fa", ".fna", ".ffn"):
        if bn.lower().endswith(ext):
            return bn[: -len(ext)]
    return os.path.splitext(bn)[0]


# -----------------------
# Load peptides & build indices
# -----------------------
def load_peptides(peptides_path: str) -> Tuple[List[str], List[str], Dict[str, List[int]]]:
    """
    Return (peptide_names, peptide_seqs, seq_to_indices).
    seq_to_indices maps exact peptide sequence -> list of indices (handles duplicate sequences).
    """
    names = []
    seqs = []
    seq_to_indices = {}
    idx = 0
    for header, seq in read_fasta_stream(peptides_path):
        names.append(header)
        seqs.append(seq.upper())
        seq_to_indices.setdefault(seq.upper(), []).append(idx)
        idx += 1
    return names, seqs, seq_to_indices


def build_kmer_index(seqs: List[str], k: int) -> Dict[str, List[int]]:
    """
    Build k-mer -> list(peptide_idx) index.
    We store each peptide's unique k-mers to avoid redundant entries.
    """
    idx = defaultdict(list)
    for i, s in enumerate(seqs):
        if len(s) < k:
            idx[s].append(i)
            continue
        seen = set()
        for j in range(len(s) - k + 1):
            kmer = s[j : j + k]
            if kmer in seen:
                continue
            seen.add(kmer)
            idx[kmer].append(i)
    return idx


# -----------------------
# Exact-first pass (writes unmatched reads to temp file)
# -----------------------
def exact_pass_and_write_unmatched(reads_path: str, seq_to_indices: Dict[str, List[int]], counts: List[int], tmp_unmatched_path: str):
    """
    Scan reads, update counts for exact matches. Writes unmatched reads to tmp_unmatched_path (FASTA).
    Returns (total_reads, matched_reads).
    """
    total_reads = 0
    matched_reads = 0

    # attempt to estimate reads for progress bar
    total = None
    try:
        total = count_headers(reads_path)
    except Exception:
        total = None

    with open(tmp_unmatched_path, "w") as unmatched_fh:
        it = read_fasta_stream(reads_path)
        if total:
            it = tqdm(it, total=total, unit="reads", desc="Exact pass")
        else:
            it = tqdm(it, unit="reads", desc="Exact pass")
        for hdr, seq in it:
            total_reads += 1
            if seq in seq_to_indices:
                matched_reads += 1
                for idx in seq_to_indices[seq]:
                    counts[idx] += 1
            else:
                # write unmatched to temp fasta
                unmatched_fh.write(f">{hdr}\n")
                for i in range(0, len(seq), 80):
                    unmatched_fh.write(seq[i : i + 80] + "\n")
        it.close()
    return total_reads, matched_reads


# -----------------------
# Fuzzy-pass helpers
# -----------------------
def adaptive_max_edits(length: int) -> int:
    """Adaptive threshold: 10% of length rounded down, min 1, cap 3."""
    return min(3, max(1, int(length * 0.1)))


def _fuzzy_worker(batch_reads: List[Tuple[str, str]], seqs: List[str], kmer_index: Dict[str, List[int]],
                  pid_to_len: List[int], top_n: int, k: int, min_shared_kmers: int = 2) -> List[Tuple[str, str, int, int, int]]:
    """
    Worker: returns a list of tuples for matched reads:
      (read_header, read_seq, best_pid, best_edit_distance, best_shared_kmer_count)
    """
    if edlib is None:
        raise RuntimeError("edlib is required for fuzzy matching. Install with `pip install edlib`")

    matched_records = []
    for hdr, seq in batch_reads:
        seq_u = seq.upper()
        L = len(seq_u)
        max_edits = adaptive_max_edits(L)

        # collect k-mer hit counts
        counter = Counter()
        if len(seq_u) < k:
            ids = kmer_index.get(seq_u, [])
            for pid in ids:
                counter[pid] += 1
        else:
            seen_kmers = set()
            for i in range(len(seq_u) - k + 1):
                kmer = seq_u[i : i + k]
                if kmer in seen_kmers:
                    continue
                seen_kmers.add(kmer)
                ids = kmer_index.get(kmer)
                if ids:
                    for pid in ids:
                        counter[pid] += 1

        if not counter:
            continue
        filtered = [(pid, score) for pid, score in counter.items() if score >= min_shared_kmers]
        if not filtered:
            continue

        candidates = []
        for pid, score in filtered:
            plen = pid_to_len[pid]
            if abs(plen - L) <= max_edits + 0:
                candidates.append((pid, score))

        if not candidates:
            continue

        candidates.sort(key=lambda x: -x[1])
        candidates = candidates[:top_n]

        best_pid = None
        best_dist = None
        best_score = None
        for pid, score in candidates:
            target = seqs[pid]
            try:
                res = edlib.align(seq_u, target, mode="HW", task="distance", k=max_edits)
            except Exception:
                continue
            dist = res.get("editDistance", -1)
            if dist == -1:
                continue
            if best_dist is None or dist < best_dist:
                best_dist = dist
                best_pid = pid
                best_score = score
                if dist == 0:
                    break

        if best_pid is not None:
            matched_records.append((hdr, seq_u, best_pid, best_dist if best_dist is not None else -1, best_score if best_score is not None else 0))

    return matched_records


# -----------------------
# Translation helper
# -----------------------
CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

def translate_dna(seq: str) -> str:
    """Translate DNA (T) to amino acids. Unknown codons -> X. Ignores trailing bases if not divisible by 3."""
    if not seq:
        return ""
    s = seq.upper().replace('U', 'T')
    aa = []
    L = (len(s) // 3) * 3
    for i in range(0, L, 3):
        codon = s[i:i+3]
        aa.append(CODON_TABLE.get(codon, 'X'))
    return "".join(aa)


# -----------------------
# Main routine
# -----------------------
def main():
    p = argparse.ArgumentParser(description="Count exact occurrences of peptides in reads FASTA with optional fuzzy fallback.")
    p.add_argument("--reads", "-r", required=True, help="Trimmed reads FASTA file (can be .gz)")
    p.add_argument("--peptides", "-p", required=True, help="Peptides FASTA file (can be .gz)")
    p.add_argument("--out", "-o", required=True, help="Output CSV path (peptide, count)")
    p.add_argument("--tolerate-errors", action="store_true", help="Enable fuzzy fallback for unmatched reads (k-mer prefilter + edlib)")
    p.add_argument("--kmer", type=int, default=8, help="k-mer size for prefilter (default 8)")
    p.add_argument("--top-candidates", type=int, default=15, help="Number of top k-mer candidates to consider per read (default 15)")
    p.add_argument("--threads", "-t", type=int, default=1, help="Number of worker processes for fuzzy pass (default 1)")
    p.add_argument("--batch-size", type=int, default=2000, help="Number of unmatched reads to send per worker batch (default 2000)")
    p.add_argument("--unmatched-out", help="If set, write reads that failed the exact pass (FASTA) to this path")
    p.add_argument("--fuzzy-troubleshoot", action="store_true", help="Write a CSV with details for each read matched by fuzzy matching (only for reads unmatched by exact pass)")
    p.add_argument("--min-shared-kmers", type=int, default=10, help="Minimum number of shared k-mers between read and peptide for candidate to be considered (default 10)")
    args = p.parse_args()

    if args.tolerate_errors and edlib is None:
        print("ERROR: --tolerate-errors requires edlib. Install with `pip install edlib`", file=sys.stderr)
        sys.exit(1)

    # load peptides
    peptide_names, peptide_seqs, seq_to_indices = load_peptides(args.peptides)
    name_to_seq = {n: s for n, s in zip(peptide_names, peptide_seqs)}

    n_peptides_entries = len(peptide_names)
    pid_to_len = [len(s) for s in peptide_seqs]

    counts = [0] * n_peptides_entries

    tmp = tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".fasta")
    tmp_name = tmp.name
    tmp.close()

    total_reads, matched_reads = exact_pass_and_write_unmatched(args.reads, seq_to_indices, counts, tmp_name)

    fuzzy_matches = 0
    if args.tolerate_errors:
        k = args.kmer
        kmer_index = build_kmer_index(peptide_seqs, k)

        def batch_iterator_reads(path, batch_size):
            it = read_fasta_stream(path)
            batch = []
            for hdr, seq in it:
                batch.append((hdr, seq))
                if len(batch) >= batch_size:
                    yield batch
                    batch = []
            if batch:
                yield batch

        try:
            estimated_unmatched = count_headers(tmp_name)
        except Exception:
            estimated_unmatched = None

        from concurrent.futures import ProcessPoolExecutor, as_completed

        top_n = args.top_candidates
        batch_size = args.batch_size
        threads = max(1, args.threads)

        matched_records = []

        if threads > 1:
            with ProcessPoolExecutor(max_workers=threads) as exe:
                futures = []
                for batch in batch_iterator_reads(tmp_name, batch_size):
                    fut = exe.submit(_fuzzy_worker, batch, peptide_seqs, kmer_index, pid_to_len, top_n, k)
                    futures.append(fut)
                if estimated_unmatched:
                    pbar = tqdm(total=estimated_unmatched, unit="reads", desc="Fuzzy pass")
                else:
                    pbar = tqdm(unit="reads", desc="Fuzzy pass")
                for fut in as_completed(futures):
                    recs = fut.result()
                    matched_records.extend(recs)
                    pbar.update(len(recs))
                pbar.close()
        else:
            if estimated_unmatched:
                pbar = tqdm(total=estimated_unmatched, unit="reads", desc="Fuzzy pass")
            else:
                pbar = tqdm(unit="reads", desc="Fuzzy pass")
            for batch in batch_iterator_reads(tmp_name, batch_size):
                recs = _fuzzy_worker(batch, peptide_seqs, kmer_index, pid_to_len, top_n, k, args.min_shared_kmers)
                matched_records.extend(recs)
                pbar.update(len(batch))
            pbar.close()

        for (_hdr, _seq, pid, _dist, _score) in matched_records:
            counts[pid] += 1

        fuzzy_matches = len(matched_records)

        if args.fuzzy_troubleshoot and matched_records:
            troubleshoot_path = os.path.splitext(args.out)[0] + "_fuzzy_troubleshoot.csv"
            with open(troubleshoot_path, "w", newline="") as tfh:
                writer = csv.writer(tfh)
                writer.writerow(["read_header", "read_seq", "peptide_name", "peptide_seq", "edit_distance", "shared_kmers"])
                for hdr, read_seq, pid, dist, score in matched_records:
                    peptide_name = peptide_names[pid]
                    peptide_seq = peptide_seqs[pid]
                    writer.writerow([hdr, read_seq, peptide_name, peptide_seq, dist, score])
            print(f"Fuzzy troubleshooting CSV written to: {troubleshoot_path}")

    if args.unmatched_out:
        try:
            outdir = os.path.dirname(os.path.abspath(args.unmatched_out))
            if outdir and not os.path.exists(outdir):
                os.makedirs(outdir, exist_ok=True)
            os.replace(tmp_name, args.unmatched_out)
            print(f"Unmatched reads (from exact pass) written to: {args.unmatched_out}")
            tmp_name = None
        except Exception as e:
            print(f"Warning: failed to write unmatched reads to {args.unmatched_out}: {e}", file=sys.stderr)

    try:
        if tmp_name:
            os.unlink(tmp_name)
    except Exception:
        pass

    # -----------------------
    # Build grouped representation (base -> suffix -> count & seq)
    # -----------------------
    groups: Dict[str, Dict[int, int]] = {}
    groups_seq: Dict[str, Dict[int, str]] = {}
    max_suffix = 0
    has_unsuffixed = False

    for name, c, seq in zip(peptide_names, counts, peptide_seqs):
        parts = name.rsplit(".", 1)
        if len(parts) == 2 and parts[1].isdigit():
            base = parts[0]
            suf = int(parts[1])
        else:
            base = name
            suf = 0
            has_unsuffixed = True
        groups.setdefault(base, {})[suf] = c
        groups_seq.setdefault(base, {})[suf] = seq
        if suf > max_suffix:
            max_suffix = suf

    grouped_peptides_count = len(groups)
    identified_grouped = sum(1 for g in groups.values() if any(count > 0 for count in g.values()))
    not_identified_grouped = grouped_peptides_count - identified_grouped

    # Compose header: Peptide, Translation, [Count (unsuffixed) maybe], Count .1..Count .max
    header = ["Peptide", "Translation"]
    if has_unsuffixed:
        header.append("Count")
    for i in range(1, max_suffix + 1):
        header.append(f"Count .{i}")

    # add Sum and Frequency columns to header (do this once, before writing the file)
    header.append("Sum")
    header.append("Frequency")

    # total of all counts (used to compute per-group frequency)
    total_sum = sum(counts)

    # -----------------------
    # Write CSV: one row per base with Translation column before counts
    # -----------------------
    with open(args.out, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        for base in sorted(groups.keys()):
            grp = groups[base]
            seq_map = groups_seq.get(base, {})
            # pick representative seq for translation: prefer .1, else lowest suffix present (incl 0)
            rep_seq = None
            if 1 in seq_map:
                rep_seq = seq_map[1]
            else:
                if seq_map:
                    rep_seq = seq_map[sorted(seq_map.keys())[0]]
            translation = translate_dna(rep_seq) if rep_seq else ""
            # build row: base, translation, then counts
            row = [base, translation]
            if has_unsuffixed:
                # include unsuffixed value or blank if absent
                if 0 in grp:
                    row.append(grp[0])
                else:
                    row.append("")
            for i in range(1, max_suffix + 1):
                row.append(grp.get(i, ""))

            # compute numeric sum across the count columns (skip blanks), append Sum and Frequency
            group_vals = [int(v) for v in row[2:] if v not in ("", None)]
            group_sum = sum(group_vals)
            row.append(group_sum)
            freq = (group_sum / total_sum) if total_sum else 0.0
            row.append(f"{freq:.6f}")
            writer.writerow(row)

    # Summary (grouped)
    unmatched_reads = total_reads - matched_reads
    total_matched_reads = matched_reads + (fuzzy_matches if args.tolerate_errors else 0)
    pct_matched = (total_matched_reads / total_reads * 100) if total_reads else 0.0
    pct_unmatched = (total_reads - total_matched_reads) / total_reads * 100 if total_reads else 0.0

    pct_identified = (identified_grouped / grouped_peptides_count * 100) if grouped_peptides_count else 0.0
    pct_not_identified = (not_identified_grouped / grouped_peptides_count * 100) if grouped_peptides_count else 0.0

    print("\n--- Summary ---")
    print(f"Sample: {basename_without_fasta(args.reads)}")
    print(f"Total reads processed: {total_reads}")
    print(f"Reads matched (exact): {matched_reads}")
    if args.tolerate_errors:
        print(f"Reads matched (fuzzy): {fuzzy_matches}")
    print(f"Reads matched (total): {total_matched_reads} ({pct_matched:.2f}%)")
    print(f"Reads unmatched: {total_reads - total_matched_reads} ({pct_unmatched:.2f}%)")
    print(f"Peptide sequences: {grouped_peptides_count}")
    print(f"Peptides identified (count > 0):   {identified_grouped} ({pct_identified:.2f}%)")
    print(f"Peptides not identified (count=0): {not_identified_grouped} ({pct_not_identified:.2f}%)")
    print("----------------")
    print(f"Counts written to {args.out}")


if __name__ == "__main__":
    main()