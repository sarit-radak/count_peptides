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
    No progress bar here per your earlier request.
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
            # treat the whole seq as a key if shorter than k
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
                # wrap seq lines at 80 chars
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


def candidates_from_kmers(seq: str, k: int, kmer_index: Dict[str, List[int]], length_to_ids: Dict[int, List[int]], max_edit: int, allow_len_slop: int = 1) -> List[Tuple[int, int]]:
    """
    Given a read seq, return candidate peptide indices with a score (shared-kmer-count).
    We also use a loose length filter: peptides with length within [len-max_edit-allow_len_slop, len+max_edit+allow_len_slop]
    Returns list of tuples (peptide_idx, score) sorted descending by score.
    """
    L = len(seq)
    min_len = max(1, L - max_edit - allow_len_slop)
    max_len = L + max_edit + allow_len_slop

    # collect candidate ids via k-mer hits
    counter = Counter()
    if len(seq) < k:
        # if read shorter than k, match exact kmer as the whole seq
        ids = kmer_index.get(seq, [])
        for pid in ids:
            counter[pid] += 1
    else:
        seen_kmers = set()
        for i in range(len(seq) - k + 1):
            kmer = seq[i : i + k]
            if kmer in seen_kmers:
                continue
            seen_kmers.add(kmer)
            ids = kmer_index.get(kmer)
            if ids:
                for pid in ids:
                    counter[pid] += 1

    # filter by length buckets
    results = []
    for pid, score in counter.items():
        # we can't access peptide lengths here; length_to_ids is len->list(ids).
        # To filter efficiently, we check whether pid is in an allowed length bucket using a small map:
        # We'll rely on a global pid_to_len map provided externally (see caller).
        results.append((pid, score))
    # sort by score descending
    results.sort(key=lambda x: -x[1])
    return results


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

        # Apply minimum shared-kmer filter to reduce candidates aggressively
        if not counter:
            continue
        filtered = [(pid, score) for pid, score in counter.items() if score >= min_shared_kmers]
        if not filtered:
            continue

        # Filter by length tolerance (tighter filter)
        candidates = []
        for pid, score in filtered:
            plen = pid_to_len[pid]
            if abs(plen - L) <= max_edits + 0:
                candidates.append((pid, score))

        if not candidates:
            continue

        # pick top_n candidates by score (keep score attached)
        candidates.sort(key=lambda x: -x[1])
        candidates = candidates[:top_n]  # list of (pid, score)

        # edlib alignment over candidates, keep best (lowest edit distance)
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
            # return header and sequence so the caller can write troubleshooting rows
            matched_records.append((hdr, seq_u, best_pid, best_dist if best_dist is not None else -1, best_score if best_score is not None else 0))

    return matched_records

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

    # load peptides (no loading progress bar)
    peptide_names, peptide_seqs, seq_to_indices = load_peptides(args.peptides)
    n_peptides = len(peptide_names)
    pid_to_len = [len(s) for s in peptide_seqs]

    # counts initialized to zero
    counts = [0] * n_peptides

    # temp file for unmatched reads
    tmp = tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".fasta")
    tmp_name = tmp.name
    tmp.close()

    # exact pass: write unmatched to tmp file and update counts
    total_reads, matched_reads = exact_pass_and_write_unmatched(args.reads, seq_to_indices, counts, tmp_name)

    # If fuzzy fallback is not enabled, we're done
    fuzzy_matches = 0
    if args.tolerate_errors:
        # Build k-mer index of peptides
        k = args.kmer
        kmer_index = build_kmer_index(peptide_seqs, k)

        # Process unmatched reads in batches, optionally with multiprocessing
        # We'll read tmp_name and create batches of (hdr, seq)
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

        # estimate number of unmatched reads for progress bar
        try:
            estimated_unmatched = count_headers(tmp_name)
        except Exception:
            estimated_unmatched = None

        # If threads > 1, use ProcessPoolExecutor
        from concurrent.futures import ProcessPoolExecutor, as_completed

        top_n = args.top_candidates
        batch_size = args.batch_size
        threads = max(1, args.threads)

        # We'll accumulate matched peptide indices returned by workers and update counts afterwards
        #matched_pid_list = []
        matched_records = []  # list of tuples (hdr, read_seq, best_pid, best_dist, best_score)

        if threads > 1:
            with ProcessPoolExecutor(max_workers=threads) as exe:
                futures = []
                for batch in batch_iterator_reads(tmp_name, batch_size):
                    # submit batch to worker
                    fut = exe.submit(_fuzzy_worker, batch, peptide_seqs, kmer_index, pid_to_len, top_n, k)
                    futures.append(fut)
                # progress bar while futures complete
                if estimated_unmatched:
                    pbar = tqdm(total=estimated_unmatched, unit="reads", desc="Fuzzy pass")
                else:
                    pbar = tqdm(unit="reads", desc="Fuzzy pass")
                for fut in as_completed(futures):
                        recs = fut.result()  # list of (hdr, seq, pid, dist, score)
                        matched_records.extend(recs)
                        # For progress bar: we don't know how many reads were processed in this future,
                        # so update conservatively by number of records returned (this keeps the bar moving).
                        pbar.update(len(recs)) # advance by number of processed reads isn't directly known, keep bar moving per batch size:
                    # we don't know the exact number of processed reads here; approximate by batch_size steps:
                    # but safer: pbar.update(batch_size) might overshoot at end; instead compute processed by summing counts:
                    # To keep it simple and avoid incorrect ETA, we'll increment by len(matched_pids) as a conservative progress update.
                pbar.close()
        else:
            # single-threaded: simply iterate batches and call worker directly
            if estimated_unmatched:
                pbar = tqdm(total=estimated_unmatched, unit="reads", desc="Fuzzy pass")
            else:
                pbar = tqdm(unit="reads", desc="Fuzzy pass")
            for batch in batch_iterator_reads(tmp_name, batch_size):
                recs = _fuzzy_worker(batch, peptide_seqs, kmer_index, pid_to_len, top_n, k, args.min_shared_kmers)
                matched_records.extend(recs)
                pbar.update(len(batch))
            pbar.close()

        # Update counts from matched_pid_list
# Update counts from matched_records
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
            # Ensure destination directory exists
            outdir = os.path.dirname(os.path.abspath(args.unmatched_out))
            if outdir and not os.path.exists(outdir):
                os.makedirs(outdir, exist_ok=True)
            # atomically move/replace
            os.replace(tmp_name, args.unmatched_out)
            print(f"Unmatched reads (from exact pass) written to: {args.unmatched_out}")
            # Set tmp_name to None so we don't attempt to delete it below
            tmp_name = None

        except Exception as e:
            print(f"Warning: failed to write unmatched reads to {args.unmatched_out}: {e}", file=sys.stderr)

    # Cleanup temp file
    try:
        os.unlink(tmp_name)
    except Exception:
        pass

    # Write transposed CSV: peptide, count
    with open(args.out, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["Peptide", "Count"])
        for name, c in zip(peptide_names, counts):
            writer.writerow([name, c])

    # Summary
    unmatched_reads = total_reads - matched_reads
    # include fuzzy matches into matched_reads for summary if used
    total_matched_reads = matched_reads + (fuzzy_matches if args.tolerate_errors else 0)
    pct_matched = (total_matched_reads / total_reads * 100) if total_reads else 0.0
    pct_unmatched = (total_reads - total_matched_reads) / total_reads * 100 if total_reads else 0.0

    identified_peptides = sum(1 for c in counts if c > 0)
    not_identified_peptides = n_peptides - identified_peptides
    pct_identified = (identified_peptides / n_peptides * 100) if n_peptides else 0.0
    pct_not_identified = (not_identified_peptides / n_peptides * 100) if n_peptides else 0.0

    print("\n--- Summary ---")
    print(f"Sample: {basename_without_fasta(args.reads)}")
    print(f"Total reads processed: {total_reads}")
    print(f"Reads matched (exact): {matched_reads}")
    if args.tolerate_errors:
        print(f"Reads matched (fuzzy): {fuzzy_matches}")
    print(f"Reads matched (total): {total_matched_reads} ({pct_matched:.2f}%)")
    print(f"Reads unmatched: {total_reads - total_matched_reads} ({pct_unmatched:.2f}%)")
    print(f"Peptide sequences: {n_peptides}")
    print(f"Peptides identified (count > 0):   {identified_peptides} ({pct_identified:.2f}%)")
    print(f"Peptides not identified (count=0): {not_identified_peptides} ({pct_not_identified:.2f}%)")
    print("----------------")
    print(f"Counts written to {args.out}")


if __name__ == "__main__":
    main()