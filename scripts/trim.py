#!/usr/bin/env python3
"""
trim.py (updated)

Behavior:
  - Two-phase processing:
      1) Direct (exact) pass: try to find both constants exactly (forward, then RC).
         This pass is fast and reported (time + reads/sec + how many failed exact).
      2) Fuzzy pass: only for reads that failed the exact pass (unless --direct-only).
  - --direct-only: skip fuzzy pass entirely (only exact matches used).
  - --fuzzy-troubleshoot: save a CSV for reads that were fuzzy-searched in ./2-trimmed/(library)_fuzzy_troubleshoot.csv.
  - A summary Excel file with counts+percentages is written to ./2-trimmed/(library)_trim_summary.xlsx.
"""
from tqdm import tqdm
import argparse
import os
import sys
from collections import Counter
from typing import Optional, Tuple, List, Dict, Any
import regex
import csv
import time
from openpyxl import Workbook

COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def read_fasta(fp):
    """Yield (header, seq) from FASTA file path or file-like object."""
    close = False
    if isinstance(fp, str):
        fh = open(fp, "r")
        close = True
    else:
        fh = fp

    header = None
    seq_parts = []
    for line in fh:
        line = line.rstrip("\n\r")
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq_parts)
            header = line[1:].split()[0]
            seq_parts = []
        else:
            seq_parts.append(line.strip())
    if header is not None:
        yield header, "".join(seq_parts)
    if close:
        fh.close()


def write_fasta_record(fh, header: str, seq: str):
    fh.write(f">{header}\n")
    for i in range(0, len(seq), 80):
        fh.write(seq[i : i + 80] + "\n")


def parse_constants_fasta(constants_path: str) -> Tuple[str, str]:
    consts = {}
    for header, seq in read_fasta(constants_path):
        consts[header] = seq.upper()
    if "5prime" not in consts or "3prime" not in consts:
        raise ValueError("Constants FASTA must contain headers exactly '5prime' and '3prime'.")
    return consts["5prime"], consts["3prime"]


def make_default_outpath(reads_path: str) -> str:
    base, ext = os.path.splitext(reads_path)
    if ext == "":
        return reads_path + "_trimmed.fasta"
    else:
        return f"{base}_trimmed{ext}"


def shared_kmers(a: str, b: str, k: int) -> int:
    if k <= 0:
        return len(set(a) & set(b))
    if len(a) < k or len(b) < k:
        return len(set(a) & set(b))
    sa = {a[i : i + k] for i in range(len(a) - k + 1)}
    sb = {b[i : i + k] for i in range(len(b) - k + 1)}
    return len(sa & sb)


def fuzzy_find_with_stats(seq: str, pattern: str, max_errors: int):
    if not pattern:
        return None, None, None, None
    if max_errors <= 0:
        idx = seq.find(pattern)
        if idx == -1:
            return None, None, None, None
        matched = seq[idx : idx + len(pattern)]
        return idx, 0, matched, None
    pat = f'({regex.escape(pattern)})' + f'{{e<={max_errors}}}'
    m = regex.search(pat, seq, regex.BESTMATCH)
    if not m:
        return None, None, None, None
    start = m.start()
    matched_region = m.group(0)
    try:
        ins, dele, sub = m.fuzzy_counts
        edit_distance = ins + dele + sub
    except Exception:
        edit_distance = _levenshtein(pattern, matched_region)
    return start, edit_distance, matched_region, m


def _levenshtein(a: str, b: str) -> int:
    if a == b:
        return 0
    la, lb = len(a), len(b)
    if la == 0:
        return lb
    if lb == 0:
        return la
    prev = list(range(lb + 1))
    for i, ca in enumerate(a, 1):
        cur = [i] + [0] * lb
        for j, cb in enumerate(b, 1):
            cost = 0 if ca == cb else 1
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + cost)
        prev = cur
    return prev[lb]


def _max_errors_for(pattern: str) -> int:
    if not pattern:
        return 0
    # adaptive: 10% of length, at least 1, cap at 3
    return min(3, max(1, int(len(pattern) * 0.1)))


def exact_search_on(seq: str, const5: str, const3: str) -> Tuple[Optional[int], Optional[int]]:
    """Return p5_end, p3_start for exact matches on seq. If order mismatch, invalidate p3."""
    p5 = seq.find(const5) if const5 else -1
    p5_idx = p5 if p5 != -1 else None
    p5_end = (p5_idx + len(const5)) if p5_idx is not None else None

    p3 = seq.find(const3) if const3 else -1
    p3_start = p3 if p3 != -1 else None

    if p5_end is not None and p3_start is not None and p3_start < p5_end:
        p3_start = None
    return p5_end, p3_start


def find_constants_fuzzy_on(seq: str, const5: str, const3: str):
    """Fuzzy search both constants on a single sequence (forward or RC)."""
    me5 = _max_errors_for(const5)
    me3 = _max_errors_for(const3)
    p5, ed5, m5, _ = fuzzy_find_with_stats(seq, const5, me5)
    p3, ed3, m3, _ = fuzzy_find_with_stats(seq, const3, me3)
    p5_end = (p5 + len(m5)) if p5 is not None else None
    p3_start = p3 if p3 is not None else None
    if p5_end is not None and p3_start is not None and p3_start < p5_end:
        p3_start = None
    return p5_end, p3_start, ed5, ed3, m5, m3


def save_examples(prefix: str, name: str, bucket: List[Tuple[str, str, str]], examples_n: int) -> Optional[str]:
    if not bucket:
        return None
    path = f"{prefix}_{name}.fasta"
    with open(path, "w") as fh:
        for hdr, seq, strand in bucket[:examples_n]:
            write_fasta_record(fh, f"{hdr} strand={strand}", seq)
    return path


def write_histogram_csv(path: str, length_counts: Counter):
    with open(path, "w") as fh:
        fh.write("length,count\n")
        for length in sorted(length_counts.keys()):
            fh.write(f"{length},{length_counts[length]}\n")


def write_summary_excel(path: str, summary_rows: List[Tuple[str, int, float]]):
    """
    summary_rows: list of (label, count, percent)
    We'll write two columns: Count and Percent next to label.
    """
    wb = Workbook()
    ws = wb.active
    ws.title = "trim_summary"
    ws.append(["metric", "count", "percent"])
    for label, count, pct in summary_rows:
        ws.append([label, count, f"{pct:.2f}%"])
    wb.save(path)


def main():
    p = argparse.ArgumentParser(description="Trim reads by 5' and 3' constants (handles reverse complement).")
    p.add_argument("-r", "--reads", required=True, help="Input reads FASTA file")
    p.add_argument("-c", "--constants", required=True, help="Constants FASTA with headers '5prime' and '3prime'")
    p.add_argument("-o", "--out", default=None, help="Output FASTA path (default: add _trimmed to reads filename)")
    p.add_argument("--troubleshoot", default=True, action="store_true", help="Produce troubleshooting report and example reads and histogram CSV")
    p.add_argument("--examples", type=int, default=10, help="How many example reads to save per category (default 10)")
    p.add_argument("--direct-only", action="store_true", help="Only use exact matches; do not run fuzzy searches on reads that fail exact matching.")
    p.add_argument("--fuzzy-troubleshoot", action="store_true", help="Save a fuzzy-troubleshoot CSV for reads that were fuzzy-searched.")
    args = p.parse_args()

    reads_path = args.reads
    const_path = args.constants
    out_path = args.out or make_default_outpath(reads_path)
    examples_n = args.examples

    try:
        const5, const3 = parse_constants_fasta(const_path)
    except Exception as e:
        print(f"Error reading constants file: {e}", file=sys.stderr)
        sys.exit(1)

    const5 = const5.upper()
    const3 = const3.upper()

    total_reads = 0

    # Counters
    passed_count = 0
    length_after = []
    missing5_count = 0
    missing3_count = 0
    missing_both_count = 0

    examples_passed = []
    examples_missing_5 = []
    examples_missing_3 = []
    examples_missing_both = []

    # Paths & dirs for fuzzy outputs and excel summary
    lib_name = os.path.splitext(os.path.basename(reads_path))[0]
    out_dir = "./2-trimmed"
    os.makedirs(out_dir, exist_ok=True)
    fuzzy_csv_path = os.path.join(out_dir, f"{lib_name}_fuzzy_troubleshoot.csv")
    excel_path = os.path.join(out_dir, f"{lib_name}_trim_summary.xlsx")

    # Prepare fuzzy CSV writer if requested
    fuzzy_csv_fh = None
    fuzzy_csv_writer = None
    if args.fuzzy_troubleshoot:
        fuzzy_csv_fh = open(fuzzy_csv_path, "w", newline="")
        fuzzy_csv_writer = csv.writer(fuzzy_csv_fh)
        fuzzy_csv_writer.writerow([
            "header",
            "full_sequence",
            "trimmed_sequence",
            "edit_distance_5",
            "shared_kmers_5",
            "edit_distance_3",
            "shared_kmers_3",
            "strand",
            "status",
            "method"
        ])

    # Estimate total reads for progress bars
    try:
        with open(reads_path, "r") as _fh:
            estimated_total = sum(1 for line in _fh if line.startswith(">"))
    except Exception:
        estimated_total = None

    # PHASE 1: Direct (exact) pass
    direct_failed_headers = set()
    direct_passed_count = 0
    t0 = time.perf_counter()
    # We'll open output FASTA early and write exact-successful trims immediately.
    with open(out_path, "w") as out_fh:
        iterator = read_fasta(reads_path)
        if estimated_total:
            iterator = tqdm(iterator, total=estimated_total, unit="reads", desc="Direct pass (exact search)")
        else:
            iterator = tqdm(iterator, unit="reads", desc="Direct pass (exact search)")

        for header, seq in iterator:
            total_reads += 1
            seq_up = seq.upper()
            p5_end, p3_start = exact_search_on(seq_up, const5, const3)
            # If exact both found on forward, accept
            if p5_end is not None and p3_start is not None:
                trimmed_seq = seq_up[p5_end:p3_start]
                direct_passed_count += 1
                passed_count += 1
                length_after.append(len(trimmed_seq))
                examples_passed.append((header, trimmed_seq, "+"))
                write_fasta_record(out_fh, header, trimmed_seq)
                continue
            # Else check exact on RC
            rc = revcomp(seq_up)
            p5_end_rc, p3_start_rc = exact_search_on(rc, const5, const3)
            if p5_end_rc is not None and p3_start_rc is not None:
                trimmed_seq = rc[p5_end_rc:p3_start_rc]
                direct_passed_count += 1
                passed_count += 1
                length_after.append(len(trimmed_seq))
                examples_passed.append((header, trimmed_seq, "-"))
                write_fasta_record(out_fh, header, trimmed_seq)
                continue
            # Not both exact on either strand -> add to set for fuzzy pass (unless direct-only)
            direct_failed_headers.add(header)
    t1 = time.perf_counter()
    direct_time = t1 - t0
    direct_failed_count = len(direct_failed_headers)
    direct_speed = (total_reads / direct_time) if direct_time > 0 else float("inf")
    print(f"Direct search: processed {total_reads} reads in {direct_time:.2f}s ({direct_speed:.1f} reads/s)")
    print(f"Direct search: {direct_passed_count} reads had both exact constants; {direct_failed_count} reads require fuzzy (if enabled).")

    # If direct-only, skip fuzzy pass entirely; we are done. (Note: trimmed FASTA already contains direct-successful trims)
    fuzzy_processed = 0
    fuzzy_passed = 0
    fuzzy_time = 0.0
    if not args.direct_only and direct_failed_count > 0:
        # PHASE 2: Fuzzy pass — iterate reads again, but only process ones in direct_failed_headers
        t2 = time.perf_counter()
        with open(out_path, "a") as out_fh:  # append fuzzy successful trims to the same output FASTA
            iterator = read_fasta(reads_path)
            if estimated_total:
                iterator = tqdm(iterator, total=estimated_total, unit="reads", desc="Fuzzy pass (for failed exacts)")
            else:
                iterator = tqdm(iterator, unit="reads", desc="Fuzzy pass (for failed exacts)")
            for header, seq in iterator:
                if header not in direct_failed_headers:
                    continue
                seq_up = seq.upper()
                fuzzy_processed += 1
                # forward fuzzy
                p5_end, p3_start, ed5, ed3, m5, m3 = find_constants_fuzzy_on(seq_up, const5, const3)
                method = "fuzzy"
                strand = "+"
                if p5_end is None and p3_start is None:
                    # try RC
                    rc = revcomp(seq_up)
                    p5_end_rc, p3_start_rc, ed5_rc, ed3_rc, m5_rc, m3_rc = find_constants_fuzzy_on(rc, const5, const3)
                    if p5_end_rc is not None or p3_start_rc is not None:
                        p5_end, p3_start = p5_end_rc, p3_start_rc
                        ed5, ed3, m5, m3 = ed5_rc, ed3_rc, m5_rc, m3_rc
                        strand = "-"
                # classify and write if both found
                if p5_end is not None and p3_start is not None:
                    fuzzy_passed += 1
                    passed_count += 1
                    if strand == "+":
                        trimmed_seq = seq_up[p5_end:p3_start]
                    else:
                        rc = revcomp(seq_up)
                        trimmed_seq = rc[p5_end:p3_start]
                    length_after.append(len(trimmed_seq))
                    examples_passed.append((header, trimmed_seq, strand))
                    write_fasta_record(out_fh, header, trimmed_seq)
                    status = "both_found"
                elif p5_end is None and p3_start is not None:
                    # 5' missing, 3' found
                    if strand == "+":
                        trimmed_seq = seq_up[:p3_start]
                    else:
                        rc = revcomp(seq_up)
                        trimmed_seq = rc[:p3_start]
                    missing5_count += 1
                    length_after.append(len(trimmed_seq))
                    examples_missing_5.append((header, trimmed_seq, strand))
                    status = "missing5"
                elif p5_end is not None and p3_start is None:
                    # 3' missing, 5' found
                    if strand == "+":
                        trimmed_seq = seq_up[p5_end:]
                    else:
                        rc = revcomp(seq_up)
                        trimmed_seq = rc[p5_end:]
                    missing3_count += 1
                    length_after.append(len(trimmed_seq))
                    examples_missing_3.append((header, trimmed_seq, strand))
                    status = "missing3"
                else:
                    # neither found
                    missing_both_count += 1
                    examples_missing_both.append((header, seq_up, "+"))
                    trimmed_seq = ""
                    status = "missing_both"

                # optionally log fuzzy-troubleshoot rows (only for reads that were fuzzied)
                if args.fuzzy_troubleshoot and fuzzy_csv_writer is not None:
                    sk5 = shared_kmers(const5, m5 if m5 else "", min(4, len(const5)))
                    sk3 = shared_kmers(const3, m3 if m3 else "", min(4, len(const3)))
                    fuzzy_csv_writer.writerow([
                        header,
                        seq_up,
                        trimmed_seq if trimmed_seq is not None else "",
                        ed5 if ed5 is not None else "",
                        sk5 if sk5 is not None else 0,
                        ed3 if ed3 is not None else "",
                        sk3 if sk3 is not None else 0,
                        strand,
                        status,
                        method
                    ])
        t3 = time.perf_counter()
        fuzzy_time = t3 - t2
        fuzzy_speed = (fuzzy_processed / fuzzy_time) if fuzzy_time > 0 else float("inf")
        print(f"Fuzzy search: processed {fuzzy_processed} reads in {fuzzy_time:.2f}s ({fuzzy_speed:.1f} reads/s)")
        print(f"Fuzzy search: {fuzzy_passed} reads had both constants after fuzzy search.")

    # Close fuzzy CSV if open
    if fuzzy_csv_fh:
        fuzzy_csv_fh.close()
        print(f"Fuzzy troubleshoot CSV written to: {fuzzy_csv_path}")

    # FINAL COUNTS: note direct_passed_count + fuzzy_passed == passed_count (should hold)
    # Prepare summary report (counts and percentages)
    print("\n--- Summary ---")
    print(f"Total reads processed: {total_reads}")
    pct_passed = (passed_count / total_reads * 100) if total_reads else 0.0
    print(f"Reads with both constants found (successful trims): {passed_count} ({pct_passed:.2f}%)")

    # compute target length breakdown among successfully trimmed reads
    target_lengths = [24, 27, 30, 33, 36]
    passed_len_counts = Counter()
    # We have examples_passed but it may not contain all passed reads if large — however we recorded lengths in length_after for every successful trim
    # We'll build counts from length_after but length_after includes partial trims lengths (we appended only successful trims to length_after for passed_count)
    # In this implementation, length_after contains lengths for successful trims and also for partial trims encountered in fuzzy pass;
    # to be safe, we'll reconstruct passed length counts by re-counting examples_passed which we appended on every successful trim.
    passed_len_counts = Counter(len(seq) for (_hdr, seq, _s) in examples_passed)

    # print target breakdown
    for L in target_lengths:
        cnt = passed_len_counts.get(L, 0)
        pct = (cnt / passed_count * 100) if passed_count else 0.0
        print(f"  Of successfully trimmed reads: {cnt} ({pct:.2f}%) are {L} bp")
    other_count = passed_count - sum(passed_len_counts.get(L, 0) for L in target_lengths)
    other_pct = (other_count / passed_count * 100) if passed_count else 0.0
    print(f"  Of successfully trimmed reads: {other_count} ({other_pct:.2f}%) are other lengths")

    def pct_fmt(count):
        return f"{count} ({(count/total_reads*100 if total_reads else 0.0):.2f}%)"

    print(f"Reads with 5' missing (but 3' found): {pct_fmt(missing5_count)}")
    print(f"Reads with 3' missing (but 5' found): {pct_fmt(missing3_count)}")
    print(f"Reads with both missing: {pct_fmt(missing_both_count)}")

    # Histogram CSV (lengths including successful trims and partials recorded)
    if args.troubleshoot:
        lengths_counter = Counter(length_after)
        hist_csv_path = os.path.splitext(out_path)[0] + "_trimmed_length_hist.csv"
        write_histogram_csv(hist_csv_path, lengths_counter)
        print(f"Histogram (length -> count) written to: {hist_csv_path}")

        prefix = os.path.splitext(out_path)[0]
        p_pass = save_examples(prefix, "examples_passed", examples_passed, examples_n)
        p_m5 = save_examples(prefix, "examples_missing_5prime", examples_missing_5, examples_n)
        p_m3 = save_examples(prefix, "examples_missing_3prime", examples_missing_3, examples_n)
        p_mb = save_examples(prefix, "examples_missing_both", examples_missing_both, examples_n)

        if p_pass:
            print(f"Examples (passed) saved to {p_pass}")
        if p_m5:
            print(f"Examples (missing 5') saved to {p_m5}")
        if p_m3:
            print(f"Examples (missing 3') saved to {p_m3}")
        if p_mb:
            print(f"Examples (missing both) saved to {p_mb}")

    print(f"\nTrimmed reads written to: {out_path}")

    # Write Excel summary with counts and percentages (shape similar to your example)
    summary_rows = []
    summary_rows.append(("Total reads processed", total_reads, 100.0))
    summary_rows.append(("Reads with both constants found (successful trims)", passed_count, pct_passed))
    for L in target_lengths:
        cnt = passed_len_counts.get(L, 0)
        pct = (cnt / passed_count * 100) if passed_count else 0.0
        summary_rows.append((f"Of successfully trimmed reads: {L} bp", cnt, pct))
    summary_rows.append(("Of successfully trimmed reads: other lengths", other_count, other_pct))
    summary_rows.append(("Reads with 5' missing (but 3' found)", missing5_count, (missing5_count/total_reads*100 if total_reads else 0.0)))
    summary_rows.append(("Reads with 3' missing (but 5' found)", missing3_count, (missing3_count/total_reads*100 if total_reads else 0.0)))
    summary_rows.append(("Reads with both missing", missing_both_count, (missing_both_count/total_reads*100 if total_reads else 0.0)))

    try:
        write_summary_excel(excel_path, summary_rows)
        print(f"Trim summary Excel written to: {excel_path}")
    except Exception as e:
        print(f"Failed to write Excel summary: {e}", file=sys.stderr)

if __name__ == "__main__":
    main()