#!/usr/bin/env python3
"""
trim.py

Trim reads by 5' and 3' constant sequences (handles reverse complement).

Usage:
    python trim.py -r reads.fasta -c constants.fasta [-o out.fasta] [--troubleshoot] [--examples N]

Behavior:
  - By default the script is tolerant (fuzzy) when searching for constants.
    Use --direct-match to require exact (direct) text matches.
  - Looks for constants on forward strand first; if not found, searches reverse-complement.
  - A read is considered "successfully trimmed" ONLY if both 5' and 3' constants are found on the SAME strand.
  - Partial trims (only one constant found) are still written to the output FASTA (unless you change behavior),
    but are NOT counted as "successfully trimmed".
  - Troubleshoot mode produces:
      * counts and percentages for categories (both found / missing 5' / missing 3' / both missing)
      * histogram CSV of lengths after trimming (counts per length)
      * example FASTA files: passed, missing 5', missing 3', missing both (up to --examples each)
      * percentage of successfully trimmed reads and, among those, percent that are 24, 27, 30, 33, 36 bp long
"""

from tqdm import tqdm
import argparse
import os
import sys
from collections import Counter
from typing import Optional, Tuple, List
import regex

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


def fuzzy_find(seq: str, pattern: str, max_errors: int):
    """
    Return the start index of the first occurrence of pattern in seq allowing up to max_errors
    edits (substitutions/insertions/deletions). Returns None if not found.
    """
    if max_errors <= 0:
        m = seq.find(pattern)
        return m if m != -1 else None
    # regex fuzzy: {e<=N} allows up to N edit operations
    # use overlapped to allow matches starting anywhere
    pat = f'({regex.escape(pattern)})' + f'{{e<={max_errors}}}'
    m = regex.search(pat, seq, regex.BESTMATCH)
    return m.start() if m else None


def find_constants_in_read(read_seq: str, const5: str, const3: str, direct_match: bool = False) -> Tuple[Optional[int], Optional[int], str]:
    """
    Search forward, then reverse-complement using fuzzy_find for tolerant matching by default.
    If direct_match=True, force exact text search (no errors allowed).
    Returns (p5_end, p3_start, strand)
      - p5_end: index immediately after 5' constant, or None
      - p3_start: index of start of 3' constant, or None
      - strand: '+' if forward, '-' if reverse-complement was used
    If both found but in wrong order (3' before end of 5'), treat as not both found (keep partials).
    """
    def _max_errors_for(pattern: str) -> int:
        if not pattern:
            return 0
        if direct_match:
            return 0  # force exact matching when direct-match mode is requested
        # adaptive: 10% of length, at least 1, cap at 3
        return min(3, max(1, int(len(pattern) * 0.1)))

    def search_on(seq):
        # Use fuzzy_find (returns start index or None). If fuzzy_find isn't found it returns None.
        if const5:
            me5 = _max_errors_for(const5)
            p5 = fuzzy_find(seq, const5, me5)
        else:
            p5 = None
        p5_end = (p5 + len(const5)) if p5 is not None else None

        if const3:
            me3 = _max_errors_for(const3)
            p3 = fuzzy_find(seq, const3, me3)
        else:
            p3 = None
        p3_start = p3 if p3 is not None else None

        # If both found but 3' occurs before end of 5', invalidate p3 (order mismatch)
        if p5_end is not None and p3_start is not None and p3_start < p5_end:
            p3_start = None
        return p5_end, p3_start

    # forward
    p5_end, p3_start = search_on(read_seq)
    if p5_end is not None or p3_start is not None:
        return p5_end, p3_start, "+"

    # reverse complement
    rc = revcomp(read_seq)
    p5_end_rc, p3_start_rc = search_on(rc)
    if p5_end_rc is not None or p3_start_rc is not None:
        return p5_end_rc, p3_start_rc, "-"

    return None, None, "+"


def save_examples(prefix: str, name: str, bucket: List[Tuple[str, str, str]], examples_n: int) -> Optional[str]:
    if not bucket:
        return None
    path = f"{prefix}_{name}.fasta"
    with open(path, "w") as fh:
        for hdr, seq, strand in bucket[:examples_n]:
            write_fasta_record(fh, f"{hdr} strand={strand}", seq)
    return path


def write_histogram_csv(path: str, length_counts: Counter):
    """
    Write a CSV with two columns: length,count
    Sorted by length ascending.
    """
    with open(path, "w") as fh:
        fh.write("length,count\n")
        for length in sorted(length_counts.keys()):
            fh.write(f"{length},{length_counts[length]}\n")


def main():
    p = argparse.ArgumentParser(description="Trim reads by 5' and 3' constants (handles reverse complement).")
    p.add_argument("-r", "--reads", required=True, help="Input reads FASTA file")
    p.add_argument("-c", "--constants", required=True, help="Constants FASTA with headers '5prime' and '3prime'")
    p.add_argument("-o", "--out", default=None, help="Output FASTA path (default: add _trimmed to reads filename)")
    p.add_argument("--troubleshoot", action="store_true", help="Produce troubleshooting report and example reads and histogram CSV")
    p.add_argument("--direct-match", action="store_true", help="Require exact direct matches for constants (disable fuzzy/tolerant matching). By default fuzzy/tolerant matching is used.")
    p.add_argument("--examples", type=int, default=10, help="How many example reads to save per category (default 10)")
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
    passed_count = 0  # both constants found on same strand -> successful trim
    length_after = []  # lengths for any trimmed sequences (including partial trims)
    missing5_count = 0  # reads where 5' is missing but 3' found
    missing3_count = 0  # reads where 3' is missing but 5' found
    missing_both_count = 0

    examples_passed = []
    examples_missing_5 = []
    examples_missing_3 = []
    examples_missing_both = []

    # Estimate total number of reads (number of FASTA headers) for an accurate progress bar.
    # This is a fast header-only scan and is much cheaper than full parsing; it requires briefly
    # reading through the file once to count '>' lines.
    try:
        with open(reads_path, "r") as _fh:
            estimated_total = sum(1 for line in _fh if line.startswith(">"))
    except Exception:
        estimated_total = None

    # We'll write trimmed sequences (including partial trims) as before
    with open(out_path, "w") as out_fh:
        iterator = read_fasta(reads_path)
        if estimated_total:
            iterator = tqdm(iterator, total=estimated_total, unit="reads", desc="Processing reads")
        else:
            iterator = tqdm(iterator, unit="reads", desc="Processing reads")
        
        for header, seq in iterator:
            total_reads += 1
            seq_up = seq.upper()
            p5_end, p3_start, strand = find_constants_in_read(seq_up, const5, const3, direct_match=args.direct_match)
            trimmed_seq = None

            # Count categories and gather examples according to presence/absence
            if p5_end is not None and p3_start is not None:
                # both found => successful trim
                if strand == "+":
                    trimmed_seq = seq_up[p5_end:p3_start]
                else:
                    rc = revcomp(seq_up)
                    trimmed_seq = rc[p5_end:p3_start]
                passed_count += 1
                length_after.append(len(trimmed_seq))
                examples_passed.append((header, trimmed_seq, strand))
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
            else:
                # neither found
                missing_both_count += 1
                examples_missing_both.append((header, seq_up, "+"))

            # Write trimmed sequence if any (partial trims are not written)
            if p5_end is not None and p3_start is not None:
                write_fasta_record(out_fh, header, trimmed_seq)

    # Troubleshooting and reporting
    # Always report % successfully trimmed and among those, % of sequences with lengths 24,27,30,33,36
    print("\n--- Summary ---")
    print(f"Total reads processed: {total_reads}")
    pct_passed = (passed_count / total_reads * 100) if total_reads else 0.0
    print(f"Reads with both constants found (successful trims): {passed_count} ({pct_passed:.2f}%)")

    # Among successfully trimmed, percentage that are exactly 24,27,30,33,36 bp
    target_lengths = [24, 27, 30, 33, 36]
    passed_lengths = []
    # To compute these we need the trimmed lengths of successful reads; they are in examples_passed or length_after but length_after includes partials.
    # We'll reconstruct passed lengths from examples_passed
    passed_len_counts = Counter()
    for (_hdr, seq, _strand) in examples_passed:
        passed_len_counts[len(seq)] += 1
        passed_lengths.append(len(seq))

    if passed_count > 0:
        for L in target_lengths:
            cnt = passed_len_counts.get(L, 0)
            pct = cnt / passed_count * 100
            print(f"  Of successfully trimmed reads: {cnt} ({pct:.2f}%) are {L} bp")
        other_count = passed_count - sum(passed_len_counts.get(L, 0) for L in target_lengths)
        print(f"  Of successfully trimmed reads: {other_count} ({other_count / passed_count * 100:.2f}%) are other lengths")
    else:
        for L in target_lengths:
            print(f"  Of successfully trimmed reads: 0 (0.00%) are {L} bp")
        other_count = passed_count - sum(passed_len_counts.get(L, 0) for L in target_lengths)
        print(f"  Of successfully trimmed reads: {other_count} ({other_count / passed_count * 100:.2f}%) are other lengths")

    # Now report the three other categories with counts and percentages (user asked explicitly to report percentages here as well)
    def pct_fmt(count):
        return f"{count} ({(count/total_reads*100 if total_reads else 0.0):.2f}%)"

    print(f"Reads with 5' missing (but 3' found): {pct_fmt(missing5_count)}")
    print(f"Reads with 3' missing (but 5' found): {pct_fmt(missing3_count)}")
    print(f"Reads with both missing: {pct_fmt(missing_both_count)}")

    # Histogram CSV of lengths after trimming (include partial as before)
    if args.troubleshoot:
        lengths_counter = Counter(length_after)
        hist_csv_path = os.path.splitext(out_path)[0] + "_trimmed_length_hist.csv"
        write_histogram_csv(hist_csv_path, lengths_counter)
        print(f"Histogram (length -> count) written to: {hist_csv_path}")

        # Save example FASTA files (up to examples_n each)
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


if __name__ == "__main__":
    main()