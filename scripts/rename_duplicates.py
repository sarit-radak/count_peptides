#!/usr/bin/env python3
"""
rename_duplicates.py

Usage:
    python3 -u ./scripts/rename_duplicates.py -i INPUT.fasta -o OUTPUT.fasta [--troubleshoot]

Behavior:
  - Translates each nucleotide sequence (standard genetic code; truncates remainder bases).
  - Groups sequences by translated peptide.
  - For each peptide group, builds a base name by joining unique trimmed headers (header before the first '.')
    with commas: e.g. KnownMHCEpep_12mer_0001,KnownMHCEpep_12mer_0002
  - Renames each sequence in that group to <base>.<N> where N is a sequential codon variant index (1-based).
  - Writes renamed nucleotide sequences to the output FASTA.
  - Prints a histogram of how many peptides occurred N times.
  - If --troubleshoot is provided, creates an XLSX file (named like the output with .xlsx appended)
    where each sheet is sequences that appeared K times (sheet name "count_K") and contains:
      new_name | original_header | trimmed_name | nucleotide_seq | aa_seq
"""

from __future__ import annotations
import argparse
import sys
from collections import defaultdict, Counter, OrderedDict
from typing import Dict, List, Tuple
import textwrap

# Standard codon table (DNA codons; stop codons mapped to '*')
CODON_TABLE = {
    # U->T replaced below manually by writing DNA codons
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

def parse_fasta(path: str) -> List[Tuple[str,str]]:
    """Return list of (header, seq) preserving input order."""
    records = []
    header = None
    seq_lines = []
    with open(path, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(seq_lines).replace(' ', '').replace('\r','').upper()))
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header is not None:
            records.append((header, ''.join(seq_lines).replace(' ', '').replace('\r','').upper()))
    return records

def translate_dna(seq: str) -> str:
    """Translate DNA sequence (A/T/G/C). Trim trailing bases if not full codon."""
    seq = seq.upper().replace('U', 'T')
    aa = []
    limit = (len(seq) // 3) * 3
    for i in range(0, limit, 3):
        codon = seq[i:i+3]
        aa.append(CODON_TABLE.get(codon, 'X'))  # unknown codons -> X
    return ''.join(aa)

def wrap_seq(seq: str, width: int = 60) -> str:
    return '\n'.join(textwrap.wrap(seq, width))

def build_new_names_and_records(records: List[Tuple[str,str]]):
    """
    Process records:
      - trimmed_name = header.split('.')[0]
      - group by translated AA
      - for each AA group, build base name = comma-joined unique trimmed_names (sorted)
      - produce new headers sequentially base.1, base.2, ...
    Returns:
      list of (new_header, nucleotide_seq, original_header, trimmed_name, aa_seq) in original input order
      plus mapping aa_seq -> list of entries (for troubleshooting)
    """
    entries = []  # preserve order: (original_header, trimmed_name, seq, aa)
    for header, seq in records:
        trimmed = header.split('.')[0]
        aa = translate_dna(seq)
        entries.append((header, trimmed, seq, aa))

    # map aa -> list of indices into entries
    aa_to_indices = defaultdict(list)
    for idx, (_, _, _, aa) in enumerate(entries):
        aa_to_indices[aa].append(idx)

    # for each aa, collect unique trimmed names and define base
    aa_base = {}
    for aa, idxs in aa_to_indices.items():
        trimmed_names = sorted({entries[i][1] for i in idxs})
        base = ','.join(trimmed_names)
        aa_base[aa] = base

    # assign new names sequentially within each aa group
    aa_counters = defaultdict(int)
    results = [None] * len(entries)
    for idx, (orig_header, trimmed, seq, aa) in enumerate(entries):
        aa_counters[aa] += 1
        new_header = f"{aa_base[aa]}.{aa_counters[aa]}"
        results[idx] = (new_header, seq, orig_header, trimmed, aa)

    return results, aa_to_indices

def write_fasta_out(path: str, renamed_records: List[Tuple[str,str, str, str, str]]):
    """
    renamed_records: list of tuples (new_header, seq, original_header, trimmed_name, aa_seq)
    """
    with open(path, 'w') as out:
        for new_header, seq, *_rest in renamed_records:
            out.write(f">{new_header}\n")
            out.write(wrap_seq(seq))
            out.write("\n")

def print_histogram(aa_to_indices: Dict[str, List[int]]):
    counts = Counter(len(idxs) for idxs in aa_to_indices.values())
    # sort by N ascending
    print("Peptide occurrence histogram (N sequences -> #peptides):")
    for n in sorted(counts):
        print(f"  {n} -> {counts[n]}")
    total_peptides = sum(counts.values())
    total_seqs = sum(n*counts[n] for n in counts)
    print(f"Total peptides: {total_peptides}; total sequences: {total_seqs}")

def write_troubleshoot_xlsx(xlsx_path: str, renamed_records: List[Tuple[str,str,str,str,str]], aa_to_indices: Dict[str,List[int]]):
    try:
        from openpyxl import Workbook
    except Exception as e:
        raise RuntimeError("Troubleshoot XLSX requires openpyxl. Install with: pip install openpyxl") from e

    # group aa sequences by count
    count_to_aas = defaultdict(list)
    for aa, idxs in aa_to_indices.items():
        count_to_aas[len(idxs)].append(aa)

    wb = Workbook()
    # remove default sheet if we'll create others
    default_sheet = wb.active
    wb.remove(default_sheet)

    # build a map aa -> list of rows (new_name, orig_header, trimmed, seq, aa)
    aa_to_rows = {}
    for new_header, seq, orig_header, trimmed, aa in renamed_records:
        aa_to_rows.setdefault(aa, []).append((new_header, orig_header, trimmed, seq, aa))

    for count, aa_list in sorted(count_to_aas.items()):
        sheet_name = f"count_{count}"
        ws = wb.create_sheet(title=sheet_name)
        ws.append(["new_name", "original_header", "trimmed_name", "nucleotide_seq", "aa_seq"])
        for aa in aa_list:
            rows = aa_to_rows.get(aa, [])
            for row in rows:
                ws.append(list(row))

    wb.save(xlsx_path)

def main(argv=None):
    ap = argparse.ArgumentParser(description="Rename codon-variant duplicates based on translated peptide.")
    ap.add_argument("-i", "--input", required=True, help="Input FASTA (nucleotide)")
    ap.add_argument("-o", "--output", required=True, help="Output FASTA with renamed headers")
    ap.add_argument("--troubleshoot", action="store_true", help="Write XLSX grouping sequences by occurrence counts")
    ap.add_argument("--wrap", type=int, default=60, help="Wrap width for output FASTA sequences (default 60)")
    args = ap.parse_args(argv)

    records = parse_fasta(args.input)
    if not records:
        print("No records found in input.", file=sys.stderr)
        sys.exit(2)

    renamed_records, aa_to_indices = build_new_names_and_records(records)

    # write output FASTA
    write_fasta_out(args.output, renamed_records)

    # histogram
    print_histogram(aa_to_indices)

    if args.troubleshoot:
        xlsx_path = args.output + ".troubleshoot.xlsx"
        try:
            write_troubleshoot_xlsx(xlsx_path, renamed_records, aa_to_indices)
            print(f"Wrote troubleshoot XLSX to: {xlsx_path}")
        except Exception as e:
            print(f"Failed to write troubleshoot XLSX: {e}", file=sys.stderr)
            sys.exit(3)

if __name__ == "__main__":
    main()