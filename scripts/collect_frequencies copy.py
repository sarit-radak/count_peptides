#!/usr/bin/env python3
"""
collect_frequencies.py

Scan ./3-counts for CSV files, extract Peptide, Translation and Frequency columns,
merge all files on Peptide, and write ./3-counts/<enclosing-dir-name>_all.csv

Output columns:
Peptide, Translation, <filename1>, <filename2>, ...
where each filename column contains the Frequency value from that file
(filename without .csv).
"""
from pathlib import Path
import csv
import pandas as pd
import sys
import logging

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def detect_delimiter(path: Path, nbytes: int = 4096) -> str:
    """Try to detect delimiter using csv.Sniffer on a sample of the file."""
    sample = path.read_text(errors="ignore")[:nbytes]
    if not sample:
        return ","
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=[",", "\t", ";", "|"])
        return dialect.delimiter
    except csv.Error:
        # fallback: prefer tab if tabs present, else comma
        if "\t" in sample:
            return "\t"
        return ","


def find_column(cols, target):
    """
    Return the actual column name from cols that matches target ignoring case and whitespace,
    or None if not found.
    """
    target_key = target.lower().strip()
    for c in cols:
        if c is None:
            continue
        if c.lower().strip() == target_key:
            return c
    return None


def read_relevant_columns(path: Path):
    """Read CSV and return DataFrame with columns Peptide, Translation (if present), Frequency (if present)."""
    delim = detect_delimiter(path)
    try:
        df = pd.read_csv(path, sep=delim, dtype=str, keep_default_na=False)
    except Exception as e:
        logging.error(f"Failed to read {path}: {e}")
        return None

    # strip column names
    df.rename(columns={c: c.strip() for c in df.columns}, inplace=True)
    cols = list(df.columns)

    peptide_col = find_column(cols, "Peptide")
    freq_col = find_column(cols, "Frequency")
    trans_col = find_column(cols, "Translation")

    if peptide_col is None:
        logging.warning(f"Skipping {path.name}: no Peptide column found.")
        return None
    if freq_col is None:
        logging.warning(f"Skipping {path.name}: no Frequency column found.")
        return None

    # Keep only needed columns (Translation may be absent)
    keep = [peptide_col, freq_col]
    if trans_col:
        keep.append(trans_col)

    sub = df[keep].copy()
    # Normalize column names to standard names
    rename_map = {peptide_col: "Peptide", freq_col: "Frequency"}
    if trans_col:
        rename_map[trans_col] = "Translation"
    sub.rename(columns=rename_map, inplace=True)

    # Ensure Peptide is string
    sub["Peptide"] = sub["Peptide"].astype(str)
    # Frequency keep as string to avoid parsing weird formats; convert later if needed
    sub["Frequency"] = sub["Frequency"].replace("", pd.NA)

    return sub


def main():
    base_dir = Path("./3-counts")
    if not base_dir.exists() or not base_dir.is_dir():
        logging.error(f"Directory {base_dir} does not exist or is not a directory.")
        sys.exit(1)

    csv_files = sorted(base_dir.glob("*.csv"))
    if not csv_files:
        logging.error(f"No .csv files found in {base_dir.resolve()}.")
        sys.exit(1)

    merged = None
    translation_reference = None
    translations_conflicts = []

    for path in csv_files:
        logging.info(f"Processing {path.name} ...")
        sub = read_relevant_columns(path)
        if sub is None:
            continue

        name = path.stem  # filename without .csv
        # set index by Peptide
        sub = sub.set_index("Peptide", drop=False)

        # handle Translation: keep only one 'Translation' column in final result
        if "Translation" in sub.columns:
            if translation_reference is None:
                # record translations for reference
                translation_reference = sub[["Peptide", "Translation"]].set_index("Peptide")["Translation"]
            else:
                # compare for conflicts; record any differences
                # align and compare where both exist
                merged_trans = pd.concat([translation_reference, sub[["Peptide", "Translation"]].set_index("Peptide")["Translation"]], axis=1)
                merged_trans.columns = ["ref", "new"]
                conflicts = merged_trans[(merged_trans["ref"].notna()) & (merged_trans["new"].notna()) & (merged_trans["ref"] != merged_trans["new"])]
                if not conflicts.empty:
                    logging.warning(f"Translation conflicts found in {path.name} for {len(conflicts)} peptides. Keeping first-seen translation.")
                    translations_conflicts.append((path.name, conflicts.head(5)))  # store a small sample

        # extract the Frequency series and name it after the file
        freq_series = sub["Frequency"].copy()
        freq_series.name = name

        # Merge into master table (outer join so we keep any peptide)
        if merged is None:
            # start with Peptide and (optional) Translation and the frequency column
            merged = pd.DataFrame(index=sub["Peptide"])
            merged.index.name = "Peptide"
            merged = merged.join(freq_series, how="outer")
            if "Translation" in sub.columns:
                # ensure translation_reference set (done above)
                merged = merged.merge(sub[["Peptide", "Translation"]].set_index("Peptide"), left_index=True, right_index=True, how="left")
            # re-order to have Translation after index later
            merged.reset_index(inplace=True)
            merged.set_index("Peptide", inplace=True)
        else:
            merged = merged.join(freq_series, how="outer")

    if merged is None:
        logging.error("No valid CSVs processed. Exiting.")
        sys.exit(1)

    # If we have a translation_reference, set Translation column based on it (keep first-seen)
    if translation_reference is not None:
        merged_trans = merged.index.to_series().map(translation_reference)
        if "Translation" in merged.columns:
            merged["Translation"] = merged["Translation"].where(merged["Translation"].notna(), merged_trans)
        else:
            merged["Translation"] = merged_trans

    # Reorder columns: Peptide, Translation, then frequency columns sorted by name
    result = merged.reset_index()
    freq_cols = [c for c in result.columns if c not in ("Peptide", "Translation")]
    freq_cols_sorted = sorted(freq_cols)
    cols_out = ["Peptide"]
    if "Translation" in result.columns:
        cols_out.append("Translation")
    cols_out.extend(freq_cols_sorted)
    result = result[cols_out]

    # ✅ Sort alphabetically by Peptide (ascending)
    result.sort_values(by="Peptide", ascending=True, inplace=True)

    # ✅ Build output filename from parent directory of ./3-counts
    parent_name = base_dir.parent.resolve().name
    out_name = f"{parent_name}_all.csv"
    out_path = base_dir / out_name

    # Write CSV
    result.to_csv(out_path, index=False)
    logging.info(f"Wrote merged output to {out_path}")

    # If conflicts were found, print a short summary
    if translations_conflicts:
        logging.warning("Translation conflicts were detected in some files. Samples (up to 5 rows per file) are shown below:")
        for fname, sample in translations_conflicts:
            logging.warning(f"File: {fname}")
            # show up to 5 sample rows from the stored conflicts
            for idx, row in sample.head(5).iterrows():
                logging.warning(f" Peptide: {idx} ref='{row['ref']}' new='{row['new']}'")

    logging.info("Done.")


if __name__ == "__main__":
    main()