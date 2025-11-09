#!/usr/bin/env python3
"""
organize_fastq_pass.py

- Reads ./barcodes.xlsx which should contain columns: Barcode, Well, Library Name
- For rows where Library Name is non-empty, rename ./fastq_pass/<Barcode> -> ./fastq_pass/<Library Name>
  (sanitizes Library Name for filesystem usage; will avoid collisions by appending suffixes if necessary)
- Move all remaining directories inside ./fastq_pass that were not renamed into ./fastq_pass_ignore

Dependencies:
    pandas (for reading Excel)
Usage:
    python3 organize_fastq_pass.py
"""
from pathlib import Path
import shutil
import logging
import pandas as pd
import re
import sys

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def sanitize_name(name: str) -> str:
    """Make a filesystem-safe name: strip whitespace and replace path separators and illegal chars."""
    name = name.strip()
    # Replace path separators and some other problematic characters with underscore
    name = re.sub(r'[\\/]+', '_', name)
    # Optionally remove control characters
    name = "".join(ch for ch in name if ch.isprintable())
    # Replace characters that are commonly problematic on filenames
    name = re.sub(r'[:\*\?"<>\|]', '_', name)
    # Collapse runs of spaces into single underscores
    name = re.sub(r'\s+', '_', name)
    return name


def unique_destination(base_dir: Path, desired_name: str, max_tries: int = 100) -> Path:
    """Return a destination path that does not already exist by appending numeric suffixes if needed."""
    dest = base_dir / desired_name
    if not dest.exists():
        return dest
    # If exists, try appending _1, _2, ...
    for i in range(1, max_tries + 1):
        candidate = base_dir / f"{desired_name}_{i}"
        if not candidate.exists():
            return candidate
    raise RuntimeError(f"Could not find a unique destination name for {desired_name} after {max_tries} tries.")


def main():
    cwd = Path.cwd()
    print (cwd)
    excel_path = cwd / "barcodes.xlsx"
    fastq_dir = cwd / "fastq_pass"
    ignore_dir = cwd / "fastq_pass_ignore"

    if not excel_path.exists():
        logging.error(f"Excel file not found: {excel_path}")
        sys.exit(1)
    if not fastq_dir.exists() or not fastq_dir.is_dir():
        logging.error(f"fastq_pass directory not found: {fastq_dir}")
        sys.exit(1)

    # Read Excel
    try:
        df = pd.read_excel(excel_path, dtype=str)
    except Exception as e:
        logging.error(f"Failed to read {excel_path}: {e}")
        sys.exit(1)
    
    # Normalize column names and ensure expected columns exist
    df.rename(columns={c: c.strip() for c in df.columns}, inplace=True)
    expected_cols = ["Barcode", "Library Name"]
    if "Barcode" not in df.columns or "Library Name" not in df.columns:
        logging.error(f"Excel must contain columns: {expected_cols}. Found: {list(df.columns)}")
        sys.exit(1)

    # Build mapping barcode -> library name (only where Library Name is non-empty)
    df["Barcode"] = df["Barcode"].astype(str).str.strip()
    df["Library Name"] = df["Library Name"].fillna("").astype(str).str.strip()

    rename_map = {}
    for _, row in df.iterrows():
        barcode = row["Barcode"]
        libname = row["Library Name"]
        if libname:
            rename_map[barcode] = sanitize_name(libname)

    # Keep track of which destination names we created (final folder names in fastq_pass)
    created_dest_names = set()

    # Perform renames
    for barcode, desired_name in rename_map.items():
        src = fastq_dir / barcode
        if not src.exists():
            logging.warning(f"Source folder for barcode not found, skipping: {src}")
            continue
        if not src.is_dir():
            logging.warning(f"Source path exists but is not a directory, skipping: {src}")
            continue

        # Determine unique destination
        dest = unique_destination(fastq_dir, desired_name)
        try:
            logging.info(f"Renaming: {src.name} -> {dest.name}")
            src.rename(dest)
            created_dest_names.add(dest.name)
        except Exception as e:
            logging.error(f"Failed to rename {src} -> {dest}: {e}")

    # Ensure ignore directory exists
    try:
        ignore_dir.mkdir(exist_ok=True)
    except Exception as e:
        logging.error(f"Failed to create ignore directory {ignore_dir}: {e}")
        sys.exit(1)

    # Move remaining directories (except the ignore folder itself) into ignore_dir
    # Do not move items that are in created_dest_names; also ignore files (only operate on directories)
    for item in fastq_dir.iterdir():
        if item.name == ignore_dir.name:
            continue
        if not item.is_dir():
            # skip files
            continue
        if item.name in created_dest_names:
            # this one was renamed to a library name; skip
            continue
        # Move this directory into ignore_dir
        target = ignore_dir / item.name
        # If target already exists, try to find a unique name
        if target.exists():
            # attempt to find unique target by suffix
            i = 1
            while True:
                candidate = ignore_dir / f"{item.name}_{i}"
                if not candidate.exists():
                    target = candidate
                    break
                i += 1
        try:
            logging.info(f"Moving to ignore: {item.name} -> {target}")
            shutil.move(str(item), str(target))
        except Exception as e:
            logging.error(f"Failed to move {item} -> {target}: {e}")

    logging.info("Done.")


if __name__ == "__main__":
    main()