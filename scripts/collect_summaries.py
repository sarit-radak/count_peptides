#!/usr/bin/env python3
"""
collect_summary.py

Collect per-library trim & count summaries into a single master ./summary.xlsx
Writes two sheets: 'trim' and 'count'. Each sheet has:
  library | <timestamp> | [all numeric columns...] | [all percent columns...]

Requires: pandas, openpyxl
    pip install pandas openpyxl
"""
import os
import glob
import argparse
from datetime import datetime
from typing import Dict, Any, Tuple, List
from pathlib import Path

try:
    import pandas as pd
except Exception:
    pd = None

# ---------- parsing helpers ----------
def parse_count(val):
    """Parse count-like value. Return int/float/string/None."""
    if pd is not None and pd.isna(val):
        return None
    if val is None:
        return None
    # already numeric
    if isinstance(val, (int, float)):
        # convert float integers to int
        if isinstance(val, float) and val.is_integer():
            return int(val)
        return val
    s = str(val).strip()
    if s == "":
        return None
    # percent-looking values handled elsewhere; allow numeric strings
    if s.endswith("%"):
        # leave percent parsing to parse_percent
        return None
    # integer?
    if s.isdigit():
        return int(s)
    try:
        fv = float(s)
        if fv.is_integer():
            return int(fv)
        return fv
    except Exception:
        return s

def parse_percent(val):
    """Return float percent (e.g. '99.67%' -> 99.67), or None if not parseable."""
    if pd is not None and pd.isna(val):
        return None
    if val is None:
        return None
    if isinstance(val, (int, float)):
        return float(val)
    s = str(val).strip()
    if s == "":
        return None
    if s.endswith("%"):
        s = s[:-1].strip()
    try:
        return float(s)
    except Exception:
        return None

def read_metric_xlsx(path: str) -> Dict[str, Tuple[Any, Any]]:
    """
    Read an xlsx file with columns 'metric','count','percent' (case-insensitive).
    Returns dict: metric_label -> (count_value, percent_value)
    """
    if pd is None:
        raise RuntimeError("pandas is required: pip install pandas openpyxl")
    # read first sheet as strings (we'll parse)
    df = pd.read_excel(path, engine="openpyxl", sheet_name=0, dtype=str)
    # normalize columns
    cols = [c.strip().lower() for c in df.columns]
    col_map = {}
    for i, c in enumerate(cols):
        if c == "metric":
            col_map["metric"] = df.columns[i]
        elif c == "count":
            col_map["count"] = df.columns[i]
        elif c == "percent":
            col_map["percent"] = df.columns[i]
    # fallback: first=metric, second=count, third=percent (if names differ)
    if "metric" not in col_map or "count" not in col_map:
        if len(df.columns) >= 2:
            col_map["metric"] = df.columns[0]
            col_map["count"] = df.columns[1]
            if len(df.columns) >= 3:
                col_map["percent"] = df.columns[2]
        else:
            raise ValueError(f"Cannot parse {path}: need at least two columns (metric,count)")
    metrics = {}
    for _, row in df.iterrows():
        metric_label = str(row[col_map["metric"]]).strip()
        raw_count = row[col_map["count"]] if col_map["count"] in row else None
        raw_pct = row[col_map["percent"]] if "percent" in col_map and col_map["percent"] in row else None
        cnt = parse_count(raw_count)
        pct = parse_percent(raw_pct)
        metrics[metric_label] = (cnt, pct)
    return metrics

def file_mtime_iso(path: str) -> str:
    ts = os.path.getmtime(path)
    return datetime.fromtimestamp(ts).isoformat(sep=" ", timespec="seconds")

# ---------- conversion helpers ----------
def metrics_to_flat(metrics: Dict[str, Tuple[Any, Any]]) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    Given metrics dict metric->(count,pct), return two dicts:
      counts: { metric__count: val }
      pcts:   { metric__percent: val }
    """
    counts = {}
    pcts = {}
    for metric, (cnt, pct) in metrics.items():
        counts[f"{metric}__count"] = cnt
        pcts[f"{metric}__percent"] = pct
    return counts, pcts

# ---------- main collect ----------
def collect_all(trim_dir: str = "./2-trimmed", count_dir: str = "./3-counts/summary", out_path: str = "./summary.xlsx"):
    if pd is None:
        raise RuntimeError("This script requires pandas and openpyxl. Install with: pip install pandas openpyxl")

    trim_pattern = os.path.join(trim_dir, "*_trim_summary.xlsx")
    trim_files = sorted(glob.glob(trim_pattern))

    trim_rows = []  # list of dicts for trim sheet
    count_rows = []  # list of dicts for count sheet
    trim_metric_names = set()
    count_metric_names = set()

    for tfile in trim_files:
        b = os.path.basename(tfile)
        if not b.endswith("_trim_summary.xlsx"):
            continue
        lib = b[:-len("_trim_summary.xlsx")]
        try:
            trim_metrics = read_metric_xlsx(tfile)
        except Exception as e:
            print(f"Warning: failed to read {tfile}: {e}")
            continue
        counts_flat, pcts_flat = metrics_to_flat(trim_metrics)
        trim_metric_names.update(k[:-len("__count")] for k in counts_flat.keys())  # store raw metric labels
        row = {}
        row["library"] = lib
        row["trim_timestamp"] = file_mtime_iso(tfile)
        # add all counts and pcts keys (raw)
        # but use flattened names as columns
        row.update(counts_flat)
        row.update(pcts_flat)
        trim_rows.append(row)

        # look for count summary counterpart
        count_file = os.path.join(count_dir, f"{lib}_count_summary.xlsx")
        if os.path.exists(count_file):
            try:
                count_metrics = read_metric_xlsx(count_file)
                c_counts_flat, c_pcts_flat = metrics_to_flat(count_metrics)
                count_metric_names.update(k[:-len("__count")] for k in c_counts_flat.keys())
                crow = {}
                crow["library"] = lib
                crow["count_timestamp"] = file_mtime_iso(count_file)
                crow.update(c_counts_flat)
                crow.update(c_pcts_flat)
                count_rows.append(crow)
            except Exception as e:
                print(f"Warning: failed to read {count_file}: {e}")
                # still add an empty row with library and timestamp None
                count_rows.append({"library": lib, "count_timestamp": None})
        else:
            # no count file: add placeholder row with library and None timestamp
            count_rows.append({"library": lib, "count_timestamp": None})

    if not trim_rows and not count_rows:
        print("No data found. Exiting.")
        return

    # Build trim dataframe:
    trim_df = pd.DataFrame(trim_rows)
    # Ensure consistent columns: library, trim_timestamp, all count cols, then all percent cols
    # Collect all metric labels from existing keys (they are like "Metric name__count")
    all_trim_count_cols = sorted([c for c in trim_df.columns if c.endswith("__count")])
    all_trim_pct_cols = sorted([c for c in trim_df.columns if c.endswith("__percent")])
    ordered_trim_cols = []
    if "library" in trim_df.columns:
        ordered_trim_cols.append("library")
    if "trim_timestamp" in trim_df.columns:
        ordered_trim_cols.append("trim_timestamp")
    ordered_trim_cols += all_trim_count_cols
    ordered_trim_cols += all_trim_pct_cols
    # include any other columns after (defensive)
    for c in trim_df.columns:
        if c not in ordered_trim_cols:
            ordered_trim_cols.append(c)
    trim_df = trim_df.reindex(columns=ordered_trim_cols)

    # Build count dataframe:
    count_df = pd.DataFrame(count_rows)
    all_count_count_cols = sorted([c for c in count_df.columns if c.endswith("__count")])
    all_count_pct_cols = sorted([c for c in count_df.columns if c.endswith("__percent")])
    ordered_count_cols = []
    if "library" in count_df.columns:
        ordered_count_cols.append("library")
    if "count_timestamp" in count_df.columns:
        ordered_count_cols.append("count_timestamp")
    ordered_count_cols += all_count_count_cols
    ordered_count_cols += all_count_pct_cols
    for c in count_df.columns:
        if c not in ordered_count_cols:
            ordered_count_cols.append(c)
    count_df = count_df.reindex(columns=ordered_count_cols)

    # Write to Excel with two sheets
    out_dir = os.path.dirname(out_path) or "."
    os.makedirs(out_dir, exist_ok=True)
    try:
        with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
            # write trim sheet
            if not trim_df.empty:
                trim_df.to_excel(writer, sheet_name="trim", index=False)
            # write count sheet
            if not count_df.empty:
                count_df.to_excel(writer, sheet_name="count", index=False)
        print(f"Wrote summary workbook to {out_path} (trim rows: {len(trim_df)}, count rows: {len(count_df)})")
    except Exception as e:
        print(f"Failed to write {out_path}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Collect trim & count summary xlsx files into a single master summary.xlsx (two sheets)")
    parser.add_argument("--trim-dir", default="./2-trimmed", help="Directory containing *_trim_summary.xlsx files")
    parser.add_argument("--count-dir", default="./3-counts/summary", help="Directory containing *_count_summary.xlsx files")
    parser.add_argument("--out", default="./summary.xlsx", help="Output master Excel file (default ./summary.xlsx)")
    args = parser.parse_args()
    collect_all(trim_dir=args.trim_dir, count_dir=args.count_dir, out_path=args.out)

if __name__ == "__main__":
    main()