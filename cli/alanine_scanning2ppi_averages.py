#!/usr/bin/env python3

import glob
import os
import re
import argparse
import zipfile
from contextlib import contextmanager
from io import TextIOWrapper
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# -----------------------------
# Parse arguments
# -----------------------------
parser = argparse.ArgumentParser(
    description="Compute average PPI probabilities from alanine scanning files and generate a barplot."
)

parser.add_argument(
    "-i", "--input_pattern",
    default="Query__Partner_AlaScan*__query_as_peptide.txt",
    help=(
        "Input glob pattern, directory, or zip file. Directories may contain "
        "unpacked .txt files and/or zip files. Default: %(default)s"
    )
)

parser.add_argument(
    "-r", "--range",
    nargs=2,
    type=int,
    metavar=("START", "END"),
    help="Position range (1-based, inclusive) to average over (default: full sequence)"
)

parser.add_argument(
    "-o", "--output",
    default="alanine_scan_barplot.png",
    help="Output plot filename (default: %(default)s)"
)

parser.add_argument(
    "-w", "--window_size",
    type=int,
    default=30,
    help="Alanine scan window size (default: %(default)s)"
)

parser.add_argument(
    "--xtick-fontsize",
    type=float,
    default=5,
    help="Font size for x-axis tick labels (default: %(default)s)"
)

parser.add_argument(
    "-m", "--model",
    choices=("peptide", "receptor"),
    default="peptide",
    help=(
        "Model results to use for plotting: peptide selects "
        "*__query_as_peptide.txt and receptor selects *__query_as_receptor.txt. "
        "Default: %(default)s"
    )
)

args = parser.parse_args()


MODEL_SUFFIXES = {
    "peptide": "__query_as_peptide.txt",
    "receptor": "__query_as_receptor.txt",
}
MODEL_LABELS = {
    "peptide": "ligand-centric",
    "receptor": "receptor-centric",
}

def collect_zip_prediction_files(zip_path, source_prefix=None):
    files = []
    with zipfile.ZipFile(zip_path) as archive:
        for info in archive.infolist():
            if info.is_dir() or not info.filename.endswith(".txt"):
                continue
            source_name = info.filename
            if source_prefix:
                source_name = f"{source_prefix}/{info.filename}"
            files.append(("zip", zip_path, info.filename, source_name))
    return files


def collect_input_files(input_pattern):
    files = []

    if os.path.isdir(input_pattern):
        for root, _, filenames in os.walk(input_pattern):
            for filename in filenames:
                path = os.path.join(root, filename)
                if filename.endswith(".txt"):
                    files.append(("file", path, None, path))
                elif zipfile.is_zipfile(path):
                    files.extend(collect_zip_prediction_files(path, path))
        return sorted(files, key=lambda item: item[3])

    if zipfile.is_zipfile(input_pattern):
        return sorted(
            collect_zip_prediction_files(input_pattern),
            key=lambda item: item[3],
        )

    for path in sorted(glob.glob(input_pattern)):
        files.append(("file", path, None, path))
    return files


def filter_files_by_model(files, model):
    suffix = MODEL_SUFFIXES[model]
    return [file_info for file_info in files if file_info[3].endswith(suffix)]


@contextmanager
def open_input_file(file_info):
    input_type, path, member_name, _ = file_info
    if input_type == "zip":
        with zipfile.ZipFile(path) as archive:
            with archive.open(member_name) as raw:
                yield TextIOWrapper(raw, encoding="utf-8", newline="")
    else:
        with open(path, newline="") as handle:
            yield handle


# -----------------------------
# Collect files
# -----------------------------
files = filter_files_by_model(collect_input_files(args.input_pattern), args.model)

if not files:
    raise ValueError(f"No {args.model} model input .txt files found")

print(f"Found {len(files)} files")


# -----------------------------
# Process files
# -----------------------------
all_probs = []
scan_labels = []

for file_info in files:
    f = file_info[3]
    match = re.search(r"AlaScan(\d+)", f)
    if not match:
        continue

    scan_start = int(match.group(1))
    scan_end = scan_start + args.window_size - 1
    label = f"{scan_start}-{scan_end}"

    with open_input_file(file_info) as handle:
        df = pd.read_csv(handle, sep="\t", comment="#")
    df.columns = ["Position", "AminoAcid", "Probability"]
    df["Position"] = df["Position"].astype(int)
    df["Probability"] = df["Probability"].astype(float)

    if args.range:
        r_start, r_end = args.range
        df_sub = df[
            (df["Position"] >= r_start) &
            (df["Position"] <= r_end)
        ].copy()

        if df_sub.empty:
            raise ValueError(
                f"No residues found in range {r_start}-{r_end} for file: {f}"
            )
    else:
        df_sub = df.copy()

    mean_prob = df_sub["Probability"].mean()

    all_probs.append(mean_prob)
    scan_labels.append(label)

# -----------------------------
# Sort by scan start
# -----------------------------
scan_data = list(zip(scan_labels, all_probs))
scan_data.sort(key=lambda x: int(x[0].split("-")[0]))

labels_sorted = [x[0] for x in scan_data]
values_sorted = [x[1] for x in scan_data]


# -----------------------------
# Plot
# -----------------------------
dpi = 300
fig_width = 3000 / dpi
fig_height = 1200 / dpi

plt.figure(figsize=(fig_width, fig_height), dpi=dpi)

x = np.arange(len(values_sorted))

plt.bar(x, values_sorted, color="blue", width=0.8)

plt.xticks(x, labels_sorted, rotation=45, ha='right', fontsize=args.xtick_fontsize)

plt.xlabel("Alanine scan window (positions)")
plt.ylabel("Avg PPI probability")

if args.range:
    plt.title(f"Average {MODEL_LABELS[args.model]} probabilities (over positions {args.range[0]}–{args.range[1]})")
else:
    plt.title(f"Average {MODEL_LABELS[args.model]} probabilities (over full sequence)")

plt.grid(axis='y', linestyle="--", alpha=0.4)

plt.tight_layout()
plt.savefig(args.output, dpi=dpi, bbox_inches="tight")
plt.close()

print(f"Bar plot saved to: {os.path.abspath(args.output)}")
