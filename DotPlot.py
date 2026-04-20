#!/usr/bin/env python3

from typing import List, Set, Optional
import argparse
import csv
import numpy as np
from Bio import SeqIO
from rapidfuzz.distance import Levenshtein
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch


# ----------------------------
# utility functions
# ----------------------------
def slide_window(x: str, w: int, step: int = 1) -> List[str]:
    if w > len(x) or w <= 0:
        return []
    starts = range(0, len(x) - w + 1, step)
    return [x[s:s + w] for s in starts]


def adist_matrix(strings: List[str]) -> np.ndarray:
    n = len(strings)
    D = np.zeros((n, n), dtype=float)
    for i in range(n):
        si = strings[i]
        for j in range(i + 1, n):
            d = Levenshtein.distance(si, strings[j])
            D[i, j] = d
            D[j, i] = d
    return D


def jaccard_matrix(sets: List[Set[int]]) -> np.ndarray:
    n = len(sets)
    J = np.zeros((n, n), dtype=float)
    for i in range(n):
        si = sets[i]
        for j in range(i, n):
            sj = sets[j]
            if len(si) == 0 and len(sj) == 0:
                v = 1.0
            else:
                inter = len(si & sj)
                uni = len(si | sj)
                v = (inter / uni) if uni > 0 else 0.0
            J[i, j] = v
            J[j, i] = v
    return J


# ----------------------------
# data loading
# ----------------------------
def load_target_sequence(fasta_path: str, target: str) -> str:
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if len(records) == 0:
        raise FileNotFoundError(f"No records in FASTA: {fasta_path}")

    tids = [rec.id for rec in records]
    seqidx = [i for i, tid in enumerate(tids) if tid.startswith(target)]

    if len(seqidx) == 0:
        raise ValueError(f"Target prefix not found in FASTA: {target}")
    if len(seqidx) > 1:
        print(f"Warning: multiple FASTA records matched target prefix '{target}'. Using the first one.")

    return str(records[seqidx[0]].seq)


def load_group_vector(drepper_path: str, target: str, group_col: int = 3) -> List[int]:
    """
    group_col:
        0-based column index.
        3 means the 4th column.
    """
    grp_vector = None

    with open(drepper_path) as f:
        reader = csv.reader(f, delimiter="\t")
        for i, row in enumerate(reader):
            if len(row) <= group_col:
                continue
            gene_id = row[0]
            if gene_id.startswith(target):
                grp_vector = [int(x) for x in row[group_col].split(",") if x != ""]
                break

    if grp_vector is None:
        raise ValueError(f"Target prefix not found in drepper TSV: {target}")

    return grp_vector


# ----------------------------
# main plotting
# ----------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Generate annotated heatmap for a target sequence."
    )
    parser.add_argument("fasta_path", help="Input FASTA file")
    parser.add_argument("drepper_path", help="Input drepper TSV file")
    parser.add_argument("target", help="Target prefix for sequence/gene ID")
    parser.add_argument("output_png", help="Output PNG filename")

    parser.add_argument("--start", type=int, default=None,
                        help="Start index for slicing (0-based, inclusive)")
    parser.add_argument("--end", type=int, default=None,
                        help="End index for slicing (0-based, exclusive)")
    parser.add_argument("--w", type=int, default=10,
                        help="Sliding window size")
    parser.add_argument("--stepsize", type=int, default=1,
                        help="Sliding window step size")
    parser.add_argument("--thresh", type=float, default=0.2,
                        help="Threshold for binarizing similarity matrix")
    parser.add_argument("--figsize", type=float, default=12,
                        help="Figure size")
    parser.add_argument("--group-col", type=int, default=3,
                        help="0-based TSV column index containing comma-separated group labels (default: 3)")

    args = parser.parse_args()

    fasta_path = args.fasta_path
    drepper_path = args.drepper_path
    target = args.target
    output_png = args.output_png

    # load data
    seq = load_target_sequence(fasta_path, target)
    grp_vector = load_group_vector(drepper_path, target, group_col=args.group_col)

    # sliding windows
    subseqs = slide_window(seq, w=args.w, step=args.stepsize)
    n = len(subseqs)

    if n == 0:
        raise ValueError("No subsequences generated. Check sequence length, window size, and step size.")

    # check length consistency
    if len(grp_vector) != n and len(grp_vector) != n + 1 and len(grp_vector) != n - 1:
        print(f"Warning: length mismatch detected: len(subseqs)={n}, len(grp_vector)={len(grp_vector)}")

    # original code had grp_vector[1:], so keep that behavior when appropriate
    if len(grp_vector) == n + 1:
        grp_use = grp_vector[1:]
    elif len(grp_vector) == n:
        grp_use = grp_vector
    else:
        min_len = min(n, len(grp_vector))
        print(f"Warning: truncating data to min length = {min_len}")
        subseqs = subseqs[:min_len]
        grp_use = grp_vector[:min_len]
        n = min_len

    # range slicing
    start = 0 if args.start is None else args.start
    end = n if args.end is None else args.end

    if start < 0 or end < 0 or start >= end or end > n:
        raise ValueError(f"Invalid range: start={start}, end={end}, valid range is 0 <= start < end <= {n}")

    subseqs_sub = subseqs[start:end]
    grp_sub = grp_use[start:end]

    # distance matrix
    D = adist_matrix(subseqs_sub)

    # similarity matrix
    W = np.exp(-D)
    np.fill_diagonal(W, 0.0)

    W_bin = np.where(W > args.thresh, 1.0, 0.0)
    np.fill_diagonal(W_bin, 1.0)

    # annotation colors
    unique_grps = np.unique(grp_sub)
    cmap_bar = plt.get_cmap("tab20", len(unique_grps))
    color_map = {g: cmap_bar(i) for i, g in enumerate(unique_grps)}

    row_colors = np.array([color_map[g] for g in grp_sub])
    col_colors = row_colors

    cmap_bin = ListedColormap(["ghostwhite", "black"])
    norm_bin = BoundaryNorm([-0.5, 0.5, 1.5], cmap_bin.N, clip=True)

    # plotting
    fig = plt.figure(figsize=(args.figsize*1.1, args.figsize), constrained_layout=True)
    gs = fig.add_gridspec(
        2, 2,
        width_ratios=(0.15, 3.0),
        height_ratios=(0.15, 3.0),
        wspace=0.05,
        hspace=0.05
    )

    ax_col = fig.add_subplot(gs[0, 1])
    ax_row = fig.add_subplot(gs[1, 0])
    ax_main = fig.add_subplot(gs[1, 1])

    # top annotation
    ax_col.imshow(col_colors[np.newaxis, :, :], aspect="auto")
    ax_col.set_axis_off()

    # left annotation
    ax_row.imshow(row_colors[:, np.newaxis, :], aspect="auto")
    ax_row.set_axis_off()

    # main heatmap
    ax_main.imshow(W_bin, cmap=cmap_bin, norm=norm_bin, aspect="auto", interpolation="nearest")
   
    length = end - start
    nticks = 10
    tick_pos = np.linspace(0, length - 1, nticks, dtype=int)
    tick_labels = start + tick_pos

    ax_main.set_xticks(tick_pos)
    ax_main.set_yticks(tick_pos)
    ax_main.set_xticklabels(tick_labels)
    ax_main.set_yticklabels(tick_labels)

    ax_main.set_title(f"{target}  [{start}:{end}]")
    ax_main.set_xlabel("Window index")
    ax_main.set_ylabel("Window index")

    handles = [
        Patch(facecolor=color_map[g], edgecolor="k", label=f"Cluster {g}")
        for g in unique_grps
    ]
    ax_main.legend(handles=handles, bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0.)

    plt.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved heatmap to: {output_png}")


if __name__ == "__main__":
    main()