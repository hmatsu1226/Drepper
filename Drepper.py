import sys
from typing import List, Tuple, Dict, Set
import numpy as np
import pandas as pd
from Bio import SeqIO
from rapidfuzz.distance import Levenshtein
import networkx as nx
from community import community_louvain

fasta_path = sys.argv[1]
output = sys.argv[2]

#parameter
w = 10
stepsize = 1
max_pairs = 50

#get subsequences with sliding window
def slide_window(x: str, w: int, step: int = 1) -> List[str]:
    if w > len(x) or w <= 0:
        return []
    starts = range(0, len(x) - w + 1, step)
    return [x[s:s + w] for s in starts]

#get Leivenshtein distance matrix among subsequences
def adist_matrix(strings: List[str]) -> np.ndarray:
    n = len(strings)
    D = np.zeros((n, n), dtype=float)
    for i in range(n):
        D[i, i] = 0.0
        si = strings[i]
        for j in range(i + 1, n):
            d = Levenshtein.distance(si, strings[j])
            D[i, j] = d
            D[j, i] = d
    return D

#Jaccard index based similarity matrix
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

#read FASTA file
records = list(SeqIO.parse(fasta_path, "fasta"))
if len(records) == 0:
    raise FileNotFoundError(f"No records in FASTA: {fasta_path}")

#names of sequences
tids = [rec.id for rec in records]

N = len(records)

result = np.zeros((N, 3), dtype=object)
index = tids

for seqidx in range(N):
    print(seqidx)

    # subseqs
    seq = str(records[seqidx].seq)
    subseqs = slide_window(seq, w=w, step=stepsize)
    n = len(subseqs)

    #
    D = adist_matrix(subseqs)

    #
    W = np.exp(-D)
    np.fill_diagonal(W, 0.0)
    thresh = 0.2
    W_bin = np.where(W > thresh, 1.0, 0.0)
    np.fill_diagonal(W_bin, 1.0)

    #skip seq if it is almost diagonal W_bin
    if W_bin.sum() < (n + 10):
        result[seqidx, 0] = 1
        result[seqidx, 1] = 0
        result[seqidx, 2] = ",".join(["0"] * n)
        continue

    #get features of each location
    features: List[Set[int]] = []
    features.append(set()) #add pseudo-node
    for i in range(W_bin.shape[0]):
        idxs = np.flatnonzero(W_bin[i] == 1)

        if len(idxs) > max_pairs:
            idxs = sorted(idxs, key=lambda v: (abs(v - i), v))

        if len(idxs) < 2:
            features.append(set())
            continue

        diffs: Set[int] = set()
        for j in range(min(len(idxs)-1,max_pairs-1)):
            for k in range(j+1, min(len(idxs),max_pairs)):
                diffs.add(abs(int(idxs[j]) - int(idxs[k])))
        features.append(diffs)

    #Similarity matrix with Jaccard index
    A = jaccard_matrix(features)

    #
    thresh = 0.5
    A_bin = np.where(A > thresh, 1.0, 0.0)
    np.fill_diagonal(A_bin, 0.0)
    G = nx.from_numpy_array(A_bin)

    # Louvain 
    parts = community_louvain.best_partition(G, weight="weight", random_state=0)
    
    # label change
    cid0 = parts[0]
    new_parts = {}
    for node, cid in parts.items():
        if cid == cid0:
            new_parts[node] = 0
        elif cid == 0:
            new_parts[node] = cid0
        else:
            new_parts[node] = cid
    parts = new_parts

    grp_vector = [parts.get(i, -1) for i in range(n+1)]
    grp_series = pd.Series(grp_vector)
    max_grp = int(grp_series.max()) if len(grp_series) else 0

    # output
    grp_csv = ",".join(map(str, grp_vector[1:]))

    result[seqidx, 0] = max_grp + 1
    result[seqidx, 1] = grp_vector[0]
    result[seqidx, 2] = grp_csv

df = pd.DataFrame(result, index=index, columns=[0, 1, 2])
df.to_csv(output, sep="\t", header=False)
print(f"[DONE] wrote: {output}")
