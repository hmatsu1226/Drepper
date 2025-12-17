# Drepper
Dot plot–based framework for repetitive sequence profiling

## Reference

## Requirements
Drepper is written with python.

### Required Python packages
- **numpy** — numerical computations
- **pandas** — tabular data handling
- **biopython** — FASTA sequence parsing (`Bio.SeqIO`)
- **rapidfuzz** — efficient Levenshtein distance computation
- **networkx** — graph construction and analysis
- **python-louvain** — Louvain community detection (`community_louvain`)

## Usage
```
python Drepper.py input.fasta output.tsv
```

## Input File
Protein or nucleotide seuqnces in standard FASTA format.

## Output File
The tool produces a tab-separated values (TSV) file, in which **each row corresponds to one input sequence**.

### Output TSV format

| Column | Name | Description |
|------|------|-------------|
| 1 | `sequence_id` | Sequence identifier (FASTA header) |
| 2 | `complexity` | Number of clusters (Complexity) |
| 3 | `pseudo_cluster` | Cluster ID of pseudo-node (expected to be 0) |
| 4 | `cluster_ids` | Comma-separated list of cluster IDs assigned to the sequence |


### Column details

- **`sequence_id`**  
  Identifier of the input sequence, taken directly from the FASTA header.

- **`complexity`**  
  Total number of clusters detected for the sequence.  
  This value represents the *Complexity* of the repeat structure.

- **`pseudo_cluster`**  
  Cluster ID of pseudo-nodes introduced during graph construction.  
  Under normal conditions, this value should be **0**.

- **`cluster_ids`**  
  A comma-separated string of cluster identifiers (e.g. `1,1,2,3,3,3`) representing the cluster assignment along the sequence.


## Example
```
python Drepper.py examples/example.fa examples/example.tsv
```

## Tutorial including visualizing Dot plot
see example.ipynb

