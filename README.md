# README: Pathway Enrichment Network Visualization

This R script aggregates enriched pathways across multiple datasets, filters by statistical significance, calculates gene set overlaps between pathways, and visualizes the resulting pathway network using igraph. Output includes a labeled network graph and a colorbar based on -log(p-value) color intensity.

---

## Purpose

- Load enrichment results from multiple text files.
- Filter for statistically significant pathways (p-value < 0.01).
- Consolidate unique pathways across datasets.
- Retrieve associated genes from a pathway reference file.
- Compute pairwise overlap scores between pathway gene sets.
- Build an undirected network graph based on gene set overlap.
- Visualize the network, coloring nodes by enrichment significance and scaling them by pathway overlap.
- Export outputs including a tab-delimited file of pathways, a PDF network plot, and a PNG colorbar.

---

## Input Files

- `enrichRpathways.txt`: Tab-delimited reference file listing pathways and their associated genes.
- Folder of pathway result files (`*.txt`), each containing enrichment results from separate datasets.

---

## Key Parameters

| Parameter         | Description                                                    |
|------------------|----------------------------------------------------------------|
| `linkageThresh`  | Minimum overlap score to create an edge between two pathways   |
| `p < 0.01`       | Filtering threshold for enrichment significance                |
| `vertex.size`    | Proportional to pathway gene overlap fraction                  |
| `vertex.color`   | Based on -log(p-value), mapped to a custom red-white palette   |
| `output.txt`     | Consolidated table of top pathways and overlap statistics      |
| `output.pdf`     | Final network graph of pathways                                |
| `colorbar.png`   | Color scale legend for pathway node significance               |

---

## Output Files

- **`output.txt`** – Consolidated and ranked list of pathways (with overlap and scores).
- **`output.pdf`** – Network plot with nodes sized by gene set overlap and colored by p-value.
- **`colorbar.png`** – Gradient legend used in the network node coloring.

---

## Dependencies

The script requires the following R packages:

```r
library(igraph)
library(bayestestR)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(fields)
