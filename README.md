# Host-Pathogen_Spatial

Code accompanying **Chapter 4** of the PhD thesis *Computational Methods for Host-Pathogen Interaction in Spatial Transcriptomics* (Federico Di Marco, University of Milano-Bicocca).

The repository contains the analysis of dual host-pathogen spatial transcriptomics (10x Visium, infected lung) with **InterScale**, a graph neural network with a local (GCN) and a global (Transformer) component trained to reconstruct masked gene expression. The gradient-weighted attention of the global component is converted into directional flows between cell labels, continuous vector fields and source/sink spatial domains. The pathogen labels (P+ / P-) are those predicted by PatHoST (Chapter 3 of the thesis).

> InterScale itself is developed in a separate repository: <https://github.com/theislab/interscale>. This repository only contains the application to the host-pathogen dataset.

## Repository structure

```
notebooks/
  Host-Pathogen.ipynb   # data preparation, model training, flows, domains, statistics, GO enrichment
scripts/
  misty.R               # MISTy multiview analysis of the flow domains (R)
configs/                # InterScale configuration used for training
```

## Workflow

The notebook is organised in sections that correspond to the analyses of Chapter 4:

| Notebook section | What it does | Thesis |
|---|---|---|
| Prepare Data | Integrates the Visium sections, adds the log1p-normalised layer and the PatHoST prediction scores, builds the spatial graph, divides each section into sliding windows (size 3500, overlap 0), random 70/20/10 split, selects 2000 highly variable genes | Methods |
| Bacterial labels | Discretises the PatHoST scores into P+ / P- (K-means with 3 clusters on the scores, the two lower clusters are merged into P-) | Methods |
| Run Local + Global | Trains InterScale (`CombinedModel`, up to 100 epochs, early stopping) | Methods |
| Evaluate | Loads the trained model and extracts the gradient-weighted attention of each window | Methods |
| Clusters and Pathogens | Composite states (host cluster + pathogen label), net flows between composite states, chord diagram, null models (permutation of the pathogen labels within host clusters and across the section) | Results 4.1 |
| Flows | Vector fields and streamplots (section `Train_1`), flow domains (K-means, K = 5) | Results 4.2 |
| Host clusters vs Domains | Sankey diagram, ARI / AMI | Results 4.2 |
| Domains spatial autocorrelation | Join-count statistics (6 nearest neighbours, 999 permutations) | Results 4.2 |
| Net balance vs divergence | Net information balance of the domains, divergence, source/sink roles, Spearman correlation | Results 4.3 |
| Domains interaction | Chord diagram between domains | Results 4.3 |
| Domains markers | Wilcoxon markers per domain and Gene Ontology enrichment, with the selection of the displayed terms | Results 4.4 |

The notebook saves the annotated section `Train_1` (`Data/dummy1.h5ad`, with the columns `flow_cluster` and `bacterial_labels`), which is the input of the R script:

`scripts/misty.R` runs MISTy on `Data/dummy1.h5ad`. The five domain indicators are the targets (used as targets only, `bypass.intra = TRUE`), the pathogen label is the only predictor, arranged in three views (own spot, local neighbours, Gaussian-weighted neighbourhood); the contribution of a spatial scale is the increase in variance explained over the own-spot view, compared with reflected maps of the pathogen as spatial null. It writes the tables and the figure to `misty/` and `Figures/`.

## Requirements

**Python** (notebook): `scanpy`, `anndata`, `squidpy`, `torch`, `scikit-learn`, `scipy`, `statsmodels`, `pandas`, `numpy`, `matplotlib`, `seaborn`, `pycirclize`, `pySankey`, `gseapy`, `joblib` and [`interscale`](https://github.com/theislab/interscale).

**R** (script): `Seurat`, `schard`, `mistyR`, `tidyverse`.



## Data

The data are not included in the repository. The notebook expects the following files in `Data/`:

| File | Content |
|---|---|
| `NBDATA.h5ad` | Visium sections (keys `Train_1`, `Train_2`, `Train_3`, `Train_5`, `Train_6`; `Train_1-3` infected, `Train_5-6` control) |
| `Bact_labels.csv` | PatHoST prediction score of each spot (`Key`, `spotID`, `Prediction`) |
| `Integrated_Clusters.csv` | Host transcriptional clusters (Leiden, mapped to Parenchyma, Inflammation, Low_quality, Airways, SubMucosa) |


The paths in the notebook and the configuration file are those of the cluster used for the analysis; adapt them before running.

## Reproducibility

- Random seeds are fixed (42) for Python, NumPy and PyTorch (`set_full_reproducibility`), for the K-means of the flow domains (42), the join-count permutations (42) and MISTy (42); the K-means of the pathogen labels uses `random_state = 0`.
- The flow domains, the divergence analysis and MISTy refer to a single section (`Train_1`); the chord diagram between composite states pools the windows of the three infected sections.
- GPU training is not bit-wise deterministic across hardware.

## Notes

- The thresholds of the MISTy views (`juxta_thr`, `para_l` in `scripts/misty.R`) are expressed in micrometres.
- The values of the net balance and divergence table used for the scatter plot are transcribed in the notebook from the output of the preceding cells.

## Contact

Federico Di Marco, University of Milano-Bicocca.
