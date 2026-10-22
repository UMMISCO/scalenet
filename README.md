# scalenet

<!-- badges: start -->
[![R-CMD-check](https://github.com/UMMISCO/scalenet/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/UMMISCO/scalenet/actions/workflows/R-CMD-check.yaml)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](LICENSE)
<!-- badges: end -->

Network reconstruction from high-dimensional data by spectral decomposition,
with the **Spectral Consensus Strategy** that fuses several embedded
reconstruction methods into one network.

> Affeldt S, Sokolovska N, Prifti E, Zucker J-D. *Spectral consensus strategy for
> accurate reconstruction of large biological networks.* BMC Bioinformatics 2016;
> 17(Suppl 16):493. [PMC5249011](https://pmc.ncbi.nlm.nih.gov/articles/PMC5249011/)

This is the reference implementation. [**xscalenet**](https://github.com/UMMISCO/xscalenet)
is a Rust port of the same core path, validated against this package's own
fixtures: on the datasets the two have been compared on they produce **identical
edge sets**, before and after the fix described below. It reconstructs a
135-variable matrix with ARACNE about 2.5× faster (0.34 s against 0.85 s), and
keeps every subgraph in memory rather than writing it to disk and reading it
back as this package does.

## Installation

```r
# CRAN
install.packages(c("doSNOW", "doParallel", "foreach", "bnlearn",
                   "fpc", "e1071", "data.table", "infotheo", "plyr", "igraph"))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("minet")

remotes::install_github("UMMISCO/scalenet")
```

## Usage

```r
library(scalenet)

pop2mat <- read.delim(system.file("extdata", "pop2mat.txt", package = "scalenet"))

scs(workspaceDir = "scalenet_results",
    argInData    = pop2mat,
    argReconsMeth = c("aracne", "bayes_hc"),
    argReconsMethInfo = list(aracne   = list(ort = "n", eweight = "epresenceScore"),
                             bayes_hc = list(ort = "y", eweight = "ecorr")),
    argEmbReconsParam = list(aracne   = list(estimator = "mi.mm", epsilon = 0.001),
                             bayes_hc = list(score = "bde", restart = 0),
                             varPerc  = 0.2),
    argPresFreqThresh = c(0.3, 0.8),
    argDiscretize = TRUE, argVerbose = TRUE)
```

Edge lists land under
`scalenet_results/consensusGraph/globalNet_presFreq<t>/edgesList.txt`, with
columns `x  y  epresenceScore  epresence  eorientScore  eorient  ecorr`.

`Readme.Rmd` carries a longer worked example, including loading an annotation
and plotting the result with igraph.

## Reproducibility

Up to version 1.2.3 the reconstruction depended on the **layout** of the input:
shuffling the columns of a matrix changed the network (ARACNE edge Jaccard
0.28–0.45 across permutations on four genus presence/absence matrices), and so
did renaming the variables, while the data itself was identical.

`Lrw = I − D⁻¹W` is not symmetric, yet `computeEigenVectVal()` decomposed it with
`eigen(symmetric = TRUE)`, which reads only the **lower triangle** — so permuting
the columns changed the matrix being decomposed, not merely its row and column
order. `LaplacianRW()` now builds the symmetric normalised Laplacian
`Lsym = I − D^-1/2 W D^-1/2` and the eigenvectors are mapped back to `Lrw`'s with
`v = D^-1/2 u`. Separately, `bnlearn::hc()` visits candidate moves in column
order, so `rMethod.hc()` now orders variables canonically by their own data and
fuses several random orderings. Both changes bring the Jaccard to **1.0000**.

`scalenet(argLegacyLrw = TRUE)` reproduces the old behaviour, so analyses
published before the fix can be regenerated. Details in [NEWS.md](NEWS.md).

## License

GPL-3. See [LICENSE](LICENSE).
