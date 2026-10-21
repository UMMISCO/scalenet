# scalenet news

## Unreleased — reconstruction no longer depends on the input layout

Shuffling the columns of an input matrix changed the reconstructed network, and
renaming the variables changed it again, while the data stayed identical. On four
genus presence/absence matrices the ARACNE edge Jaccard across 10 column
permutations was **0.28–0.45**; it is now **1.0000**.

### The spectral step decomposed the wrong matrix

`LaplacianRW()` built `decompose.mat = diag(n) - D^-1 %*% W`, the random-walk
Laplacian `Lrw`, which is **not symmetric** — and `computeEigenVectVal()` then
called `eigen(decompose.mat, symmetric = TRUE)`. That call reads only the **lower
triangle**, so entry `(i, j)` was taken as `-W[i,j]/D[i]` with `i` the *larger*
index. Permuting the input columns therefore changed the matrix being
decomposed, not merely its row and column order: the eigenvectors changed, the
`G+,m` / `G-,m` variable subsets changed, and every reconstructed edge changed
with them.

Verified rather than assumed:
`eigen(matrix(c(1,99,7,1), 2, 2, byrow = TRUE), symmetric = TRUE)` returns the
eigenvalues of the lower-triangle symmetrisation (`8, -6`), not of the upper one
(`100, -98`).

`LaplacianRW()` now decomposes the symmetric normalised Laplacian
`Lsym = I - D^-1/2 W D^-1/2`, which is symmetric by construction, and
`computeEigenVectVal()` maps its eigenvectors back to `Lrw`'s with
`v = D^-1/2 u` (same eigenvalues, since `Lrw = D^-1/2 Lsym D^1/2`). Checked
numerically against a direct multiplication by the asymmetric `Lrw`: relative
residual `2e-15`, where the old path gives `1.5e-1` — it was not computing `Lrw`'s
eigenvectors at all. Eigenvector signs are canonicalised, since a sign flip swaps
an eigenvector's `G+` and `G-` subsets.

`scalenet(argLegacyLrw = TRUE)` restores the old behaviour so previously
published results can be regenerated.

**Cached matrices are keyed by variant.** `LaplacianRW()` and
`computeEigenVectVal()` cache `decompose.mat` and the eigenpairs as `.rds` files
beside the input and reload them when present. The file names now carry the
variant (`_decompLsym` vs `_decompLrwLegacy`), because a cache written before this
change would otherwise be reloaded after it and silently reinstate the bug.

### Variable subsets were cut on unstable ties

`createSubInputData.spectral()` ranked eigenvector elements with
`sort(..., method = "quick")`, which is **not stable**, so tied elements were
ordered unpredictably and the `G+,m` / `G-,m` cut could fall differently on the
same data. It now orders by value and then by variable name.

### Hill climbing followed the column order

`bnlearn::hc()` keeps the first best-scoring move it finds, so the order in which
it visits candidates decides every tie and steers the greedy path — and it visits
them in the data frame's column order, which was the input column order. Measured
on `BRUVS_U1`: permutation Jaccard **0.74–0.88**, renaming **0.90–0.96**.

`rMethod.hc()` now does two things. It sorts the variables into a canonical order
derived from **their own data** (the value vector, ties broken by name), which is
invariant to renaming — the data does not move — and to column permutation — the
sorted sequence is unchanged. Names and column positions are equally arbitrary
labels, so ordering by either is invariant to one input change and not the other;
only a data-derived order escapes both. It then repeats the search under
`argHcOrders` random orderings and keeps the arcs that recur in at least
`argHcOrderFreq` of them, so the result does not rest on one arbitrary scan order
either. Both renaming and permutation now give Jaccard **1.0000**.

Arcs are read with `bnlearn::arcs()` instead of walking the `graphNEL` by
position, which also fixes a latent mismatch: the old loop indexed node names
positionally against the *unfiltered* column list, so a single-level variable
dropped before the search could shift every subsequent name.

New arguments: `scalenet(argHcOrders = 24, argHcOrderFreq = 0.5)`. Hill-climbing
restarts and ordering consensus both diversify the greedy search, so under
consensus the restarts are redundant; set `restart = 0` in `argReconsParam` when
you turn orderings on. Unlike the Rust port's CLI, this function does not override
a restart count you passed.

### Also fixed

`createSubOutputGraphs()` read the reconstruction parameters flat
(`recons.param[["estimator"]]`) on the `spectral` path but nested
(`recons.param[["aracne"]][["estimator"]]`) on the clustering paths — and
`scalenet()`'s own documented example uses the nested form, so the documented
spectral + ARACNE call passed an empty estimator to `minet` and died with
"unknown estimator". Both forms are now accepted.

### Cross-validation against the Rust port

For both the legacy and the fixed path, this package and
[xscalenet](https://github.com/gmtsciencedev/xscalenet) produce **identical edge
sets** on `BRUVS_own` and `BRUVS_U1` (95/95 and 77/77 legacy, 203/203 and 74/74
fixed) — two independent implementations agreeing edge for edge, before and after.
