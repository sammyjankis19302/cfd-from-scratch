# Original scripts

Everything here is my own hand-written work, kept unchanged as the record of
where each module started. Nothing has been tidied, renamed or corrected.

They are worth keeping for two reasons. A repository that shows growth is more
convincing than one that pretends the author was fluent from the start; and the
mistakes documented in the module notes above are only meaningful if the code
that contained them is still visible.

| file | what it is |
|---|---|
| `first_attempt_laplace.py` | The Jacobi Laplace solver analysed in `../README.md` |
| `2D_laplas_Jacobi.py` | Laplace with Jacobi iteration |
| `2D_Poisson.py` | Poisson with a non-zero source — the step the cavity needs |
| `102.py`, `103.py` | Sparse-matrix formulations (`scipy.sparse`), a direct alternative to iterating |

## Known issue, left unfixed

`2D_Poisson.py` sets `itt=5`. Five Jacobi sweeps cannot converge a 101×101
Poisson problem — it needs a few thousand, since Jacobi is `O(N²)` — so it
always prints `diverged heheheh`. Almost certainly `5000` was intended.

It is still the most important script in this folder: it is the first with a
**non-zero source term** (`o[20:31,20:31]=1`), which is the step from Laplace to
Poisson and exactly what the cavity's `∇²ψ = −ω` requires.
