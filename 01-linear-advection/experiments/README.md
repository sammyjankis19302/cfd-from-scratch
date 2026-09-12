# Original scripts

Everything here is my own hand-written work, kept unchanged as the record of
where each module started. Nothing has been tidied, renamed or corrected.

They are worth keeping for two reasons. A repository that shows growth is more
convincing than one that pretends the author was fluent from the start; and the
mistakes documented in the module notes above are only meaningful if the code
that contained them is still visible.

| file | what it is |
|---|---|
| `first_attempt.py` | The very first 1D advection script — the one with the frozen inflow boundary documented in `../BUGS.md` |
| `1D_advection_FTBS.py` | First-order upwind |
| `1D_advection_FTCS.py` | Centred space — unstable, kept deliberately |
| `1D_advection_FTFS.py` | Downwind — unstable, kept deliberately |
| `r1.py`, `r2.py`, `experiment_102.py` | Scratch exploration files |
| `101.py`, `104.py`, `105.py` | Scheme comparisons, including the FTCS instability |
