# Original scripts

Everything here is my own hand-written work, kept unchanged as the record of
where each module started. Nothing has been tidied, renamed or corrected.

They are worth keeping for two reasons. A repository that shows growth is more
convincing than one that pretends the author was fluent from the start; and the
mistakes documented in the module notes above are only meaningful if the code
that contained them is still visible.

| file | what it is |
|---|---|
| `first_attempt_cavity.py` | `2D_Ghia8.py` — the final version, analysed in `../README.md` |
| `2D_Ghia.py` … `2D_Ghia7.py` | The parameter study: grids 51/101/201, Re 100/400/1000, t* 50/100/200. Same solver throughout; only the parameters differ. |
| `ghiii.py`, `t12.py` | Earlier cavity attempts |

`../figures/original_runs/` holds the streamfunction, vorticity and velocity
plots those runs produced. They are the original output, kept alongside the
regenerated figures so the two can be compared.

## Known issue, left unfixed

`2D_Ghia.py` **does not parse**: its first line reads `mport numpy as np`. One
missing character. The later versions in the series are unaffected.
