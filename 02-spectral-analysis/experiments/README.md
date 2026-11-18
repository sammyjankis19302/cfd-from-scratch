# Original scripts

Everything here is my own hand-written work, kept unchanged as the record of
where each module started. Nothing has been tidied, renamed or corrected.

They are worth keeping for two reasons. A repository that shows growth is more
convincing than one that pretends the author was fluent from the start; and the
mistakes documented in the module notes above are only meaningful if the code
that contained them is still visible.

| file | what it is |
|---|---|
| `first_attempt_order1.py` | First-derivative symbols — `k_num/k` for FD-1, BD-1, CD-2, CD-4, plotted as real and imaginary parts |
| `first_attempt_order2.py` | Second-derivative symbols — `k_num²/k²` for CD-2, CD-4, FD-2, BD-2 |
| `make_figures.py` | Regenerates this module's `figures/` from `../src/spectral.py` |

## A note on these two

Unlike the other modules, there is no bug to document here. Every one of the
eight closed-form expressions in these two scripts was checked by applying the
stencil **directly** to `exp(ikx)` and comparing — and all eight agreed to
`1e-14` or better, first time.

That is why `src/spectral.py` is a tidied version of this work rather than a
correction of it. The derivations were done on paper before either script was
written, which is the whole reason they came out right.

`tests/test_spectral.py` performs that stencil-versus-formula check on every
run, so the agreement is re-established every time the suite is executed rather
than asserted once here.
