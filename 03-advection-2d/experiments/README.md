# Original scripts

Everything here is my own hand-written work, kept unchanged as the record of
where each module started. Nothing has been tidied, renamed or corrected.

They are worth keeping for two reasons. A repository that shows growth is more
convincing than one that pretends the author was fluent from the start; and the
mistakes documented in the module notes above are only meaningful if the code
that contained them is still visible.

| file | what it is |
|---|---|
| `2D_advection_linear.py` | 2D linear advection |
| `t1_2D_Advection.py` | Earlier 2D advection attempt |

## Known issues, left unfixed

`t1_2D_Advection.py` **does not parse**: line 55 reads `for j in range(1, ny):i`,
and that stray `i` turns the loop into a one-liner, making the nested `for k`
block an IndentationError. `2D_advection_linear.py` is the working version.

Both use `range(1, ny)` and `range(1, nx)`, so the `x = 0` and `y = 0` inflow
boundaries are never updated — the same frozen-inflow pattern documented in
`../../01-linear-advection/BUGS.md`. It happens to be harmless here, because
`sin(πx)sin(πy)` is exactly zero along both those edges and stays zero. As in
1D, the code is right for a reason that has nothing to do with the loop bounds.
