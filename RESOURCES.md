# Resources

The sources behind this repository, in case they help someone learning the same
material. The suggested route through the modules is in the
[README](README.md#start-here); this page is the reference list.

> **Nothing copyrighted is hosted here.** These are citations and links. Please
> obtain books and papers through your library or institutional access.

---

## Papers

| Paper | Used for |
|---|---|
| **Ghia, Ghia & Shin (1982)**, *J. Comput. Phys.* **48**, 387–411 · [doi:10.1016/0021-9991(82)90058-4](https://doi.org/10.1016/0021-9991(82)90058-4) | The lid-driven cavity benchmark. Tables I and II are the centreline velocities module 06 validates against, Table IV the wall vorticity, Table V the primary- and secondary-vortex properties. Their coupled strongly-implicit multigrid method is the deliberate contrast to this repository's Jacobi + explicit Euler. |
| **Nallasamy & Krishna Prasad (1977)**, *J. Fluid Mech.* **79**(2), 391–414 | Cavity flow to Re = 50 000; growth and decay of the secondary eddies. |
| **Tokuda (1968)**, *J. Fluid Mech.* **33**(4), 657–672 | Impulsive motion of a flat plate — the unsteady startup problem. |
| **Agarwal (2012)**, third-order accurate upwind scheme for Navier–Stokes at high Re | Directly relevant to module 06's central finding: first-order upwind's numerical viscosity is *why* higher-order upwinding exists. The natural next step from this work. |
| **Jeong & Hussain**, on the identification of a vortex (λ₂ criterion) | Objective vortex identification, for locating cavity vortex centres without eyeballing contours. |

Transcribed benchmark data from Ghia et al. lives in
[`06-lid-driven-cavity/reference/`](06-lid-driven-cavity/reference/) as CSV,
with the OCR caveats documented in the file headers.

---

## Software

Deliberately minimal — the point was to write the numerics, not to assemble them.

- **NumPy** and **Matplotlib** — the only two dependencies
- **Python 3** on Ubuntu
- **Git** — version control, learned while building this repository

---

<!-- TODO Uday — two sections still to add, whenever you have ten minutes:

## Books
A short honest list of what you ACTUALLY read, one line each on what it was
useful FOR. Three real entries beat fifteen aspirational ones; anyone in the
field can tell the difference.

## Courses and lectures
NPTEL, lecture series, university notes — same rule.

## Handwritten notes
If you scan the original derivations, link them here. theory/ holds the typed
versions; the handwritten originals are evidence the work was done by hand
first, which is unusual and worth showing. Scan to PDF, keep the files small.
-->
