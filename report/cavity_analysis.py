"""
Reproduces every cavity number quoted in Chapters 10-12 of the report from the saved
run folders.   Usage:  python3 cavity_analysis.py /path/to/Cavityyyy   (or "." from inside it)
Reads only saved outputs (parameters.py, run_info.txt, solution_fields.npz,
Grig study/grid_independence_summary.csv).  Needs numpy and scipy.
Benchmark: Ghia, Ghia & Shin (1982), Tables I, II and V.  Two printed values in Table I
cannot be right and are replaced, as in my compare_ghia.py, by -0.08664 (Re=3200, y=0.4531;
printed -0.86636) and -0.03111 (Re=10000, y=0.5; printed +0.03111); a check below shows that
excluding them instead changes no quoted number.  v(400, x=0.9063) (printed -0.23827) is
excluded from every comparison.  v(1000, x=0.9063) is -0.51500 in compare_ghia.py, a typing
slip; Ghia print -0.51550, which is used here.
Run g3 (Re = 100 on 257^2, "Grid study/g3") was finished on the internship computer; the runs
repeated there are identical, array for array, to the ones in this folder.
"""
import sys, os, glob, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
from scipy.optimize import brentq
ROOT = sys.argv[1] if len(sys.argv) > 1 else "Cavityyyy"

# ---- Ghia et al. (1982) centreline data (from compare_ghia.py) -------------
ghia_y = np.array(
    [
        1.0000,
        0.9766,
        0.9688,
        0.9609,
        0.9531,
        0.8516,
        0.7344,
        0.6172,
        0.5000,
        0.4531,
        0.2813,
        0.1719,
        0.1016,
        0.0703,
        0.0625,
        0.0547,
        0.0000,
    ]
)
ghia_u = {
    100: np.array(
        [
            1.00000,
            0.84123,
            0.78871,
            0.73722,
            0.68717,
            0.23151,
            0.00332,
            -0.13641,
            -0.20581,
            -0.21090,
            -0.15662,
            -0.10150,
            -0.06434,
            -0.04775,
            -0.04192,
            -0.03717,
            0.00000,
        ]
    ),
    400: np.array(
        [
            1.00000,
            0.75837,
            0.68439,
            0.61756,
            0.55892,
            0.29093,
            0.16256,
            0.02135,
            -0.11477,
            -0.17119,
            -0.32726,
            -0.24299,
            -0.14612,
            -0.10338,
            -0.09266,
            -0.08186,
            0.00000,
        ]
    ),
    1000: np.array(
        [
            1.00000,
            0.65928,
            0.57492,
            0.51117,
            0.46604,
            0.33304,
            0.18719,
            0.05702,
            -0.06080,
            -0.10648,
            -0.27805,
            -0.38289,
            -0.29730,
            -0.22220,
            -0.20196,
            -0.18109,
            0.00000,
        ]
    ),
    3200: np.array(
        [
            1.00000,
            0.53236,
            0.48296,
            0.46547,
            0.46101,
            0.34682,
            0.19791,
            0.07156,
            -0.04272,
            -0.08664,
            -0.24427,
            -0.34323,
            -0.41933,
            -0.37827,
            -0.35344,
            -0.32407,
            0.00000,
        ]
    ),
    5000: np.array(
        [
            1.00000,
            0.48223,
            0.46120,
            0.45992,
            0.46036,
            0.33556,
            0.20087,
            0.08183,
            -0.03039,
            -0.07404,
            -0.22855,
            -0.33050,
            -0.40435,
            -0.43643,
            -0.42901,
            -0.41165,
            0.00000,
        ]
    ),
    7500: np.array(
        [
            1.00000,
            0.47244,
            0.47048,
            0.47323,
            0.47167,
            0.34228,
            0.20591,
            0.08342,
            -0.03800,
            -0.07503,
            -0.23176,
            -0.32393,
            -0.38324,
            -0.43025,
            -0.43590,
            -0.43154,
            0.00000,
        ]
    ),
    10000: np.array(
        [
            1.00000,
            0.47221,
            0.47783,
            0.48070,
            0.47804,
            0.34635,
            0.20673,
            0.08344,
            -0.03111,
            -0.07540,
            -0.23186,
            -0.32709,
            -0.38000,
            -0.41657,
            -0.42537,
            -0.42735,
            0.00000,
        ]
    ),
}
ghia_x = np.array(
    [
        1.0000,
        0.9688,
        0.9609,
        0.9531,
        0.9453,
        0.9063,
        0.8594,
        0.8047,
        0.5000,
        0.2344,
        0.2266,
        0.1563,
        0.0938,
        0.0781,
        0.0703,
        0.0625,
        0.0000,
    ]
)
ghia_v = {
    100: np.array(
        [
            0.00000,
            -0.05906,
            -0.07391,
            -0.08864,
            -0.10313,
            -0.16914,
            -0.22445,
            -0.24533,
            0.05454,
            0.17527,
            0.17507,
            0.16077,
            0.12317,
            0.10890,
            0.10091,
            0.09233,
            0.00000,
        ]
    ),
    400: np.array(
        [
            0.00000,
            -0.12146,
            -0.15663,
            -0.19254,
            -0.22847,
            -0.38,
            -0.44993,
            -0.38598,
            0.05186,
            0.30174,
            0.30203,
            0.28124,
            0.22965,
            0.20920,
            0.19713,
            0.18360,
            0.00000,
        ]
    ),
    1000: np.array(
        [
            0.00000,
            -0.21388,
            -0.27669,
            -0.33714,
            -0.39188,
            -0.51550,
            -0.42665,
            -0.31966,
            0.02526,
            0.32235,
            0.33075,
            0.37095,
            0.32627,
            0.30353,
            0.29012,
            0.27485,
            0.00000,
        ]
    ),
    3200: np.array(
        [
            0.00000,
            -0.39017,
            -0.47425,
            -0.52357,
            -0.54053,
            -0.44307,
            -0.37401,
            -0.31184,
            0.00999,
            0.28188,
            0.29030,
            0.37119,
            0.42768,
            0.41906,
            0.40917,
            0.39560,
            0.00000,
        ]
    ),
    5000: np.array(
        [
            0.00000,
            -0.49774,
            -0.55069,
            -0.55408,
            -0.52876,
            -0.41442,
            -0.36214,
            -0.30018,
            0.00945,
            0.27280,
            0.28066,
            0.35368,
            0.42951,
            0.43648,
            0.43329,
            0.42447,
            0.00000,
        ]
    ),
    7500: np.array(
        [
            0.00000,
            -0.53858,
            -0.55216,
            -0.52347,
            -0.48590,
            -0.41050,
            -0.36213,
            -0.30448,
            0.00824,
            0.27348,
            0.28117,
            0.35060,
            0.41824,
            0.43564,
            0.44030,
            0.43979,
            0.00000,
        ]
    ),
    10000: np.array(
        [
            0.00000,
            -0.54302,
            -0.52987,
            -0.49099,
            -0.45863,
            -0.41496,
            -0.36737,
            -0.30719,
            0.00831,
            0.27224,
            0.28003,
            0.35070,
            0.41487,
            0.43124,
            0.43733,
            0.43983,
            0.00000,
        ]
    ),
}



ghia_v[400] = ghia_v[400].copy(); ghia_v[400][5] = -0.23827   # restore printed value; point is excluded below
PSI_GHIA = {100: -0.103423, 400: -0.113909, 1000: -0.117929, 3200: -0.120377, 5000: -0.118966,
            7500: -0.119976, 10000: -0.119731}   # Ghia et al. (1982), Table V

def load_runs():
    out = {}
    for pf in sorted(glob.glob(os.path.join(ROOT, "*", "*", "parameters.py"))):
        folder = os.path.dirname(pf); ns = {}; exec(open(pf).read(), ns)
        npz = glob.glob(folder + "/pictures/*/solution_fields.npz")
        if not npz: continue
        info = {}
        for line in open(glob.glob(folder + "/pictures/*/run_info.txt")[0]):
            if "=" in line: k, v = line.split("=", 1); info.setdefault(k.strip(), v.strip())
        z = np.load(npz[0])
        out[os.path.basename(folder)] = dict(Re=int(ns["Re"]), n=ns["nx"], dt=ns["dt"], info=info, **{k: z[k] for k in z.files})
    return out

def exact_psi(o, n):
    m, h = n - 2, 1 / (n - 1)
    T = sp.diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(m, m)); I = sp.identity(m)
    A = ((sp.kron(I, T) + sp.kron(T, I)) / h**2).tocsc()
    p = np.zeros((n, n)); p[1:-1, 1:-1] = spl.spsolve(A, -o[1:-1, 1:-1].ravel()).reshape(m, m); return p

R = load_runs()
CASES = ["P1","Q1","V1","g3","P2","Q2","V2","P3","Q3","V3","P4","Q4","V4","P5","Q5","V5","P6","Q6","P7","Q7","V6"]   # g3: Re=100 on 257^2 (internship PC)
print("nu_r = nu_num/nu = Re*h/2 (steady state: the time-step part of the numerical viscosity vanishes)")
print(f"{'run':4}{'Re':>6}{'N':>5}{'conv':>6}{'d_omega':>10}{'D':>7}{'nu_r':>7}{'psi_min':>11}{'vs Ghia%':>9}{'max|du|':>9}{'max|dv|':>9}{'area+%':>7}{'Poisson%':>9}")
for t in CASES:
    d = R[t]; Re, n = d["Re"], d["n"]; h = 1 / (n - 1)
    iu = np.arange(1, 16); iv = iu[iu != 5] if Re == 400 else iu
    du = np.interp(ghia_y, d["y"], d["u"][:, n // 2])[iu] - ghia_u[Re][iu]
    dv = np.interp(ghia_x, d["x"], d["v"][n // 2, :])[iv] - ghia_v[Re][iv]
    s = d["s"]; pe = exact_psi(d["o"], n)
    g = PSI_GHIA.get(Re); dev = f"{(s.min()-g)/g*100:9.2f}" if g else "      ---"
    print(f"{t:4}{Re:6d}{n:5d}{d['info']['Converged']:>6}{float(d['info']['Final error']):10.2e}"
          f"{(1/Re)*d['dt']*2/h**2:7.3f}{Re*h/2:7.2f}{s.min():11.6f}{dev}"
          f"{np.abs(du).max():9.4f}{np.abs(dv).max():9.4f}{(s[1:-1,1:-1]>1e-6).mean()*100:7.1f}"
          f"{abs(s.min()-pe.min())/abs(pe.min())*100:9.3f}")

import json
res={}
for t in CASES:
    d=R[t]; Re,n=d["Re"],d["n"]; h=1/(n-1); s=d["s"]
    res[t]=dict(Re=Re,n=n,dx=h,psi=float(s.min()),psi_g=PSI_GHIA.get(Re),area=float((s[1:-1,1:-1]>1e-6).mean()*100),ratio=Re*h/2)   # steady state: nu_num = U h / 2
json.dump(res,open("results.json","w"),indent=1)

print("\nRichardson extrapolation of psi_min (grids 51, 101, 129):")
for Re, tr in [(100, ["P1","Q1","V1"]), (400, ["P2","Q2","V2"]), (1000, ["P3","Q3","V3"])]:
    f = np.array([R[t]["s"].min() for t in tr]); hh = np.array([1/50, 1/100, 1/128])
    g = lambda p: (f[0]-f[1])/(hh[0]**p-hh[1]**p) - (f[1]-f[2])/(hh[1]**p-hh[2]**p)
    p = brentq(g, 0.05, 5); C = (f[0]-f[1])/(hh[0]**p-hh[1]**p)
    ps_p = f[0] - C*hh[0]**p; ps_1 = f[2] + (f[2]-f[1])/(hh[1]/hh[2]-1)
    print(f"  Re={Re:5d}: p={p:.2f}  psi*(p)={ps_p:.6f}  psi*(p=1)={ps_1:.6f}  Ghia={PSI_GHIA[Re]:.6f}")

print("\nRe = 100 with the fourth grid (257^2, run g3):")
f4 = {n: R[t]["s"].min() for n, t in [(51, "P1"), (101, "Q1"), (129, "V1"), (257, "g3")]}
for a, b, c in [(51, 101, 257), (51, 129, 257), (101, 129, 257)]:
    ha, hb, hc = 1/(a-1), 1/(b-1), 1/(c-1)
    g = lambda p: (f4[a]-f4[b])/(ha**p-hb**p) - (f4[b]-f4[c])/(hb**p-hc**p)
    p = brentq(g, 0.05, 6); C = (f4[a]-f4[b])/(ha**p-hb**p); ps = f4[a] - C*ha**p
    print(f"  grids {a}, {b}, {c}: p = {p:.2f}, psi* = {ps:.6f} ({(ps-PSI_GHIA[100])/PSI_GHIA[100]*100:+.2f} % from Ghia)")
ps2 = f4[257] + (f4[257]-f4[129])/((256/128)**2-1)
print(f"  129^2 and 257^2 with p = 2: psi* = {ps2:.6f} ({(ps2-PSI_GHIA[100])/PSI_GHIA[100]*100:+.2f} % from Ghia)")

# ---------------------------------------------------------------------------------------------
print("\nFirst-order extrapolation (p = 1) from the two finest grids at Re >= 5000:")
for Re, (a, b) in [(5000, ("Q5", "V5")), (10000, ("Q7", "V6"))]:
    fa, fb = R[a]["s"].min(), R[b]["s"].min(); ha, hb = 1/(R[a]["n"]-1), 1/(R[b]["n"]-1)
    ps = fb + (fb - fa)/(ha/hb - 1)
    print(f"  Re={Re:5d}: {R[a]['n']}^2 and {R[b]['n']}^2 -> psi* = {ps:.6f}, vs Ghia {(ps-PSI_GHIA[Re])/PSI_GHIA[Re]*100:+.1f} %")

print("\nPrimary-vortex centre (location of psi_min) and share of the cavity with psi > 0:")
for t in CASES:
    d = R[t]; s = d["s"]; j, i = np.unravel_index(np.argmin(s), s.shape)
    print(f"  {t:3} Re={d['Re']:5d} {d['n']:3d}^2: centre ({d['x'][i]:.3f}, {d['y'][j]:.3f})"
          f"   psi>0: {(s[1:-1,1:-1] > 1e-6).mean()*100:5.1f} %")

print("\nVorticity inside the vortex core (psi <= 0.7 psi_min) and at the lid, 129^2:")
for t in ["V1", "V2", "V3"]:
    d = R[t]; s, o = d["s"], d["o"]; n = d["n"]; core = s <= 0.7*s.min()
    top = o[-1, 1:-1]
    print(f"  Re={d['Re']:5d}: core mean {o[core].mean():.2f}, std {o[core].std():.2f};"
          f"  lid: mid-lid {o[-1, n//2]:.2f}, next to the top-left corner {top[0]:.1f}")
print("  lid vorticity next to the top-left / top-right corner at Re=100 on 51^2, 101^2, 129^2, 257^2:",
      ", ".join(f"{R[t]['o'][-1,1]:.0f} / {R[t]['o'][-1,-2]:.0f}" for t in ["P1", "Q1", "V1", "g3"]),
      "(grows with refinement, roughly as 2/h: corner singularity)")
d = R["g3"]; s_, o_ = d["s"], d["o"]; core = s_ <= 0.7*s_.min()
print(f"  Re=100 on 257^2: core mean {o_[core].mean():.2f}, std {o_[core].std():.2f}; mid-lid {o_[-1, d['n']//2]:.2f}")

print("\nGrid study, 51^2 against 101^2 (my grid_independence_summary.csv) and psi_min on both grids:")
gs = np.loadtxt(os.path.join(ROOT, "Grig study", "grid_independence_summary.csv"), delimiter=",", skiprows=1)
for row, (p, q) in zip(gs, [("P1","Q1"),("P2","Q2"),("P3","Q3"),("P4","Q4"),("P5","Q5"),("P6","Q6"),("P7","Q7")]):
    a, b = R[p]["s"].min(), R[q]["s"].min()
    print(f"  Re={int(row[0]):5d}: max|du| = {row[2]:.3f}, max|dv| = {row[4]:.3f};"
          f"  psi_min {a:.4f} vs {b:.4f} (differ by {abs(a-b)/abs(b)*100:.0f} %)")

print("\nExcluding the two corrected benchmark points instead of correcting them:")
for t in ["P4", "Q4", "V4", "P7", "Q7", "V6"]:
    d = R[t]; Re, n = d["Re"], d["n"]; iu = np.arange(1, 16)
    du = np.interp(ghia_y, d["y"], d["u"][:, n // 2]) - ghia_u[Re]
    k = int(np.argmin(abs(ghia_y - (0.4531 if Re == 3200 else 0.5))))
    print(f"  {t}: max|du| {np.abs(du[iu]).max():.4f} -> {np.abs(du[iu[iu != k]]).max():.4f}")

print("\nWall-clock time of each run (run_info.txt):")
for t in CASES:
    sec = R[t]["info"].get("Computational time", "")      # first entry is in seconds
    txt = f"{float(sec.split()[0])/3600:.1f} h" if sec else "n/a"
    print(f"  {t:3} Re={R[t]['Re']:5d} {R[t]['n']:3d}^2: {txt}")

# ---------------------------------------------------------------------------------------------
print("\nPoisson residual of the final fields, max |lap_h(psi) + omega| over the interior points")
print("(the check at the end of the final solver.py, applied to the saved fields):")
for t in CASES:
    d = R[t]; s_, o_ = d["s"], d["o"]; hx = d["x"][1] - d["x"][0]; hy = d["y"][1] - d["y"][0]
    lap = ((s_[1:-1, 2:] - 2*s_[1:-1, 1:-1] + s_[1:-1, :-2]) / hx**2
           + (s_[2:, 1:-1] - 2*s_[1:-1, 1:-1] + s_[:-2, 1:-1]) / hy**2)
    r = np.abs(lap + o_[1:-1, 1:-1]); j, i = np.unravel_index(np.argmax(r), r.shape)
    print(f"  {t:3} Re={d['Re']:5d} {d['n']:3d}^2: max residual {r.max():.2e} at (x, y) = "
          f"({d['x'][i+1]:.3f}, {d['y'][j+1]:.3f});  one sweep would change psi by h^2 r/4 = {r.max()*hx**2/4:.1e}")

