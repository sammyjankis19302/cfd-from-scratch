"""Verification for the spectral analysis module.

The closed-form expressions are only useful if they really are what the
stencils do. Test 1 applies each stencil DIRECTLY to exp(i k x) and compares
against the formula -- no algebra trusted, everything checked.

  1. every closed form matches its stencil to machine precision
  2. central stencils have zero imaginary part (no dissipation, ever)
  3. BD-1 damps and FD-1 amplifies -- opposite signs, same magnitude
  4. all schemes are exact as K -> 0 (consistency)
  5. CD-2 group velocity goes negative beyond 4 points per wavelength

Run:
    python3 tests/test_spectral.py
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src.spectral import (FIRST, SECOND, group_velocity,
                          points_per_wavelength)

K = np.array([0.05, 0.3, 0.8, 1.5, 2.5, 3.0, np.pi])


def _symbol(coefs, offsets, order):
    """Apply the stencil to exp(i K j) directly, with dx = 1."""
    s = sum(c * np.exp(1j * K * o) for c, o in zip(coefs, offsets))
    return (-1j * s) / K if order == 1 else (-s) / K**2


STENCILS_1 = {
    "BD-1": ([1, -1], [0, -1]),
    "FD-1": ([1, -1], [1, 0]),
    "CD-2": ([0.5, -0.5], [1, -1]),
    "CD-4": ([-1/12, 8/12, -8/12, 1/12], [2, 1, -1, -2]),
}
STENCILS_2 = {
    "CD-2": ([1, -2, 1], [1, 0, -1]),
    "CD-4": ([-1/12, 16/12, -30/12, 16/12, -1/12], [2, 1, 0, -1, -2]),
    "FD-2": ([1, -2, 1], [2, 1, 0]),
    "BD-2": ([1, -2, 1], [0, -1, -2]),
}


def test_formulas_match_stencils():
    """The closed forms must BE the stencils, not merely resemble them."""
    print("  first derivative:")
    for name, (c, o) in STENCILS_1.items():
        d = np.abs(_symbol(c, o, 1) - FIRST[name](K)).max()
        print(f"    {name}  max difference = {d:.2e}")
        assert d < 1e-14, f"{name} formula does not match its stencil"

    print("  second derivative:")
    for name, (c, o) in STENCILS_2.items():
        d = np.abs(_symbol(c, o, 2) - SECOND[name](K)).max()
        print(f"    {name}  max difference = {d:.2e}")
        # 1e-13 not 1e-14: the second-derivative symbol divides by K^2, so at
        # the smallest K in the sample (0.05) round-off is amplified by 400.
        assert d < 1e-13, f"{name} formula does not match its stencil"

    print("  all eight closed forms verified against the stencils: PASS\n")


def test_central_stencils_have_no_dissipation():
    """Symmetric stencils are purely real: they never damp, at any wavenumber.

    Which is exactly why CD-2 is useless on its own for advection -- it has no
    mechanism at all to control the 2-point wave. Damping must then come from
    somewhere else (the time scheme, or added artificial viscosity).
    """
    for name in ("CD-2", "CD-4"):
        m = np.abs(FIRST[name](K).imag).max()
        print(f"  {name} max|Im(k_num/k)| = {m:.2e}")
        assert m < 1e-15, f"{name} should be purely real"
    print("  central stencils are non-dissipative: PASS\n")


def test_upwind_damps_downwind_amplifies():
    """BD-1 and FD-1 differ ONLY in the sign of the imaginary part.

    Same formal order, identical phase behaviour, opposite stability. The
    entire difference between a scheme that works and one that explodes lives
    in that sign.
    """
    bd, fd = FIRST["BD-1"](K), FIRST["FD-1"](K)
    print(f"  max|Re(BD-1) - Re(FD-1)| = {np.abs(bd.real - fd.real).max():.2e}")
    print(f"  max|Im(BD-1) + Im(FD-1)| = {np.abs(bd.imag + fd.imag).max():.2e}")
    print(f"  at K = pi:  Im(BD-1) = {bd.imag[-1]:+.4f}   Im(FD-1) = {fd.imag[-1]:+.4f}")
    assert np.abs(bd.real - fd.real).max() < 1e-15, "real parts should be identical"
    assert np.abs(bd.imag + fd.imag).max() < 1e-15, "imaginary parts should be opposite"
    assert bd.imag[-1] < 0 < fd.imag[-1], "BD-1 must damp, FD-1 must amplify"
    print("  identical phase, opposite dissipation: PASS\n")


def test_consistency_at_long_wavelength():
    """Every scheme is exact as K -> 0, AND approaches it at its formal order.

    The first version of this test asserted |k_num/k - 1| < 1e-7 at K = 1e-4
    and failed for BD-1, which gave 5e-5. That was the test being wrong, not
    the scheme: BD-1 is FIRST order, so its error is O(K) = 1e-4 -- exactly
    what it should be. An absolute tolerance cannot express "first order".

    So measure the RATE instead: halve K and see how far the error falls.
    That is the same idea as a grid convergence study, in wavenumber space,
    and it ties this module back to the Taylor analysis in module 01.
    """
    expected = {"BD-1": 1, "FD-1": 1, "CD-2": 2, "CD-4": 4}
    Ks = np.array([1e-2, 5e-3])

    print("  " + "stencil".ljust(9) + "err(K)".rjust(12) + "err(K/2)".rjust(12)
          + "measured p".rjust(13) + "formal".rjust(9))
    for name, p_formal in expected.items():
        e = np.abs(FIRST[name](Ks) - 1.0)
        p = np.log2(e[0] / e[1])
        print("  " + name.ljust(9) + f"{e[0]:12.3e}" + f"{e[1]:12.3e}"
              + f"{p:13.2f}" + f"{p_formal:9d}")
        assert abs(p - p_formal) < 0.05, f"{name}: measured p={p:.2f}"

    print("  each stencil approaches exactness at its formal order: PASS\n")


def test_negative_group_velocity():
    """CD-2 energy travels BACKWARDS below 4 points per wavelength.

    V_g/c = cos(K), which changes sign at K = pi/2, i.e. PPW = 4. Under-resolved
    modes carry energy upstream while the exact solution carries it downstream.
    This is not instability and not smearing -- it is signal going the wrong way,
    and it is invisible to a convergence study.
    """
    Kc = np.linspace(1e-6, np.pi, 200001)
    vg = group_velocity("CD-2", Kc)
    K_zero = Kc[vg < 0].min()
    ppw = 2 * np.pi / K_zero

    print(f"  V_g < 0 for K > {K_zero:.4f} = {K_zero/np.pi:.3f} pi")
    print(f"  i.e. fewer than {ppw:.1f} points per wavelength")
    print(f"  at K = pi:  V_g/c = {vg[-1]:+.4f}  (full speed, wrong direction)")
    assert abs(K_zero - np.pi / 2) < 1e-3, "sign change should be at K = pi/2"
    assert abs(ppw - 4.0) < 0.01, "threshold should be 4 PPW"
    print("  spurious upstream propagation confirmed: PASS\n")


def test_points_per_wavelength():
    """Resolving power, which is NOT the same as order of accuracy."""
    print("  " + "stencil".ljust(9) + "1% phase error".rjust(16)
          + "0.1%".rjust(10))
    res = {}
    for name in FIRST:
        a = points_per_wavelength(name, 0.01)
        b = points_per_wavelength(name, 0.001)
        res[name] = a
        print("  " + name.ljust(9) + f"{a:16.1f}" + f"{b:10.1f}")

    assert abs(res["BD-1"] - res["CD-2"]) < 0.1, \
        "BD-1 and CD-2 share a real part, so PPW must match"
    assert res["CD-4"] < res["CD-2"] / 2, "CD-4 should need far fewer points"
    print("\n  BD-1, FD-1 and CD-2 need the SAME 25.6 points per wavelength")
    print("  despite BD-1 being first order and CD-2 second order: PASS\n")


if __name__ == "__main__":
    tests = [
        ("closed forms match the stencils", test_formulas_match_stencils),
        ("central stencils do not dissipate", test_central_stencils_have_no_dissipation),
        ("upwind damps, downwind amplifies", test_upwind_damps_downwind_amplifies),
        ("consistency as K -> 0", test_consistency_at_long_wavelength),
        ("negative group velocity", test_negative_group_velocity),
        ("points per wavelength", test_points_per_wavelength),
    ]
    failed = 0
    for name, fn in tests:
        print(f"[{name}]")
        try:
            fn()
        except AssertionError as exc:
            print(f"  FAIL: {exc}\n")
            failed += 1
    if failed:
        print(f"{failed} test(s) failed.")
        sys.exit(1)
    print("All tests passed.")
