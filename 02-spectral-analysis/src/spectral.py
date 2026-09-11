"""Global spectral analysis of finite-difference stencils.

Order of accuracy tells you how error behaves as dx -> 0. It does NOT tell you
which scheme is more accurate on the grid you can actually afford. This module
answers that second question.

THE IDEA

Feed a single Fourier mode exp(i k x) into a stencil. The exact first
derivative returns i*k times the mode. A stencil returns i*k_num times it,
where k_num is the NUMERICAL WAVENUMBER -- generally complex, and generally
not equal to k.

Everything follows from comparing k_num/k against 1:

    Re(k_num/k) != 1   ->  the mode travels at the wrong SPEED   -> dispersion
    Im(k_num/k) != 0   ->  the mode grows or decays              -> dissipation

Module 01 established the baseline this is measured against: for the exact
advection equation every mode has |G| = 1 and travels at exactly c.

Write K = k*dx, the non-dimensional wavenumber. K ranges over (0, pi]; K = pi
is the shortest representable wave, two points per wavelength. Points per
wavelength is PPW = 2*pi/K.
"""

import numpy as np


# ---------------------------------------------------------------------------
# First derivative: k_num / k
# ---------------------------------------------------------------------------
def bd1_first(K):
    """Backward difference. Real part = sin(K)/K, imaginary part negative."""
    return np.sin(K) / K - 1j * (1 - np.cos(K)) / K


def fd1_first(K):
    """Forward difference. Same real part as BD-1, imaginary part POSITIVE.

    That sign is the whole story: for c > 0 the negative imaginary part of BD-1
    damps, and the positive one of FD-1 amplifies. Same formal order, same
    phase behaviour, opposite stability.
    """
    return np.sin(K) / K + 1j * (1 - np.cos(K)) / K


def cd2_first(K):
    """Second-order central. Purely real -- NO dissipation at any wavenumber."""
    return np.sin(K) / K + 0j


def cd4_first(K):
    """Fourth-order central. Also purely real."""
    return (8 * np.sin(K) - np.sin(2 * K)) / (6 * K) + 0j


# ---------------------------------------------------------------------------
# Second derivative: k_num^2 / k^2
# ---------------------------------------------------------------------------
def cd2_second(K):
    """Second-order central second derivative. Real: (sin(K/2)/(K/2))^2."""
    return (np.sin(K / 2) / (K / 2)) ** 2 + 0j


def cd4_second(K):
    """Fourth-order central second derivative."""
    return (30 - 32 * np.cos(K) + 2 * np.cos(2 * K)) / (12 * K**2) + 0j


def fd2_second(K):
    """One-sided forward second derivative: u[i], u[i+1], u[i+2].

    The one-sided stencils acquire an imaginary part that the symmetric ones
    do not -- the same asymmetry that makes them dissipative for the first
    derivative.
    """
    s = (np.sin(K / 2) / (K / 2)) ** 2
    return s * np.cos(K) + 1j * s * np.sin(K)


def bd2_second(K):
    """One-sided backward second derivative: u[i], u[i-1], u[i-2]."""
    s = (np.sin(K / 2) / (K / 2)) ** 2
    return s * np.cos(K) - 1j * s * np.sin(K)


FIRST = {"BD-1": bd1_first, "FD-1": fd1_first, "CD-2": cd2_first, "CD-4": cd4_first}
SECOND = {"CD-2": cd2_second, "CD-4": cd4_second,
          "FD-2": fd2_second, "BD-2": bd2_second}


# ---------------------------------------------------------------------------
# Derived quantities
# ---------------------------------------------------------------------------
def group_velocity(name, K):
    """V_g / c = d(Re(k_num) dx) / d(K).

    The PHASE speed says how fast an individual crest moves. The GROUP velocity
    says how fast energy -- a wave packet's envelope -- moves. For the exact
    equation both equal c for every wavenumber.

    Where V_g < 0, energy at that wavenumber travels UPSTREAM while the exact
    solution carries it downstream. Those are the spurious upstream q-waves:
    they are not instability (the amplitude may be perfectly bounded) and not
    smearing. They are signal going the wrong way.
    """
    analytic = {
        "CD-2": lambda K: np.cos(K),
        "BD-1": lambda K: np.cos(K),
        "FD-1": lambda K: np.cos(K),
        "CD-4": lambda K: (8 * np.cos(K) - 2 * np.cos(2 * K)) / 6,
    }
    return analytic[name](K)


def points_per_wavelength(name, tol=0.01, n=200001):
    """Smallest PPW at which the phase error stays below tol.

    This is the practically useful number: how finely must I resolve a feature
    before this scheme represents it faithfully? It is a statement about the
    grid you can afford, not about the limit dx -> 0.
    """
    K = np.linspace(1e-9, np.pi, n)
    err = np.abs(FIRST[name](K).real - 1.0)
    bad = err > tol
    K_max = K[bad].min() if bad.any() else np.pi
    return 2 * np.pi / K_max
