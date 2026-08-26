#!/usr/bin/env python3
"""Exact fixed-frame projections of the two generalized metrics.

The doubled-coordinate order is
    (tilde x_+, tilde x_-, tilde y; x^+, x^-, y),
and the frame columns are
    (oplus, ominus, y), (bar oplus, bar ominus, bar y).
All calculations are exact SymPy algebra.
"""

from __future__ import annotations

import sympy as sp


u = sp.symbols("u", positive=True)
lp, lm = sp.symbols("L_+ L_-", positive=True)
W = sp.symbols("W", real=True)
Pi = lp * lm
F = 1 / u + Pi * u

J = sp.zeros(6)
J[:3, 3:] = sp.eye(3)
J[3:, :3] = sp.eye(3)

Hinf = sp.Matrix(
    [
        [0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, -1, 0],
        [0, 0, 1, 0, 0, 0],
        [1, 0, 0, 0, 0, 0],
        [0, -1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1],
    ]
)

rt2 = sp.sqrt(2)
Vlow = sp.Matrix(
    [
        [1 / rt2, 0, 0],
        [0, 0, 0],
        [0, 0, -1 / rt2],
        [0, -rt2, 0],
        [0, 0, 0],
        [0, 0, -1 / rt2],
    ]
)
Vbarlow = sp.Matrix(
    [
        [0, 0, 0],
        [0, 1 / rt2, 0],
        [0, 0, -1 / rt2],
        [0, 0, 0],
        [rt2, 0, 0],
        [0, 0, 1 / rt2],
    ]
)
Vup = J * Vlow
Vbarup = J * Vbarlow


def project(H: sp.Matrix) -> sp.Matrix:
    return sp.simplify(Vup.T * (H - Hinf) * Vbarup)


# Riemannian Banados generalized metric.
g = sp.Matrix([[2 * lp, -F, 0], [-F, 2 * lm, 0], [0, 0, 1]])
B = sp.Matrix([[0, -F, 0], [F, 0, 0], [0, 0, 0]])
ginv = sp.simplify(g.inv())
HR = sp.Matrix.vstack(
    sp.Matrix.hstack(ginv, -ginv * B),
    sp.Matrix.hstack(B * ginv, g - B * ginv * B),
)
proj_R = project(sp.simplify(HR)).applyfunc(sp.factor)


# Non-Riemannian exact generalized metric.  C and S are kept independent
# here; the exact solution later imposes C=cosh(chi), S=sinh(chi).
C, S, a = sp.symbols("C S a", nonzero=True)
HNR = sp.Matrix(
    [
        [0, 0, 0, C, -S / a, 0],
        [0, 0, 0, a * S, -C, 0],
        [0, 0, 1, 0, 0, 0],
        [C, a * S, 0, -W * a * S, W * C, 0],
        [-S / a, -C, 0, W * C, -W * S / a, 0],
        [0, 0, 0, 0, 0, 1],
    ]
)
proj_NR = project(HNR).applyfunc(sp.factor)


# Substitute the exact solution's variables:
# q = u*sqrt(Pi/2), chi=2*sqrt(2)*atanh(q), a=sqrt(lp/lm).
q = sp.symbols("q")
Cq = sp.cosh(2 * sp.sqrt(2) * sp.atanh(q))
Sq = sp.sinh(2 * sp.sqrt(2) * sp.atanh(q))
proj_NR_q = proj_NR.subs({C: Cq, S: Sq, a: sp.sqrt(lp / lm)})


def series_q(expr: sp.Expr, order: int = 4) -> sp.Expr:
    return sp.series(expr, q, 0, order).removeO().expand()


proj_NR_series = proj_NR_q.applyfunc(lambda z: series_q(z, 4))
proj_NR_u_series = proj_NR_series.subs(q, u * sp.sqrt(Pi / 2)).applyfunc(sp.simplify)


def assert_zero_matrix(M: sp.Matrix, label: str) -> None:
    if any(sp.simplify(x) != 0 for x in M):
        raise AssertionError(f"{label} failed:\n{sp.pretty(M)}")


# Expected leading falloffs in the fixed limiting frame.
expected_R_leading = sp.Matrix(
    [
        [2 * lp * u, 2 * Pi * u, 0],
        [2 * u, 2 * lm * u, 0],
        [0, 0, 0],
    ]
)
expected_NR_leading = sp.Matrix(
    [
        [2 * lp * u, W / 2, 0],
        [0, 2 * lm * u, 0],
        [0, 0, 0],
    ]
)

expected_R_exact = (
    u * (1 + Pi * u**2) / (1 - Pi * u**2) ** 2
) * sp.Matrix([[2 * lp, 2 * Pi, 0], [2, 2 * lm, 0], [0, 0, 0]])
expected_NR_exact = sp.Matrix(
    [[a * S, C * W / 2, 0], [0, S / a, 0], [0, 0, 0]]
)
assert_zero_matrix(proj_R - expected_R_exact, "R exact projection")
assert_zero_matrix(proj_NR - expected_NR_exact, "NR exact projection")

# The R matrix is exact.  Its leading series must match the EDFE modes.
proj_R_leading = proj_R.applyfunc(lambda z: sp.series(z, u, 0, 2).removeO())
print("R exact projected matrix:")
sp.print_latex(proj_R)
print("R leading projected matrix:")
sp.print_latex(proj_R_leading)
assert_zero_matrix(proj_R_leading - expected_R_leading, "R leading projection")

# The NR matrix has W in the hair channel.  When W=W0+W1*u+O(u^2),
# the leading state-data sector W0=0 gives W1*u in that entry.
W0, W1 = sp.symbols("W_0 W_1")
nr_leading = proj_NR_u_series.subs(W, W0 + W1 * u).applyfunc(
    lambda z: sp.series(z, u, 0, 2).removeO()
)
expected_NR_modes = expected_NR_leading.subs(W, W0 + W1 * u)
print("NR exact projected matrix (C=cosh chi, S=sinh chi):")
sp.print_latex(proj_NR)
print("NR projected matrix through O(u^3):")
sp.print_latex(proj_NR_u_series)
print("NR leading projected matrix:")
sp.print_latex(nr_leading)
assert_zero_matrix(nr_leading - expected_NR_modes, "NR leading projection")

# Negative controls: the R type-changing coefficient and the NR hair
# coefficient must not be dropped or sign-flipped.
if sp.simplify(proj_R_leading[1, 0]) == 0:
    raise AssertionError("R negative control failed: missing type-changing mode")
if sp.simplify(nr_leading[0, 1] - (-W0 - W1 * u) / 2) == 0:
    raise AssertionError("NR negative control failed: hair sign was not fixed")
if sp.simplify(nr_leading[1, 0]) != 0:
    raise AssertionError("NR negative control failed: spurious type-changing mode")

print("[PASS] R and NR fixed-frame projections match the linearized EDFE modes")
