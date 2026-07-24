"""Leading DFT Noether potential for constant boundary translations.

This evaluates the state-dependent part of the modified potential

    Khat^{AB}=K^{AB}+2 X^{[A} B^{B]}

using Eqs. (3.34), (A.4) of arXiv:1507.07545.  Constant translation
parameters obey the restriction needed by that global-charge formula,
so this check is narrower but cleaner than applying it directly to the
full Virasoro parameters.
"""

from __future__ import annotations

import sys

try:
    import sympy as sp
except ModuleNotFoundError:
    sys.path.insert(0, r"C:\tmp\sympy_packages")
    import sympy as sp


N = 6
y, ell = sp.symbols("y ell", positive=True)
lp, lm, w1 = sp.symbols("Lplus Lminus W1")
ep, em = sp.symbols("eps_plus eps_minus")


def odd_metric() -> sp.Matrix:
    eye = sp.eye(3)
    zero = sp.zeros(3)
    return zero.row_join(eye).col_join(eye.row_join(zero))


J = odd_metric()


def partial(expr: sp.Expr, index: int) -> sp.Expr:
    return sp.diff(expr, y) if index == 5 else sp.S.Zero


def leading_metric() -> sp.Matrix:
    u = sp.exp(-2 * y / ell)
    h = sp.zeros(N)
    h[0, 3] = h[3, 0] = 1
    h[1, 4] = h[4, 1] = -1
    h[2, 2] = 1
    h[5, 5] = 1
    h[0, 4] = h[4, 0] = -2 * lm * u
    h[1, 3] = h[3, 1] = 2 * lp * u
    h[3, 4] = h[4, 3] = w1 * u
    return h


def boundary_vector(H_down: sp.Matrix, d: sp.Expr) -> sp.Matrix:
    H_up = J * H_down * J
    result = sp.zeros(N, 1)
    for a in range(N):
        result[a] = sp.simplify(
            4 * sum(H_up[a, b] * partial(d, b) for b in range(N))
            - sum(partial(H_up[a, b], b) for b in range(N))
        )
    return result


def noether_potential(H_down: sp.Matrix, X_up: sp.Matrix) -> sp.Matrix:
    H_up = J * H_down * J
    H_first_up = J * H_down
    H_second_up = H_down * J
    X_down = J * X_up
    result = sp.zeros(N)

    for a in range(N):
        for b in range(N):
            term1 = sp.S.Zero
            term2 = sp.S.Zero
            term3 = sp.S.Zero
            term4 = sp.S.Zero
            for c in range(N):
                term1 += -(
                    H_up[c, a]
                    * (
                        partial(X_up[b], c)
                        + sum(J[b, f] * partial(X_down[c], f) for f in range(N))
                    )
                    - H_up[c, b]
                    * (
                        partial(X_up[a], c)
                        + sum(J[a, f] * partial(X_down[c], f) for f in range(N))
                    )
                )
                for dd in range(N):
                    for e in range(N):
                        term2 += -(
                            H_up[c, a] * H_up[b, dd]
                            - H_up[c, b] * H_up[a, dd]
                        ) * partial(H_down[dd, e], c) * X_up[e]
                        term4 += -sp.Rational(1, 2) * H_second_up[e, c] * (
                            H_first_up[a, dd] * partial(H_up[b, dd], c)
                            - H_first_up[b, dd] * partial(H_up[a, dd], c)
                        ) * X_up[e]
            for e in range(N):
                term3 += sum(
                    J[a, f] * partial(H_first_up[b, e], f)
                    - J[b, f] * partial(H_first_up[a, e], f)
                    for f in range(N)
                ) * X_up[e]
            result[a, b] = sp.simplify(term1 + term2 + term3 + term4)
    return result


def main() -> None:
    h = leading_metric()
    u = sp.exp(-2 * y / ell)
    tangent_check = (h * J * h - J).applyfunc(
        lambda entry: sp.limit(sp.expand(entry) / u, y, sp.oo)
    )
    assert tangent_check == sp.zeros(N)
    d = -y / ell
    x = sp.Matrix([0, 0, 0, ep, em, 0])
    k = noether_potential(h, x)
    bvec = boundary_vector(h, d)
    khat = sp.zeros(N)
    for a in range(N):
        for bb in range(N):
            khat[a, bb] = sp.simplify(k[a, bb] + x[a] * bvec[bb] - x[bb] * bvec[a])

    for component in ((3, 5), (4, 5)):
        expr = sp.simplify(sp.exp(2 * y / ell) * khat[component[0], component[1]])
        print(f"e^(-2d) Khat^{component} =")
        print(expr)
        vacuum_parameter = ep if component == (3, 5) else em
        renormalized = expr + 4 * vacuum_parameter * sp.exp(2 * y / ell) / ell
        print("finite boundary coefficient after vacuum subtraction =")
        print(sp.limit(renormalized, y, sp.oo))


if __name__ == "__main__":
    main()

