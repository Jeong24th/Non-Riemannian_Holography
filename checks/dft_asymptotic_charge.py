"""Leading modified DFT Noether potential for a chiral epsilon^+(x^+).

The Park--Rey--Rim--Sakatani global-charge formula strictly assumes
partial_B partial_[C X_D]=0.  A Virasoro parameter does not obey that
restriction, so this script is diagnostic rather than a completed
covariant-phase-space charge construction.  It determines which state
data occur in the leading surface integrand before the missing AdS
Wald/counterterm improvement is supplied.
"""

from __future__ import annotations

import sys

try:
    import sympy as sp
except ModuleNotFoundError:
    sys.path.insert(0, r"C:\tmp\sympy_packages")
    import sympy as sp


N = 6
xp, xm, y, ell = sp.symbols("xplus xminus y ell", positive=True)
u = sp.exp(-2 * y / ell)
lp = sp.Function("Lplus")(xp)
lm = sp.Function("Lminus")(xm)
w1 = sp.Function("W1")(xp, xm)
eps = sp.Function("eps")(xp)


def odd_metric() -> sp.Matrix:
    eye = sp.eye(3)
    zero = sp.zeros(3)
    return zero.row_join(eye).col_join(eye.row_join(zero))


J = odd_metric()


def partial(expr: sp.Expr, index: int) -> sp.Expr:
    if index == 3:
        return sp.diff(expr, xp)
    if index == 4:
        return sp.diff(expr, xm)
    if index == 5:
        return sp.diff(expr, y)
    return sp.S.Zero


def leading_metric() -> sp.Matrix:
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
    h_up = J * H_down * J
    result = sp.zeros(N, 1)
    for a in range(N):
        result[a] = sp.simplify(
            4 * sum(h_up[a, b] * partial(d, b) for b in range(N))
            - sum(partial(h_up[a, b], b) for b in range(N))
        )
    return result


def noether_potential(H_down: sp.Matrix, X_up: sp.Matrix) -> sp.Matrix:
    h_up = J * H_down * J
    h_first_up = J * H_down
    h_second_up = H_down * J
    x_down = J * X_up
    result = sp.zeros(N)
    for a in range(N):
        for b in range(N):
            value = sp.S.Zero
            for c in range(N):
                value -= h_up[c, a] * (
                    partial(X_up[b], c)
                    + sum(J[b, f] * partial(x_down[c], f) for f in range(N))
                )
                value += h_up[c, b] * (
                    partial(X_up[a], c)
                    + sum(J[a, f] * partial(x_down[c], f) for f in range(N))
                )
                for dd in range(N):
                    for e in range(N):
                        value -= (
                            h_up[c, a] * h_up[b, dd]
                            - h_up[c, b] * h_up[a, dd]
                        ) * partial(H_down[dd, e], c) * X_up[e]
                        value -= sp.Rational(1, 2) * h_second_up[e, c] * (
                            h_first_up[a, dd] * partial(h_up[b, dd], c)
                            - h_first_up[b, dd] * partial(h_up[a, dd], c)
                        ) * X_up[e]
            for e in range(N):
                value += sum(
                    J[a, f] * partial(h_first_up[b, e], f)
                    - J[b, f] * partial(h_first_up[a, e], f)
                    for f in range(N)
                ) * X_up[e]
            result[a, b] = sp.simplify(value)
    return result


def main() -> None:
    h = leading_metric()
    d = -y / ell
    x = sp.Matrix(
        [
            0,
            ell**2 * u * lm * sp.diff(eps, xp, 2) / 2,
            -ell * sp.diff(eps, xp) / 2,
            eps,
            0,
            -ell * sp.diff(eps, xp) / 2,
        ]
    )
    k = noether_potential(h, x)
    bvec = boundary_vector(h, d)
    khat = sp.zeros(N)
    for a in range(N):
        for b in range(N):
            khat[a, b] = sp.simplify(k[a, b] + x[a] * bvec[b] - x[b] * bvec[a])

    # The plus translation appeared in the (-,y) physical component in
    # the constant-parameter check.
    expr = sp.expand(sp.exp(2 * y / ell) * khat[4, 5])
    print("e^(-2d) Khat^(-,y) for epsilon^+(x^+) =")
    print(expr)
    print("boundary limit =")
    print(sp.limit(expr, y, sp.oo))


if __name__ == "__main__":
    main()

