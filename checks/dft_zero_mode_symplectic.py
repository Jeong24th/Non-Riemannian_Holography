"""Symbolic DFT zero-mode symplectic-current check.

The calculation uses the torsionless semi-covariant DFT connection and
symplectic potential of arXiv:1507.07545.  Indices are ordered as
(tilde +, tilde -, tilde y; +, -, y), with J = [[0,I],[I,0]].

We perturb the type-(1,1) non-Riemannian vacuum by the two zero-mode tangent
directions used in NR_Holography.tex:

    r = delta H^{+-},        w = delta H_{+-}.

Only y dependence is retained.  The script computes the bilinear radial
symplectic current Omega^y(r,w) directly from

    Theta^A = 2 (P^{AC}P^{BD} - Pbar^{AC}Pbar^{BD}) delta Gamma_{BCD}.
"""

from __future__ import annotations

import sys

try:
    import sympy as sp
except ModuleNotFoundError:
    sys.path.insert(0, r"C:\tmp\sympy_packages")
    import sympy as sp


N = 6
D_PHYSICAL = 3
y, ell = sp.symbols("y ell", positive=True, finite=True)


def odd_metric() -> sp.Matrix:
    eye = sp.eye(3)
    zero = sp.zeros(3)
    return zero.row_join(eye).col_join(eye.row_join(zero))


J = odd_metric()


def nr_vacuum_metric() -> sp.Matrix:
    h = sp.zeros(N)
    h[0, 3] = h[3, 0] = 1
    h[1, 4] = h[4, 1] = -1
    h[2, 2] = 1
    h[5, 5] = 1
    return h


H0 = nr_vacuum_metric()


def partial(expr: sp.Expr, index: int) -> sp.Expr:
    """Physical-section derivative; the reduced problem depends only on y."""
    if index == 5:
        return sp.diff(expr, y)
    return sp.S.Zero


def matrix_partial(matrix: sp.Matrix, index: int) -> sp.Matrix:
    return matrix.applyfunc(lambda value: partial(value, index))


def connection(H_down: sp.Matrix, d: sp.Expr) -> list[list[list[sp.Expr]]]:
    """Return Gamma[C][A][B] with every index lowered.

    This implements Eq. (2.19) of arXiv:1507.07545 (v2), using
    antisymmetrization with weight one half.
    """
    P_down = (J + H_down) / 2
    Pb_down = (J - H_down) / 2
    P_mixed = P_down * J
    Pb_mixed = Pb_down * J
    P_right = J * P_down
    Pb_right = J * Pb_down

    dP = [matrix_partial(P_down, c) for c in range(N)]

    # T[E,D] = (P partial^E P Pbar)_{E D}; partial^E=J^{EF}partial_F.
    raised_dP: list[sp.Matrix] = []
    for e in range(N):
        raised_dP.append(sum((J[e, f] * dP[f] for f in range(N)), sp.zeros(N)))

    t_ed = sp.zeros(N)
    for e in range(N):
        block = P_mixed * raised_dP[e] * Pb_right
        for dd in range(N):
            t_ed[e, dd] = block[e, dd]
    p_dp_pb_antisym = (t_ed - t_ed.T) / 2

    gamma = [[[sp.S.Zero for _ in range(N)] for _ in range(N)] for _ in range(N)]
    for c in range(N):
        first_block = P_mixed * dP[c] * Pb_right
        for a in range(N):
            for b in range(N):
                first = first_block[a, b] - first_block[b, a]

                second = sp.S.Zero
                for dd in range(N):
                    for e in range(N):
                        projector_factor = (
                            Pb_mixed[a, dd] * Pb_mixed[b, e]
                            - Pb_mixed[b, dd] * Pb_mixed[a, e]
                            - P_mixed[a, dd] * P_mixed[b, e]
                            + P_mixed[b, dd] * P_mixed[a, e]
                        )
                        # Overall 2 times weight-1/2 antisymmetrization.
                        second += projector_factor * dP[dd][e, c]

                third = sp.S.Zero
                for dd in range(N):
                    projected_ca_b = (
                        Pb_down[c, a] * Pb_mixed[b, dd]
                        - Pb_down[c, b] * Pb_mixed[a, dd]
                        + P_down[c, a] * P_mixed[b, dd]
                        - P_down[c, b] * P_mixed[a, dd]
                    ) / 2
                    # Einstein summation over E in
                    # (P partial^E P Pbar)_[E D].
                    dilaton_vector = partial(d, dd) + sum(
                        p_dp_pb_antisym[ee, dd] for ee in range(N)
                    )
                    third += projected_ca_b * dilaton_vector

                gamma[c][a][b] = sp.simplify(first + second - 4 * third / (D_PHYSICAL - 1))
    return gamma


def theta(H_background: sp.Matrix, d_background: sp.Expr, variation: sp.Matrix) -> list[sp.Expr]:
    """Compute Theta^A for a tangent variation of H at fixed d."""
    tau = sp.symbols("tau")
    gamma_tau = connection(H_background + tau * variation, d_background)
    delta_gamma = [
        [
            [sp.diff(gamma_tau[c][a][b], tau).subs(tau, 0) for b in range(N)]
            for a in range(N)
        ]
        for c in range(N)
    ]

    P_up = J * ((J + H_background) / 2) * J
    Pb_up = J * ((J - H_background) / 2) * J
    result: list[sp.Expr] = []
    for a in range(N):
        value = sp.S.Zero
        for b in range(N):
            for c in range(N):
                for dd in range(N):
                    value += 2 * (
                        P_up[a, c] * P_up[b, dd]
                        - Pb_up[a, c] * Pb_up[b, dd]
                    ) * delta_gamma[b][c][dd]
        result.append(sp.simplify(value))
    return result


def gamma_square_lagrangian(H_down: sp.Matrix, d: sp.Expr) -> sp.Expr:
    """The surface-improved (Gamma^2) DFT Lagrangian density.

    This is Eq. (3.31) of arXiv:1507.07545, before the overall
    gravitational factor 1/(16 pi G).
    """
    gamma = connection(H_down, d)
    P_up = J * ((J + H_down) / 2) * J
    Pb_up = J * ((J - H_down) / 2) * J

    value = sp.S.Zero
    for a in range(N):
        for b in range(N):
            for c in range(N):
                for dd in range(N):
                    projector = P_up[a, c] * P_up[b, dd] - Pb_up[a, c] * Pb_up[b, dd]
                    if projector == 0:
                        continue
                    contraction = sp.S.Zero
                    for e in range(N):
                        gamma_ac_up_e = sum(gamma[a][c][f] * J[f, e] for f in range(N))
                        gamma_ab_up_e = sum(gamma[a][b][f] * J[f, e] for f in range(N))
                        gamma_up_e_ab = sum(J[e, f] * gamma[f][a][b] for f in range(N))
                        contraction += (
                            gamma_ac_up_e * gamma[b][dd][e]
                            - gamma_ab_up_e * gamma[dd][c][e]
                            + sp.Rational(1, 2) * gamma_up_e_ab * gamma[e][c][dd]
                        )
                    value += projector * contraction
    return sp.simplify(sp.exp(-2 * d) * value)


def main() -> None:
    r = sp.Function("r")(y)
    w = sp.Function("w")(y)
    h_r = sp.zeros(N)
    h_w = sp.zeros(N)
    h_r[0, 1] = h_r[1, 0] = r
    h_w[3, 4] = h_w[4, 3] = w

    a, b = sp.symbols("a b")
    theta_w_on_r = theta(H0 + a * h_r, -y / ell, h_w)
    theta_r_on_w = theta(H0 + b * h_w, -y / ell, h_r)

    cross_wr = sp.diff(theta_w_on_r[5], a).subs(a, 0)
    cross_rw = sp.diff(theta_r_on_w[5], b).subs(b, 0)
    omega_y = sp.simplify(sp.exp(2 * y / ell) * (cross_wr - cross_rw))

    print("Theta_y cross term: delta_r Theta(h_w) =")
    print(sp.simplify(cross_wr))
    print("Theta_y cross term: delta_w Theta(h_r) =")
    print(sp.simplify(cross_rw))
    print("e^{-2d} Omega^y(r,w) =")
    print(omega_y)

    r0, r2, w0, w2 = sp.symbols("r0 r2 w0 w2")
    falloffs = {
        r: r0 + r2 * sp.exp(-2 * y / ell),
        w: w0 + w2 * sp.exp(-2 * y / ell),
        sp.diff(r, y): sp.diff(r0 + r2 * sp.exp(-2 * y / ell), y),
        sp.diff(w, y): sp.diff(w0 + w2 * sp.exp(-2 * y / ell), y),
        sp.diff(r, y, 2): sp.diff(r0 + r2 * sp.exp(-2 * y / ell), y, 2),
        sp.diff(w, y, 2): sp.diff(w0 + w2 * sp.exp(-2 * y / ell), y, 2),
    }
    print("On-shell falloff substitution =")
    print(sp.simplify(omega_y.xreplace(falloffs)))

    # Complete the two tangent directions through second order so that
    # H J H = J + O(a^2 b, a b^2, a^3, b^3).  With S=H J and
    # {S0,delta S}=0, the minimal completion is
    # S2=-1/2 S0 (delta S)^2.
    h1 = a * h_r + b * h_w
    h2 = -sp.Rational(1, 2) * H0 * J * h1 * J * h1
    h_completed = H0 + h1 + h2
    constraint = sp.expand(h_completed * J * h_completed - J)
    constraint_ab = constraint.applyfunc(
        lambda entry: sp.diff(sp.diff(entry, a), b).subs({a: 0, b: 0})
    )
    assert constraint_ab == sp.zeros(N)

    lhat = gamma_square_lagrangian(h_completed, -y / ell)
    lhat_cross = sp.simplify(sp.diff(sp.diff(lhat, a), b).subs({a: 0, b: 0}))
    print("Gamma^2 action density, coefficient of a*b =")
    print(lhat_cross)
    print("Gamma^2 density on shell =")
    print(sp.simplify(lhat_cross.xreplace(falloffs)))


if __name__ == "__main__":
    main()

