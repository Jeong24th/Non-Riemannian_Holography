"""Renormalized asymptotic DFT surface-charge one-form for the NR saddle.

This implements the Park--Rey--Rim--Sakatani covariant-phase-space
identity (arXiv:1507.07545, Sec. 3.2) for the source-free asymptotic
phase space of NR_Holography.tex.  The doubled index order is

    (tilde +, tilde -, tilde y; +, -, y).

For field-dependent asymptotic parameters the surface-charge density is

  k_X^{AB} = delta(e^{-2d} Khat_X^{AB})
             - e^{-2d} Khat_{delta X}^{AB}
             + 2 e^{-2d} X^{[A} Thetahat^{B]}.

The last two terms are absent from the earlier diagnostic calculation.
The generalized metric and dilaton are retained through z^2, where
z=exp(-2y/l); this is sufficient for every finite boundary term.
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
xp, xm, z, ell, tau = sp.symbols("xplus xminus z ell tau", positive=True)

Lp = sp.Function("Lplus")(xp)
Lm = sp.Function("Lminus")(xm)
W = sp.Function("W1")(xp, xm)
dLp = sp.Function("dLplus")(xp)
dLm = sp.Function("dLminus")(xm)
dW = sp.Function("dW1")(xp, xm)
ep = sp.Function("eps_plus")(xp)
em = sp.Function("eps_minus")(xm)
hp = sp.Function("eta_plus")(xp)
hm = sp.Function("eta_minus")(xm)


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
        return -2 * z * sp.diff(expr, z) / ell
    return sp.S.Zero


def matrix_partial(matrix: sp.Matrix, index: int) -> sp.Matrix:
    return matrix.applyfunc(lambda value: partial(value, index))


def nr_metric(lp: sp.Expr, lm: sp.Expr, w: sp.Expr) -> sp.Matrix:
    """Exact NR family expanded through z^2."""
    product = lp * lm
    h = sp.zeros(N)
    h[0, 3] = h[3, 0] = 1 + 2 * product * z**2
    h[1, 4] = h[4, 1] = -1 - 2 * product * z**2
    h[2, 2] = h[5, 5] = 1
    h[0, 4] = h[4, 0] = -2 * lm * z
    h[1, 3] = h[3, 1] = 2 * lp * z
    h[3, 3] = -2 * lp * w * z**2
    h[3, 4] = h[4, 3] = w * z
    h[4, 4] = -2 * lm * w * z**2
    return h


def nr_dilaton(lp: sp.Expr, lm: sp.Expr) -> sp.Expr:
    return sp.log(z) / 2 + lp * lm * z**2 / 4


def density(lp: sp.Expr, lm: sp.Expr) -> sp.Expr:
    return (1 - lp * lm * z**2 / 2) / z


def connection(h_down: sp.Matrix, d: sp.Expr) -> list[list[list[sp.Expr]]]:
    """Torsionless semi-covariant connection Gamma[C][A][B]."""
    p_down = (J + h_down) / 2
    pb_down = (J - h_down) / 2
    p_mixed = p_down * J
    pb_mixed = pb_down * J
    p_right = J * p_down
    pb_right = J * pb_down
    dp = [matrix_partial(p_down, c) for c in range(N)]

    raised_dp = []
    for e in range(N):
        raised_dp.append(sum((J[e, f] * dp[f] for f in range(N)), sp.zeros(N)))

    t_ed = sp.zeros(N)
    for e in range(N):
        block = p_mixed * raised_dp[e] * pb_right
        for dd in range(N):
            t_ed[e, dd] = block[e, dd]
    antisym = (t_ed - t_ed.T) / 2

    gamma = [[[sp.S.Zero for _ in range(N)] for _ in range(N)] for _ in range(N)]
    for c in range(N):
        first_block = p_mixed * dp[c] * pb_right
        for a in range(N):
            for b in range(N):
                first = first_block[a, b] - first_block[b, a]
                second = sp.S.Zero
                for dd in range(N):
                    for e in range(N):
                        projector = (
                            pb_mixed[a, dd] * pb_mixed[b, e]
                            - pb_mixed[b, dd] * pb_mixed[a, e]
                            - p_mixed[a, dd] * p_mixed[b, e]
                            + p_mixed[b, dd] * p_mixed[a, e]
                        )
                        second += projector * dp[dd][e, c]

                third = sp.S.Zero
                for dd in range(N):
                    projected = (
                        pb_down[c, a] * pb_mixed[b, dd]
                        - pb_down[c, b] * pb_mixed[a, dd]
                        + p_down[c, a] * p_mixed[b, dd]
                        - p_down[c, b] * p_mixed[a, dd]
                    ) / 2
                    dilaton_vector = partial(d, dd) + sum(
                        antisym[ee, dd] for ee in range(N)
                    )
                    third += projected * dilaton_vector
                gamma[c][a][b] = sp.expand(first + second - 4 * third / (D_PHYSICAL - 1))
    return gamma


def boundary_vector(h_down: sp.Matrix, d: sp.Expr) -> sp.Matrix:
    h_up = J * h_down * J
    result = sp.zeros(N, 1)
    for a in range(N):
        result[a] = sp.expand(
            4 * sum(h_up[a, b] * partial(d, b) for b in range(N))
            - sum(partial(h_up[a, b], b) for b in range(N))
        )
    return result


def noether_component(h_down: sp.Matrix, x_up: sp.Matrix, a: int, b: int) -> sp.Expr:
    """Component K^{ab}[X], Eq. (A.4) of arXiv:1507.07545."""
    h_up = J * h_down * J
    h_first_up = J * h_down
    h_second_up = h_down * J
    x_down = J * x_up
    value = sp.S.Zero
    for c in range(N):
        value -= h_up[c, a] * (
            partial(x_up[b], c)
            + sum(J[b, f] * partial(x_down[c], f) for f in range(N))
        )
        value += h_up[c, b] * (
            partial(x_up[a], c)
            + sum(J[a, f] * partial(x_down[c], f) for f in range(N))
        )
        for dd in range(N):
            for e in range(N):
                value -= (
                    h_up[c, a] * h_up[b, dd]
                    - h_up[c, b] * h_up[a, dd]
                ) * partial(h_down[dd, e], c) * x_up[e]
                value -= sp.Rational(1, 2) * h_second_up[e, c] * (
                    h_first_up[a, dd] * partial(h_up[b, dd], c)
                    - h_first_up[b, dd] * partial(h_up[a, dd], c)
                ) * x_up[e]
    for e in range(N):
        value += sum(
            J[a, f] * partial(h_first_up[b, e], f)
            - J[b, f] * partial(h_first_up[a, e], f)
            for f in range(N)
        ) * x_up[e]
    return sp.expand(value)


def khat_component(
    h_down: sp.Matrix, d: sp.Expr, x_up: sp.Matrix, a: int, b: int
) -> sp.Expr:
    bvec = boundary_vector(h_down, d)
    return sp.expand(
        noether_component(h_down, x_up, a, b)
        + x_up[a] * bvec[b]
        - x_up[b] * bvec[a]
    )


def theta_hat_density_components(
    h: sp.Matrix,
    d: sp.Expr,
    dh: sp.Matrix,
    dd: sp.Expr,
    e_density: sp.Expr,
    components: tuple[int, ...],
) -> dict[int, sp.Expr]:
    """Return e^{-2d} Thetahat^A for a tangent phase-space variation."""
    gamma = connection(h, d)
    h_up = J * h * J
    dh_up = J * dh * J

    theta: dict[int, sp.Expr] = {}
    for a in components:
        value = 4 * sum(h_up[a, b] * partial(dd, b) for b in range(N))
        # - nabla_B delta H^{AB}; raised indices acquire minus connection terms.
        for b in range(N):
            value -= partial(dh_up[a, b], b)
            for c in range(N):
                for f in range(N):
                    value += gamma[b][c][f] * J[f, a] * dh_up[c, b]
                    value += gamma[b][c][f] * J[f, b] * dh_up[a, c]
        theta[a] = sp.expand(value)

    h_tau = nr_metric(Lp + tau * dLp, Lm + tau * dLm, W + tau * dW)
    d_tau = nr_dilaton(Lp + tau * dLp, Lm + tau * dLm)
    e_tau = density(Lp + tau * dLp, Lm + tau * dLm)
    b_tau = boundary_vector(h_tau, d_tau)
    result = {}
    for a in components:
        delta_e_b = sp.diff(e_tau * b_tau[a], tau).subs(tau, 0)
        result[a] = sp.expand(e_density * theta[a] - delta_e_b)
    return result


def x_plus_for(parameter: sp.Expr, lm: sp.Expr) -> sp.Matrix:
    return sp.Matrix(
        [
            0,
            ell**2 * z * lm * sp.diff(parameter, xp, 2) / 2,
            -ell * sp.diff(parameter, xp) / 2,
            parameter,
            0,
            -ell * sp.diff(parameter, xp) / 2,
        ]
    )


def x_minus_for(parameter: sp.Expr, lp: sp.Expr) -> sp.Matrix:
    return sp.Matrix(
        [
            -ell**2 * z * lp * sp.diff(parameter, xm, 2) / 2,
            0,
            ell * sp.diff(parameter, xm) / 2,
            0,
            parameter,
            -ell * sp.diff(parameter, xm) / 2,
        ]
    )


def x_plus(lm: sp.Expr) -> sp.Matrix:
    return x_plus_for(ep, lm)


def x_minus(lp: sp.Expr) -> sp.Matrix:
    return x_minus_for(em, lp)


def c_bracket(x: sp.Matrix, q: sp.Matrix) -> sp.Matrix:
    """DFT C-bracket of two doubled vectors."""
    x_down = J * x
    q_down = J * q
    result = sp.zeros(N, 1)
    for a in range(N):
        value = sum(
            x[b] * partial(q[a], b) - q[b] * partial(x[a], b)
            for b in range(N)
        )
        value += sp.Rational(1, 2) * sum(
            q_down[b] * sum(J[a, c] * partial(x[b], c) for c in range(N))
            - x_down[b] * sum(J[a, c] * partial(q[b], c) for c in range(N))
            for b in range(N)
        )
        result[a] = sp.expand(value)
    return result


def vector_is_zero(vector: sp.Matrix) -> bool:
    return all(sp.simplify(entry) == 0 for entry in vector)


def finite(expr: sp.Expr) -> sp.Expr:
    return sp.simplify(sp.limit(sp.expand(expr), z, 0, dir="+"))


def charge_one_form(
    chirality: str,
    h: sp.Matrix,
    d: sp.Expr,
    e_density: sp.Expr,
    theta_hat: dict[int, sp.Expr],
) -> sp.Expr:
    if chirality == "+":
        a, b = 4, 5
        x = x_plus(Lm)
        x_tau = x_plus(Lm + tau * dLm)
        theta_term = x[a] * theta_hat[b] - x[b] * theta_hat[a]
    else:
        a, b = 3, 5
        x = x_minus(Lp)
        x_tau = x_minus(Lp + tau * dLp)
        theta_term = x[a] * theta_hat[b] - x[b] * theta_hat[a]

    h_tau = nr_metric(Lp + tau * dLp, Lm + tau * dLm, W + tau * dW)
    d_tau = nr_dilaton(Lp + tau * dLp, Lm + tau * dLm)
    e_tau = density(Lp + tau * dLp, Lm + tau * dLm)
    varied = sp.diff(e_tau * khat_component(h_tau, d_tau, x_tau, a, b), tau).subs(
        tau, 0
    )
    delta_x = x_tau.diff(tau).subs(tau, 0)
    field_dependent = e_density * khat_component(h, d, delta_x, a, b)
    return finite(sp.expand(varied - field_dependent + theta_term))


def arbitrary_charge_one_form(
    x: sp.Matrix,
    a: int,
    b: int,
    h: sp.Matrix,
    d: sp.Expr,
    e_density: sp.Expr,
    theta_hat: dict[int, sp.Expr],
) -> sp.Expr:
    """Charge one-form for a field-independent doubled parameter."""
    h_tau = nr_metric(Lp + tau * dLp, Lm + tau * dLm, W + tau * dW)
    d_tau = nr_dilaton(Lp + tau * dLp, Lm + tau * dLm)
    e_tau = density(Lp + tau * dLp, Lm + tau * dLm)
    varied = sp.diff(e_tau * khat_component(h_tau, d_tau, x, a, b), tau).subs(
        tau, 0
    )
    theta_term = x[a] * theta_hat[b] - x[b] * theta_hat[a]
    return finite(sp.expand(varied + theta_term))


def main() -> None:
    h = nr_metric(Lp, Lm, W)
    d = nr_dilaton(Lp, Lm)
    e_density = density(Lp, Lm)

    # The truncation must obey H J H = J through the retained order.
    constraint = sp.expand(h * J * h - J)
    for entry in constraint:
        assert sp.expand(entry).coeff(z, 0) == 0
        assert sp.expand(entry).coeff(z, 1) == 0
        assert sp.expand(entry).coeff(z, 2) == 0

    h_tau = nr_metric(Lp + tau * dLp, Lm + tau * dLm, W + tau * dW)
    d_tau = nr_dilaton(Lp + tau * dLp, Lm + tau * dLm)
    dh = h_tau.diff(tau).subs(tau, 0)
    dd = sp.diff(d_tau, tau).subs(tau, 0)

    theta_hat = theta_hat_density_components(
        h, d, dh, dd, e_density, components=(3, 4, 5)
    )
    print("Finite e^{-2d} Thetahat components:")
    for a in (3, 4, 5):
        print(a, finite(theta_hat[a]))

    k_plus = charge_one_form("+", h, d, e_density, theta_hat)
    k_minus = charge_one_form("-", h, d, e_density, theta_hat)
    print("Finite k_plus^(minus,y) =")
    print(k_plus)
    print("Finite k_minus^(plus,y) =")
    print(k_minus)

    # Same-chirality brackets close without a field-dependent correction.
    alpha_plus = ep * sp.diff(hp, xp) - hp * sp.diff(ep, xp)
    alpha_minus = em * sp.diff(hm, xm) - hm * sp.diff(em, xm)
    plus_bracket = c_bracket(x_plus_for(ep, Lm), x_plus_for(hp, Lm))
    minus_bracket = c_bracket(x_minus_for(em, Lp), x_minus_for(hm, Lp))
    plus_difference = (plus_bracket - x_plus_for(alpha_plus, Lm)).applyfunc(sp.simplify)
    minus_difference = (minus_bracket - x_minus_for(alpha_minus, Lp)).applyfunc(sp.simplify)
    print("plus C-bracket closes:", vector_is_zero(plus_difference))
    print("plus bracket difference:", list(plus_difference))
    print("minus C-bracket closes:", vector_is_zero(minus_difference))
    print("minus bracket difference:", list(minus_difference))

    # The raw mixed C-bracket is cancelled by the Barnich--Troessaert
    # field-dependence terms: [X,Y]_*=[X,Y]_C-delta_X Y+delta_Y X.
    x_p = x_plus_for(ep, Lm)
    x_m = x_minus_for(em, Lp)
    delta_p_lp = ep * sp.diff(Lp, xp) + 2 * sp.diff(ep, xp) * Lp
    delta_m_lm = em * sp.diff(Lm, xm) + 2 * sp.diff(em, xm) * Lm
    delta_p_xm = x_minus_for(em, Lp + tau * delta_p_lp).diff(tau).subs(tau, 0)
    delta_m_xp = x_plus_for(ep, Lm + tau * delta_m_lm).diff(tau).subs(tau, 0)
    mixed_adjusted = c_bracket(x_p, x_m) - delta_p_xm + delta_m_xp
    mixed_adjusted = mixed_adjusted.applyfunc(sp.simplify)
    print("mixed adjusted C-bracket vanishes:", vector_is_zero(mixed_adjusted))
    print("mixed adjusted bracket:", list(mixed_adjusted))
    mixed_q_plus = finite(e_density * khat_component(h, d, mixed_adjusted, 4, 5))
    mixed_q_minus = finite(e_density * khat_component(h, d, mixed_adjusted, 3, 5))
    print("mixed residual plus potential:", mixed_q_plus)
    print("mixed residual minus potential:", mixed_q_minus)

    # The same-chirality bracket differs from the chosen representative by
    # a closed B-gauge parameter.  Test whether this reducibility parameter
    # has a nonzero surface charge (a possible winding central element).
    zeta_p = sp.Function("zeta_plus")(xp)
    zeta_m = sp.Function("zeta_minus")(xm)
    z_p = sp.Matrix([zeta_p, 0, 0, 0, 0, 0])
    z_m = sp.Matrix([0, zeta_m, 0, 0, 0, 0])
    k_zp = arbitrary_charge_one_form(z_p, 4, 5, h, d, e_density, theta_hat)
    k_zm = arbitrary_charge_one_form(z_m, 3, 5, h, d, e_density, theta_hat)
    q_zp = finite(e_density * khat_component(h, d, z_p, 4, 5))
    q_zm = finite(e_density * khat_component(h, d, z_m, 3, 5))
    print("closed B-gauge plus potential:", q_zp)
    print("closed B-gauge minus potential:", q_zm)
    print("closed B-gauge plus charge one-form:", k_zp)
    print("closed B-gauge minus charge one-form:", k_zm)

    plus_reducibility = all(plus_difference[i] == 0 for i in range(1, N))
    minus_reducibility = all(
        minus_difference[i] == 0 for i in range(N) if i != 1
    )
    assert plus_reducibility and minus_reducibility
    assert q_zp == 0 and q_zm == 0 and k_zp == 0 and k_zm == 0
    assert mixed_q_plus == 0 and mixed_q_minus == 0
    print("same-chirality charge brackets close modulo zero-charge reducibility: True")
    print("mixed charge bracket vanishes modulo zero-charge residuals: True")

    # The possible algebra cocycle is a total derivative on the boundary circle.
    delta_hp_lp = hp * sp.diff(Lp, xp) + 2 * sp.diff(hp, xp) * Lp
    plus_cocycle_density = sp.expand(ep * delta_hp_lp - alpha_plus * Lp)
    delta_hm_lm = hm * sp.diff(Lm, xm) + 2 * sp.diff(hm, xm) * Lm
    minus_cocycle_density = sp.expand(em * delta_hm_lm - alpha_minus * Lm)
    assert sp.simplify(
        plus_cocycle_density - sp.diff(ep * hp * Lp, xp)
    ) == 0
    assert sp.simplify(
        minus_cocycle_density - sp.diff(em * hm * Lm, xm)
    ) == 0
    print("central cocycles are boundary total derivatives: True")


if __name__ == "__main__":
    main()

