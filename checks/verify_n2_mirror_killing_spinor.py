"""Verify the missing primed half of the type-II SDFT Killing system.

The N=1 calculation in ``verify_hairy_killing_spinor.py`` retains only

    D_{\bar p} epsilon = 0,       gamma^p D_p epsilon = 0.

The N=2, D=10 SDFT of arXiv:1210.5078 contains in addition

    D'_p epsilon' = 0,            bar-gamma^{\bar p} D'_{\bar p} epsilon' = 0.

This script evaluates the latter pair on the most general one-sided
non-Riemannian branch

    L_+ = 0,
    L_- = L_-(x^-),
    W = W_0(x^+,x^-) + exp(-2y/l) W_1(x^+,x^-),

keeping all first derivatives of L_-, W_0 and W_1 independent.  It uses the
same ten-dimensional NS-NS uplift and the same torsionless semi-covariant
connection as the unprimed calculation.  The barred Clifford matrices obey
{bar-gamma,bar-gamma}=2 bar-eta with bar-eta=-eta.
"""

from __future__ import annotations

import sympy as sp

from verify_10d_killing_spinor import (
    block_diag,
    clifford_matrices,
    riemannian_vielbeins,
    sector_dft_metric,
    semi_covariant_connection,
    spin_connection,
)


def extremal_minus_general_vielbeins(
    u: sp.Symbol, lm: sp.Symbol, w0: sp.Symbol, w1: sp.Symbol
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Polynomial L_+=0 double vielbeins, with W=W_0+u W_1."""
    root2 = sp.sqrt(2)
    w = w0 + u * w1
    v = sp.Matrix(
        [
            [-1, 0, 0],
            [0, 0, 0],
            [0, 0, 1 / root2],
            [0, 1, 0],
            [-w / 2, -lm * u, 0],
            [0, 0, 1 / root2],
        ]
    )
    vbar = sp.Matrix(
        [
            [lm * u, 0, 0],
            [1, 0, 0],
            [0, 0, 1 / root2],
            [-w / 2, 0, 0],
            [w * lm * u / 2, 1, 0],
            [0, 0, -1 / root2],
        ]
    )
    eta = sp.Matrix([[0, -1, 0], [-1, 0, 0], [0, 0, 1]])
    return v, vbar, eta


def matrix_connection(
    phi: list[sp.Matrix], gammas: list[sp.Matrix]
) -> list[sp.Matrix]:
    bivectors = [
        [
            (gammas[left] * gammas[right]
             - gammas[right] * gammas[left]) / 2
            for right in range(10)
        ]
        for left in range(10)
    ]
    result = []
    for a_index in range(20):
        value = sp.zeros(32)
        for left in range(10):
            for right in range(10):
                coefficient = phi[a_index][left, right]
                if coefficient != 0:
                    value += coefficient * bivectors[left][right] / 4
        result.append(sp.simplify(value))
    return result


def coefficient_in(matrix: sp.Matrix, symbol: sp.Symbol) -> sp.Matrix:
    return matrix.applyfunc(lambda value: sp.factor(sp.diff(value, symbol)))


def primed_sphere_integrability_check() -> bool:
    """Check the flatness of the primed S3 parallel-spinor connection."""
    theta = sp.symbols("theta", positive=True)
    sigma1 = sp.Matrix([[0, 1], [1, 0]])
    sigma2 = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sigma3 = sp.diag(1, -1)
    bmatrix = sp.sin(theta) * sigma3 + sp.cos(theta) * sigma2
    connection = [
        -sp.I * sigma1 / 2,
        -sp.I * bmatrix / 2,
        -sp.I * bmatrix / 2,
    ]
    integrable = (
        sp.simplify(
            sp.diff(connection[1], theta)
            + connection[1] * connection[0]
            - connection[0] * connection[1]
        )
        == sp.zeros(2)
        and connection[1] == connection[2]
    )
    return integrable


def main() -> None:
    sphere_integrable = primed_sphere_integrability_check()
    l, u, lm, w0, w1 = sp.symbols(
        "l u Lm W0 W1", positive=True, nonzero=True
    )
    lmm = sp.symbols("Lm_m")
    w0p, w0m, w1p, w1m = sp.symbols("W0_p W0_m W1_p W1_m")
    theta = sp.symbols("theta", positive=True)
    theta_point = sp.pi / 6

    v_nr, vbar_nr, eta_nr = extremal_minus_general_vielbeins(
        u, lm, w0, w1
    )
    derivative_v_nr = {
        3: sp.diff(v_nr, w0) * w0p + sp.diff(v_nr, w1) * w1p,
        4: (
            sp.diff(v_nr, lm) * lmm
            + sp.diff(v_nr, w0) * w0m
            + sp.diff(v_nr, w1) * w1m
        ),
        5: sp.diff(v_nr, u) * (-2 * u / l),
    }
    derivative_vbar_nr = {
        3: sp.diff(vbar_nr, w0) * w0p + sp.diff(vbar_nr, w1) * w1p,
        4: (
            sp.diff(vbar_nr, lm) * lmm
            + sp.diff(vbar_nr, w0) * w0m
            + sp.diff(vbar_nr, w1) * w1m
        ),
        5: sp.diff(vbar_nr, u) * (-2 * u / l),
    }
    p_nr = sp.simplify(v_nr * eta_nr * v_nr.T)
    derivative_p_nr = {
        3: sp.diff(p_nr, w0) * w0p + sp.diff(p_nr, w1) * w1p,
        4: (
            sp.diff(p_nr, lm) * lmm
            + sp.diff(p_nr, w0) * w0m
            + sp.diff(p_nr, w1) * w1m
        ),
        5: sp.diff(p_nr, u) * (-2 * u / l),
    }

    metric_s = sp.diag(
        l**2, l**2 * sp.cos(theta) ** 2, l**2 * sp.sin(theta) ** 2
    )
    bfield_s = sp.zeros(3)
    bfield_s[1, 2] = l**2 * sp.cos(theta) ** 2
    bfield_s[2, 1] = -bfield_s[1, 2]
    frame_s = sp.diag(l, l * sp.cos(theta), l * sp.sin(theta))
    v_s, vbar_s, eta_s = riemannian_vielbeins(
        metric_s, bfield_s, frame_s
    )
    v_r, vbar_r, eta_r = riemannian_vielbeins(
        sp.eye(4), sp.zeros(4), sp.eye(4)
    )

    v_symbolic = block_diag(v_nr, v_s, v_r)
    vbar_symbolic = block_diag(vbar_nr, vbar_s, vbar_r)
    eta = block_diag(eta_nr, eta_s, eta_r)
    bareta = -eta
    jmetric = block_diag(
        sector_dft_metric(3), sector_dft_metric(3), sector_dft_metric(4)
    )
    p_symbolic = sp.simplify(v_symbolic * eta * v_symbolic.T)
    pbar_symbolic = sp.simplify(vbar_symbolic * bareta * vbar_symbolic.T)

    derivative_v_symbolic = {
        index: block_diag(value, sp.zeros(6, 3), sp.zeros(8, 4))
        for index, value in derivative_v_nr.items()
    }
    derivative_vbar_symbolic = {
        index: block_diag(value, sp.zeros(6, 3), sp.zeros(8, 4))
        for index, value in derivative_vbar_nr.items()
    }
    derivative_p_symbolic = {
        index: block_diag(value, sp.zeros(6), sp.zeros(8))
        for index, value in derivative_p_nr.items()
    }
    derivative_v_symbolic[9] = sp.diff(v_symbolic, theta)
    derivative_vbar_symbolic[9] = sp.diff(vbar_symbolic, theta)
    derivative_p_symbolic[9] = sp.diff(p_symbolic, theta)

    substitutions = {theta: theta_point}
    v = sp.simplify(v_symbolic.subs(substitutions))
    vbar = sp.simplify(vbar_symbolic.subs(substitutions))
    p_lower = sp.simplify(p_symbolic.subs(substitutions))
    pbar_lower = sp.simplify(pbar_symbolic.subs(substitutions))
    derivatives_v = {
        key: sp.simplify(value.subs(substitutions))
        for key, value in derivative_v_symbolic.items()
    }
    derivatives_vbar = {
        key: sp.simplify(value.subs(substitutions))
        for key, value in derivative_vbar_symbolic.items()
    }
    derivatives_p = {
        key: sp.simplify(value.subs(substitutions))
        for key, value in derivative_p_symbolic.items()
    }

    derivative_d = sp.zeros(20, 1)
    derivative_d[5] = -1 / l
    derivative_d[9] = -sp.Rational(1, 2) * (
        sp.cot(theta_point) - sp.tan(theta_point)
    )
    affine = semi_covariant_connection(
        p_lower, pbar_lower, jmetric, derivatives_p, derivative_d, 10
    )
    phibar = spin_connection(
        vbar, bareta, jmetric, affine, derivatives_vbar
    )

    gamma, chirality = clifford_matrices()
    bargamma = [sp.I * value for value in gamma]
    omega_bar = matrix_connection(phibar, bargamma)
    v_up = sp.simplify(jmetric * v * eta)
    vbar_up = sp.simplify(jmetric * vbar * bareta)

    # Solve the three primed S3 gravitino equations D'_p epsilon'=0 for
    # the physical-coordinate spinor derivatives.  This avoids assuming the
    # sign of the torsionful S3 connection and tests it from the SDFT data.
    physical_sphere = (9, 10, 11)
    flat_sphere = (3, 4, 5)
    coefficient = sp.Matrix(
        [
            [v_up[a_index, p_index] for a_index in physical_sphere]
            for p_index in flat_sphere
        ]
    )
    inverse = sp.simplify(coefficient.inv())
    nonderivative = []
    for p_index in flat_sphere:
        value = sp.zeros(32)
        for a_index in range(20):
            factor = v_up[a_index, p_index]
            if factor != 0:
                value += factor * omega_bar[a_index]
        nonderivative.append(sp.simplify(value))
    sphere_derivatives = {}
    for mu_index, a_index in enumerate(physical_sphere):
        value = sp.zeros(32)
        for p_local in range(3):
            value -= inverse[mu_index, p_local] * nonderivative[p_local]
        sphere_derivatives[a_index] = sp.simplify(value)

    # Extract the two-by-two S3 factor when the result is of the expected
    # I_NR x K_S3 x I_aux x I_R4 form.
    sphere_factor = {}
    for a_index, value in sphere_derivatives.items():
        factor = sp.Matrix(
            2, 2, lambda row, column: value[row * 8, column * 8]
        )
        reconstructed = sp.kronecker_product(
            sp.eye(2), factor, sp.eye(2), sp.eye(4)
        )
        sphere_factor[a_index] = (
            factor,
            sp.simplify(value - reconstructed) == sp.zeros(32),
        )

    external_symbols = {
        lm, w0, w1, lmm, w0p, w0m, w1p, w1m, u
    }
    sphere_external_dependence = set()
    for value in sphere_derivatives.values():
        sphere_external_dependence |= value.free_symbols & external_symbols

    e0, e1, ep0, ep1, em0, em1, ey0, ey1 = sp.symbols(
        "e0 e1 ep0 ep1 em0 em1 ey0 ey1"
    )
    variables = [ep0, ep1, em0, em1, ey0, ey1, e0, e1]
    internal = sp.eye(16)
    spinor = sp.kronecker_product(sp.Matrix([e0, e1]), internal)
    partials = {
        3: sp.kronecker_product(sp.Matrix([ep0, ep1]), internal),
        4: sp.kronecker_product(sp.Matrix([em0, em1]), internal),
        5: sp.kronecker_product(sp.Matrix([ey0, ey1]), internal),
    }
    covariant = []
    for a_index in range(20):
        derivative = partials.get(a_index, sp.zeros(32, 16))
        if a_index in sphere_derivatives:
            derivative += sphere_derivatives[a_index] * spinor
        covariant.append(sp.simplify(derivative + omega_bar[a_index] * spinor))

    residuals = []
    for p_index in range(10):
        equation = sp.zeros(32, 16)
        for a_index in range(20):
            factor = v_up[a_index, p_index]
            if factor != 0:
                equation += factor * covariant[a_index]
        residuals.append(sp.simplify(equation))
    dilatino = sp.zeros(32, 16)
    for bar_index in range(10):
        projected = sp.zeros(32, 16)
        for a_index in range(20):
            factor = vbar_up[a_index, bar_index]
            if factor != 0:
                projected += factor * covariant[a_index]
        dilatino += bargamma[bar_index] * projected
    residuals.append(sp.simplify(dilatino))

    # Test epsilon'=(F(x^-),0) for arbitrary F.  The coefficients of F and
    # F' must annihilate the same internal polarization space.
    candidate_constraints = []
    for residual in residuals:
        for variable in (e0, em0):
            value = coefficient_in(residual, variable)
            if value != sp.zeros(32, 16):
                candidate_constraints.append(value)
    stacked = sp.Matrix.vstack(*candidate_constraints)
    raw_rank = stacked.rank()

    # The barred chirality convention can differ by an overall sign when
    # bar-gamma=i gamma.  Report both choices; each physical choice should
    # retain four real internal polarizations.
    chirality_internal = chirality.extract(range(16), range(16))
    results = []
    for chirality_sign in (1, -1):
        constrained = stacked.col_join(
            chirality_internal - chirality_sign * sp.eye(16)
        )
        results.append(
            (chirality_sign, constrained.rank(), constrained.nullspace())
        )

    # Verify the full first-order row space after restriction to every
    # surviving internal basis, not merely the proposed candidate fields.
    print("=== N=2 primed full-D=10 mirror check ===")
    print("primed S3 connection integrable for all theta:", sphere_integrable)
    print("P+Pbar=J:", p_lower + pbar_lower == jmetric)
    print("sphere derivative external dependence:", sphere_external_dependence)
    print("primed S3 derivative factors at theta=pi/6:")
    for a_index in physical_sphere:
        print(" ", a_index, sphere_factor[a_index])
    print("candidate equations:", stacked.rows)
    print("candidate rank/nullity before Weyl:", raw_rank, 16 - raw_rank)
    for chirality_sign, rank, basis in results:
        print(
            "bar-chirality", chirality_sign,
            "rank/nullity:", rank, 16 - rank,
            "basis columns:",
            [
                [index for index, entry in enumerate(vector) if entry != 0]
                for vector in basis
            ],
        )

    for chirality_sign, _, basis in results:
        if not basis:
            continue
        candidate = sp.Matrix.hstack(*basis)
        nonzero = []
        substitution = {
            e1: 0,
            ep0: 0,
            ep1: 0,
            em1: 0,
            ey0: 0,
            ey1: 0,
        }
        for equation_index, residual in enumerate(residuals):
            value = sp.simplify(residual.subs(substitution) * candidate)
            if value != sp.zeros(32, len(basis)):
                nonzero.append((equation_index, value))
        print(
            "bar-chirality", chirality_sign,
            "residual count on epsilon'=(F(x^-),0):", len(nonzero),
        )


if __name__ == "__main__":
    main()

