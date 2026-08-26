"""Killing-spinor system on the constant-(L_+,L_-,W_1) NR saddle.

This extends ``verify_10d_killing_spinor.py`` away from the asymptotic vacuum.
The ten-dimensional S3 x R4 factor supplies the same fermion shift as in the
vacuum calculation, so the nontrivial radial problem can first be exposed in
the three-dimensional equations

    D_bar{p} E = 0,
    gamma^p D_p E = -E/(sqrt(2) l).

The script constructs the exact non-Riemannian double vielbein, derives the
semi-covariant spin connection, and prints the reduced first-order system for
constant positive L_+, L_- and W = W_1 h/sqrt(L_+ L_-), h=chi/2.
"""

from __future__ import annotations

import sys
import itertools

import sympy as sp

from verify_10d_killing_spinor import (
    block_diag,
    build_ten_dimensional,
    clifford_matrices,
    dft_metric,
    riemannian_vielbeins,
    rotate_nr_vielbeins_to_kim,
    semi_covariant_connection,
    sector_dft_metric,
    spin_connection,
)


def exact_nr_vielbeins(
    h: sp.Symbol, a: sp.Symbol, w: sp.Symbol
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Return Eq. (SM21)'s NR double vielbeins in (dual; physical) order."""
    root2 = sp.sqrt(2)
    ch = sp.cosh(h)
    sh = sp.sinh(h)
    v = sp.Matrix(
        [
            [-ch, 0, 0],
            [-a * sh, 0, 0],
            [0, 0, 1 / root2],
            [w * a * sh / 2, ch, 0],
            [-w * ch / 2, -sh / a, 0],
            [0, 0, 1 / root2],
        ]
    )
    vbar = sp.Matrix(
        [
            [sh / a, 0, 0],
            [ch, 0, 0],
            [0, 0, 1 / root2],
            [-w * ch / 2, -a * sh, 0],
            [w * sh / (2 * a), ch, 0],
            [0, 0, -1 / root2],
        ]
    )
    eta = sp.Matrix([[0, -1, 0], [-1, 0, 0], [0, 0, 1]])
    v, vbar = rotate_nr_vielbeins_to_kim(v, vbar)
    return v, vbar, eta


def extremal_plus_vielbeins(
    u: sp.Symbol, lp: sp.Symbol, w1: sp.Symbol
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Nonsingular L_-=0 limit, with W=W_1 u and u=e^{-2y/l}."""
    root2 = sp.sqrt(2)
    w = w1 * u
    v = sp.Matrix(
        [
            [-1, 0, 0],
            [-lp * u, 0, 0],
            [0, 0, 1 / root2],
            [w * lp * u / 2, 1, 0],
            [-w / 2, 0, 0],
            [0, 0, 1 / root2],
        ]
    )
    vbar = sp.Matrix(
        [
            [0, 0, 0],
            [1, 0, 0],
            [0, 0, 1 / root2],
            [-w / 2, -lp * u, 0],
            [0, 1, 0],
            [0, 0, -1 / root2],
        ]
    )
    eta = sp.Matrix([[0, -1, 0], [-1, 0, 0], [0, 0, 1]])
    v, vbar = rotate_nr_vielbeins_to_kim(v, vbar)
    return v, vbar, eta


def extremal_plus_general_vielbeins(
    u: sp.Symbol, lp: sp.Symbol, w0: sp.Symbol, w1: sp.Symbol
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Most general nonsingular L_-=0 homogeneous-hair vielbein."""
    root2 = sp.sqrt(2)
    w = w0 + w1 * u
    v = sp.Matrix(
        [
            [-1, 0, 0],
            [-lp * u, 0, 0],
            [0, 0, 1 / root2],
            [w * lp * u / 2, 1, 0],
            [-w / 2, 0, 0],
            [0, 0, 1 / root2],
        ]
    )
    vbar = sp.Matrix(
        [
            [0, 0, 0],
            [1, 0, 0],
            [0, 0, 1 / root2],
            [-w / 2, -lp * u, 0],
            [0, 1, 0],
            [0, 0, -1 / root2],
        ]
    )
    eta = sp.Matrix([[0, -1, 0], [-1, 0, 0], [0, 0, 1]])
    v, vbar = rotate_nr_vielbeins_to_kim(v, vbar)
    return v, vbar, eta


def extremal_minus_vielbeins(
    u: sp.Symbol, lm: sp.Symbol, w1: sp.Symbol
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Nonsingular L_+=0 limit, with W=W_1 u and u=e^{-2y/l}."""
    root2 = sp.sqrt(2)
    w = w1 * u
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
    v, vbar = rotate_nr_vielbeins_to_kim(v, vbar)
    return v, vbar, eta


def coefficient_matrix(expr: sp.Matrix, variables: list[sp.Symbol]) -> sp.Matrix:
    """Return the linear coefficients of a matrix expression."""
    rows = []
    for entry in list(expr):
        rows.append([sp.factor(sp.diff(entry, variable)) for variable in variables])
    return sp.Matrix(rows)


def build_reduced_system() -> None:
    l, h, a, p, w1 = sp.symbols("l h a p w1", positive=True, nonzero=True)
    q = sp.tanh(h / sp.sqrt(2))
    h_y = -2 * sp.sqrt(2) * q / (l * (1 - q**2))
    w = w1 * h / p
    d_y = -1 / l + h_y * q / sp.sqrt(2)

    v, vbar, eta = exact_nr_vielbeins(h, a, w)
    jmetric = dft_metric(3)
    p_lower = sp.simplify(v * eta * v.T)
    pbar_lower = sp.simplify(vbar * (-eta) * vbar.T)
    # w=w(h), so SymPy's total h derivative already includes w_y.
    derivative_v_y = sp.diff(v, h) * h_y
    derivative_p_y = sp.diff(p_lower, h) * h_y
    derivative_d = sp.zeros(6, 1)
    derivative_d[5] = d_y
    gamma = semi_covariant_connection(
        p_lower, pbar_lower, jmetric, {5: derivative_p_y}, derivative_d, 3
    )
    phi = spin_connection(v, eta, jmetric, gamma, {5: derivative_v_y})

    root2 = sp.sqrt(2)
    gamma_flat = [
        sp.Matrix([[0, 0], [root2, 0]]),
        sp.Matrix([[0, -root2], [0, 0]]),
        sp.diag(1, -1),
    ]
    gamma_bivectors = [
        [
            (gamma_flat[left] * gamma_flat[right]
             - gamma_flat[right] * gamma_flat[left]) / 2
            for right in range(3)
        ]
        for left in range(3)
    ]
    omega = []
    for a_index in range(6):
        matrix = sp.zeros(2)
        for left in range(3):
            for right in range(3):
                matrix += phi[a_index][left, right] * gamma_bivectors[left][right] / 4
        omega.append(sp.simplify(matrix))

    v_up = sp.simplify(jmetric * v * eta)
    vbar_up = sp.simplify(jmetric * vbar * (-eta))
    e0, e1 = sp.symbols("e0 e1")
    ep0, ep1 = sp.symbols("ep0 ep1")
    em0, em1 = sp.symbols("em0 em1")
    ey0, ey1 = sp.symbols("ey0 ey1")
    spinor = sp.Matrix([e0, e1])
    partials = {
        3: sp.Matrix([ep0, ep1]),
        4: sp.Matrix([em0, em1]),
        5: sp.Matrix([ey0, ey1]),
    }
    covariant = [partials.get(index, sp.zeros(2, 1)) + omega[index] * spinor for index in range(6)]

    equations = []
    labels = []
    for bar_index in range(3):
        residual = sp.zeros(2, 1)
        for a_index in range(6):
            residual += vbar_up[a_index, bar_index] * covariant[a_index]
        equations.append(sp.simplify(residual))
        labels.append(f"gravitino[{bar_index}]")

    dirac = sp.zeros(2, 1)
    for flat_index in range(3):
        projected = sp.zeros(2, 1)
        for a_index in range(6):
            projected += v_up[a_index, flat_index] * covariant[a_index]
        dirac += gamma_flat[flat_index] * projected
    equations.append(sp.simplify(dirac + spinor / (root2 * l)))
    labels.append("shifted dilatino")

    variables = [ep0, ep1, em0, em1, ey0, ey1, e0, e1]
    full = sp.Matrix.vstack(*equations)
    system = coefficient_matrix(full, variables)

    print("=== Exact reduced Killing system ===")
    print("variables:", variables)
    for label, equation in zip(labels, equations):
        print(label)
        for component in equation:
            print(" ", sp.factor(component))

    # Evaluate the pointwise differential system at generic regular values.
    # h=log(2) makes sinh(h)=3/4 and cosh(h)=5/4 exactly.
    generic = {
        l: 1,
        h: sp.log(2),
        a: sp.Rational(3, 2),
        p: sp.Rational(5, 4),
        w1: sp.Rational(7, 5),
    }
    matrix_generic = sp.N(system.subs(generic), 50)
    print("pointwise equation rank at generic nonzero L+/L-/W1 =", matrix_generic.rank())

    # The algebraic compatibility conditions are the left-null vectors of the
    # six derivative columns.  Acting on the last two columns gives constraints
    # that a nonzero spinor must satisfy independently of its first derivatives.
    derivative_block = matrix_generic[:, :6]
    algebraic_block = matrix_generic[:, 6:]
    left_kernel = derivative_block.T.nullspace()
    obstruction = sp.Matrix.vstack(*[(vector.T * algebraic_block) for vector in left_kernel])
    print("number of pointwise derivative-independent constraints =", obstruction.rows)
    print("rank on the two spinor components =", obstruction.rank())
    print("constraint matrix (numeric) =")
    print(obstruction.evalf(12))

    first_constraint = sp.factor(equations[0][0] / e1)
    second_after_first = sp.factor(
        equations[1][1].subs({e1: 0, ep1: 0, em1: 0, ey1: 0}) / e0
    )
    print("first smooth-field constraint coefficient on e1 =", first_constraint)
    print("subsequent constraint coefficient on e0 =", second_after_first)
    print(
        "generic non-extremal conclusion: e1=0 and then e0=0 for h != 0; "
        "W1 does not remove the obstruction"
    )


def build_hair_only_system() -> None:
    """Check the separate L_+=L_-=0, W=W_1 e^{-2y/l} limiting branch."""
    l, w = sp.symbols("l w", positive=True, nonzero=True)
    w_y = -2 * w / l
    v, vbar, eta = exact_nr_vielbeins(sp.Integer(0), sp.Integer(1), w)
    jmetric = dft_metric(3)
    p_lower = sp.simplify(v * eta * v.T)
    pbar_lower = sp.simplify(vbar * (-eta) * vbar.T)
    derivative_v_y = sp.diff(v, w) * w_y
    derivative_p_y = sp.diff(p_lower, w) * w_y
    derivative_d = sp.zeros(6, 1)
    derivative_d[5] = -1 / l
    gamma = semi_covariant_connection(
        p_lower, pbar_lower, jmetric, {5: derivative_p_y}, derivative_d, 3
    )
    phi = spin_connection(v, eta, jmetric, gamma, {5: derivative_v_y})

    root2 = sp.sqrt(2)
    gamma_flat = [
        sp.Matrix([[0, 0], [root2, 0]]),
        sp.Matrix([[0, -root2], [0, 0]]),
        sp.diag(1, -1),
    ]
    gamma_bivectors = [
        [
            (gamma_flat[left] * gamma_flat[right]
             - gamma_flat[right] * gamma_flat[left]) / 2
            for right in range(3)
        ]
        for left in range(3)
    ]
    omega = []
    for a_index in range(6):
        matrix = sp.zeros(2)
        for left in range(3):
            for right in range(3):
                matrix += phi[a_index][left, right] * gamma_bivectors[left][right] / 4
        omega.append(sp.simplify(matrix))

    v_up = sp.simplify(jmetric * v * eta)
    vbar_up = sp.simplify(jmetric * vbar * (-eta))
    e0, e1 = sp.symbols("e0 e1")
    ep0, ep1 = sp.symbols("ep0 ep1")
    em0, em1 = sp.symbols("em0 em1")
    ey0, ey1 = sp.symbols("ey0 ey1")
    spinor = sp.Matrix([e0, e1])
    partials = {
        3: sp.Matrix([ep0, ep1]),
        4: sp.Matrix([em0, em1]),
        5: sp.Matrix([ey0, ey1]),
    }
    covariant = [partials.get(index, sp.zeros(2, 1)) + omega[index] * spinor for index in range(6)]
    equations = []
    labels = []
    for bar_index in range(3):
        residual = sp.zeros(2, 1)
        for a_index in range(6):
            residual += vbar_up[a_index, bar_index] * covariant[a_index]
        equations.append(sp.simplify(residual))
        labels.append(f"gravitino[{bar_index}]")
    dirac = sp.zeros(2, 1)
    for flat_index in range(3):
        projected = sp.zeros(2, 1)
        for a_index in range(6):
            projected += v_up[a_index, flat_index] * covariant[a_index]
        dirac += gamma_flat[flat_index] * projected
    equations.append(sp.simplify(dirac + spinor / (root2 * l)))
    labels.append("shifted dilatino")

    print("\n=== Hair-only limiting branch ===")
    for label, equation in zip(labels, equations):
        print(label)
        for component in equation:
            print(" ", sp.factor(component))


def full_ten_dimensional_hair_only() -> None:
    """Check the hair-only no-go directly in the full D=10 product."""
    l, w = sp.symbols("l w", positive=True, nonzero=True)
    theta = sp.symbols("theta", positive=True)
    theta_point = sp.pi / 6
    w_y = -2 * w / l

    v_nr, vbar_nr, eta_nr = exact_nr_vielbeins(sp.Integer(0), sp.Integer(1), w)
    metric_s = sp.diag(l**2, l**2 * sp.cos(theta) ** 2, l**2 * sp.sin(theta) ** 2)
    bfield_s = sp.zeros(3)
    bfield_s[1, 2] = l**2 * sp.cos(theta) ** 2
    bfield_s[2, 1] = -bfield_s[1, 2]
    frame_s = sp.diag(l, l * sp.cos(theta), l * sp.sin(theta))
    v_s, vbar_s, eta_s = riemannian_vielbeins(metric_s, bfield_s, frame_s)
    v_r, vbar_r, eta_r = riemannian_vielbeins(sp.eye(4), sp.zeros(4), sp.eye(4))

    v_symbolic = block_diag(v_nr, v_s, v_r)
    vbar_symbolic = block_diag(vbar_nr, vbar_s, vbar_r)
    eta = block_diag(eta_nr, eta_s, eta_r)
    jmetric = block_diag(sector_dft_metric(3), sector_dft_metric(3), sector_dft_metric(4))
    p_symbolic = sp.simplify(v_symbolic * eta * v_symbolic.T)
    pbar_symbolic = sp.simplify(vbar_symbolic * (-eta) * vbar_symbolic.T)

    derivative_v_y = sp.diff(v_symbolic, w) * w_y
    derivative_p_y = sp.diff(p_symbolic, w) * w_y
    derivative_v_theta = sp.diff(v_symbolic, theta)
    derivative_p_theta = sp.diff(p_symbolic, theta)
    substitutions = {theta: theta_point}
    v = sp.simplify(v_symbolic.subs(substitutions))
    vbar = sp.simplify(vbar_symbolic.subs(substitutions))
    p_lower = sp.simplify(p_symbolic.subs(substitutions))
    pbar_lower = sp.simplify(pbar_symbolic.subs(substitutions))
    derivative_v_y = sp.simplify(derivative_v_y.subs(substitutions))
    derivative_p_y = sp.simplify(derivative_p_y.subs(substitutions))
    derivative_v_theta = sp.simplify(derivative_v_theta.subs(substitutions))
    derivative_p_theta = sp.simplify(derivative_p_theta.subs(substitutions))
    derivative_d = sp.zeros(20, 1)
    derivative_d[5] = -1 / l
    derivative_d[9] = -sp.Rational(1, 2) * (
        sp.cot(theta_point) - sp.tan(theta_point)
    )
    gamma = semi_covariant_connection(
        p_lower,
        pbar_lower,
        jmetric,
        {5: derivative_p_y, 9: derivative_p_theta},
        derivative_d,
        10,
    )
    phi = spin_connection(
        v,
        eta,
        jmetric,
        gamma,
        {5: derivative_v_y, 9: derivative_v_theta},
    )

    gammas, _ = clifford_matrices()
    gamma_bivectors = [
        [
            (gammas[left] * gammas[right] - gammas[right] * gammas[left]) / 2
            for right in range(10)
        ]
        for left in range(10)
    ]
    omega = []
    for a_index in range(20):
        matrix = sp.zeros(32)
        for left in range(10):
            for right in range(10):
                if phi[a_index][left, right] != 0:
                    matrix += phi[a_index][left, right] * gamma_bivectors[left][right] / 4
        omega.append(sp.simplify(matrix))

    sigma1 = sp.Matrix([[0, 1], [1, 0]])
    sigma2 = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sigma3 = sp.diag(1, -1)
    candidate = sp.zeros(16, 4)
    for column, row in enumerate((5, 6, 13, 14)):
        candidate[row, column] = 1
    a_sphere = sp.sin(theta_point) * sigma3 - sp.cos(theta_point) * sigma2
    sphere_derivatives = {
        9: sp.I * sigma1 / 2,
        10: -sp.I * a_sphere / 2,
        11: sp.I * a_sphere / 2,
    }

    e0, e1 = sp.symbols("e0 e1")
    ep0, ep1 = sp.symbols("ep0 ep1")
    em0, em1 = sp.symbols("em0 em1")
    ey0, ey1 = sp.symbols("ey0 ey1")
    variables = [ep0, ep1, em0, em1, ey0, ey1, e0, e1]
    spinor = sp.kronecker_product(sp.Matrix([e0, e1]), candidate)
    partial_nr = {
        3: sp.kronecker_product(sp.Matrix([ep0, ep1]), candidate),
        4: sp.kronecker_product(sp.Matrix([em0, em1]), candidate),
        5: sp.kronecker_product(sp.Matrix([ey0, ey1]), candidate),
    }
    covariant = []
    for a_index in range(20):
        derivative = partial_nr.get(a_index, sp.zeros(32, 4))
        if a_index in sphere_derivatives:
            internal_derivative = sp.kronecker_product(
                sphere_derivatives[a_index], sp.eye(2), sp.eye(4)
            )
            derivative += sp.kronecker_product(
                sp.Matrix([e0, e1]), internal_derivative * candidate
            )
        covariant.append(derivative + omega[a_index] * spinor)

    v_up = sp.simplify(jmetric * v * eta)
    vbar_up = sp.simplify(jmetric * vbar * (-eta))
    residuals = []
    for bar_index in range(10):
        equation = sp.zeros(32, 4)
        for a_index in range(20):
            if vbar_up[a_index, bar_index] != 0:
                equation += vbar_up[a_index, bar_index] * covariant[a_index]
        residuals.append(equation)
    dilatino = sp.zeros(32, 4)
    for flat_index in range(10):
        projected = sp.zeros(32, 4)
        for a_index in range(20):
            if v_up[a_index, flat_index] != 0:
                projected += v_up[a_index, flat_index] * covariant[a_index]
        dilatino += gammas[flat_index] * projected
    residuals.append(dilatino)

    rows = []
    for residual in residuals:
        for entry in list(residual):
            if entry != 0:
                rows.append([sp.factor(sp.diff(entry, variable)) for variable in variables])
    system = sp.Matrix(rows)
    print("\n=== Full D=10 hair-only check at theta=pi/6 ===")
    print("nonzero component equations =", system.rows)
    print("symbolic coefficient-matrix rank =", system.rank())
    reduced_rref, pivots = system.rref()
    print("pivot columns =", pivots)
    print("independent equations (RREF) =")
    for row in reduced_rref.tolist():
        if any(value != 0 for value in row):
            print(" ", row)


def full_ten_dimensional_nonextremal_point() -> None:
    """Direct D=10 check at one regular L_+ L_- W_1 != 0 bulk point."""
    l, h, a, p, w1 = sp.symbols("l h a p w1", positive=True, nonzero=True)
    theta = sp.symbols("theta", positive=True)
    q = sp.tanh(h / sp.sqrt(2))
    h_y = -2 * sp.sqrt(2) * q / (l * (1 - q**2))
    w = w1 * h / p
    d_y = -1 / l + h_y * q / sp.sqrt(2)

    v_nr, vbar_nr, eta_nr = exact_nr_vielbeins(h, a, w)
    derivative_v_nr_y = sp.diff(v_nr, h) * h_y
    p_nr = sp.simplify(v_nr * eta_nr * v_nr.T)
    derivative_p_nr_y = sp.diff(p_nr, h) * h_y

    metric_s = sp.diag(l**2, l**2 * sp.cos(theta) ** 2, l**2 * sp.sin(theta) ** 2)
    bfield_s = sp.zeros(3)
    bfield_s[1, 2] = l**2 * sp.cos(theta) ** 2
    bfield_s[2, 1] = -bfield_s[1, 2]
    frame_s = sp.diag(l, l * sp.cos(theta), l * sp.sin(theta))
    v_s, vbar_s, eta_s = riemannian_vielbeins(metric_s, bfield_s, frame_s)
    v_r, vbar_r, eta_r = riemannian_vielbeins(sp.eye(4), sp.zeros(4), sp.eye(4))

    v_symbolic = block_diag(v_nr, v_s, v_r)
    vbar_symbolic = block_diag(vbar_nr, vbar_s, vbar_r)
    eta = block_diag(eta_nr, eta_s, eta_r)
    jmetric = block_diag(sector_dft_metric(3), sector_dft_metric(3), sector_dft_metric(4))
    p_symbolic = sp.simplify(v_symbolic * eta * v_symbolic.T)
    pbar_symbolic = sp.simplify(vbar_symbolic * (-eta) * vbar_symbolic.T)
    derivative_v_y_symbolic = block_diag(derivative_v_nr_y, sp.zeros(6, 3), sp.zeros(8, 4))
    derivative_p_y_symbolic = block_diag(derivative_p_nr_y, sp.zeros(6), sp.zeros(8))
    derivative_v_theta_symbolic = sp.diff(v_symbolic, theta)
    derivative_p_theta_symbolic = sp.diff(p_symbolic, theta)

    # q=1/2 is a regular bulk point.  Together with a=3/2 and p=5/4 it
    # corresponds to L_+=15/8 and L_-=5/6.
    h_point = sp.sqrt(2) * sp.atanh(sp.Rational(1, 2))
    substitutions = {
        l: 1,
        h: h_point,
        a: sp.Rational(3, 2),
        p: sp.Rational(5, 4),
        w1: sp.Rational(7, 5),
        theta: sp.pi / 6,
    }
    v = sp.simplify(v_symbolic.subs(substitutions))
    vbar = sp.simplify(vbar_symbolic.subs(substitutions))
    p_lower = sp.simplify(p_symbolic.subs(substitutions))
    pbar_lower = sp.simplify(pbar_symbolic.subs(substitutions))
    derivative_v_y = sp.simplify(derivative_v_y_symbolic.subs(substitutions))
    derivative_p_y = sp.simplify(derivative_p_y_symbolic.subs(substitutions))
    derivative_v_theta = sp.simplify(derivative_v_theta_symbolic.subs(substitutions))
    derivative_p_theta = sp.simplify(derivative_p_theta_symbolic.subs(substitutions))
    derivative_d = sp.zeros(20, 1)
    derivative_d[5] = sp.simplify(d_y.subs(substitutions))
    derivative_d[9] = -sp.Rational(1, 2) * (
        sp.cot(sp.pi / 6) - sp.tan(sp.pi / 6)
    )
    gamma = semi_covariant_connection(
        p_lower,
        pbar_lower,
        jmetric,
        {5: derivative_p_y, 9: derivative_p_theta},
        derivative_d,
        10,
    )
    phi = spin_connection(
        v,
        eta,
        jmetric,
        gamma,
        {5: derivative_v_y, 9: derivative_v_theta},
    )

    gammas, _ = clifford_matrices()
    gamma_bivectors = [
        [
            (gammas[left] * gammas[right] - gammas[right] * gammas[left]) / 2
            for right in range(10)
        ]
        for left in range(10)
    ]
    omega = []
    for a_index in range(20):
        matrix = sp.zeros(32)
        for left in range(10):
            for right in range(10):
                if phi[a_index][left, right] != 0:
                    matrix += phi[a_index][left, right] * gamma_bivectors[left][right] / 4
        omega.append(sp.simplify(matrix))

    sigma1 = sp.Matrix([[0, 1], [1, 0]])
    sigma2 = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sigma3 = sp.diag(1, -1)
    candidate = sp.zeros(16, 4)
    for column, row in enumerate((5, 6, 13, 14)):
        candidate[row, column] = 1
    theta_point = sp.pi / 6
    a_sphere = sp.sin(theta_point) * sigma3 - sp.cos(theta_point) * sigma2
    sphere_derivatives = {
        9: sp.I * sigma1 / 2,
        10: -sp.I * a_sphere / 2,
        11: sp.I * a_sphere / 2,
    }

    e0, e1 = sp.symbols("e0 e1")
    ep0, ep1 = sp.symbols("ep0 ep1")
    em0, em1 = sp.symbols("em0 em1")
    ey0, ey1 = sp.symbols("ey0 ey1")
    variables = [ep0, ep1, em0, em1, ey0, ey1, e0, e1]
    spinor = sp.kronecker_product(sp.Matrix([e0, e1]), candidate)
    partial_nr = {
        3: sp.kronecker_product(sp.Matrix([ep0, ep1]), candidate),
        4: sp.kronecker_product(sp.Matrix([em0, em1]), candidate),
        5: sp.kronecker_product(sp.Matrix([ey0, ey1]), candidate),
    }
    covariant = []
    for a_index in range(20):
        derivative = partial_nr.get(a_index, sp.zeros(32, 4))
        if a_index in sphere_derivatives:
            internal_derivative = sp.kronecker_product(
                sphere_derivatives[a_index], sp.eye(2), sp.eye(4)
            )
            derivative += sp.kronecker_product(
                sp.Matrix([e0, e1]), internal_derivative * candidate
            )
        covariant.append(derivative + omega[a_index] * spinor)

    v_up = sp.simplify(jmetric * v * eta)
    vbar_up = sp.simplify(jmetric * vbar * (-eta))
    residuals = []
    for bar_index in range(10):
        equation = sp.zeros(32, 4)
        for a_index in range(20):
            if vbar_up[a_index, bar_index] != 0:
                equation += vbar_up[a_index, bar_index] * covariant[a_index]
        residuals.append(equation)
    dilatino = sp.zeros(32, 4)
    for flat_index in range(10):
        projected = sp.zeros(32, 4)
        for a_index in range(20):
            if v_up[a_index, flat_index] != 0:
                projected += v_up[a_index, flat_index] * covariant[a_index]
        dilatino += gammas[flat_index] * projected
    residuals.append(dilatino)

    rows = []
    for residual in residuals:
        for entry in list(residual):
            if entry != 0:
                rows.append([sp.factor(sp.diff(entry, variable)) for variable in variables])
    system = sp.Matrix(rows)
    print("\n=== Full D=10 non-extremal point check ===")
    print("point: q=1/2, L+=15/8, L-=5/6, W1=7/5, theta=pi/6")
    print("nonzero component equations =", system.rows)
    print("coefficient-matrix rank =", system.rank())
    reduced_rref, pivots = system.rref()
    print("pivot columns =", pivots)
    print("independent equations (RREF) =")
    for row in reduced_rref.tolist():
        if any(value != 0 for value in row):
            print(" ", [sp.N(value, 10) for value in row])


def asymptotic_general_profile() -> None:
    """First subleading equations for arbitrary L_+(x+), L_-(x-), W_1(x)."""
    l, u = sp.symbols("l u", positive=True, nonzero=True)
    lp, lm, wh = sp.symbols("Lp Lm W1")
    lpp, lmm, whp, whm = sp.symbols("Lp_p Lm_m W1_p W1_m")
    root2 = sp.sqrt(2)
    eta = sp.Matrix([[0, -1, 0], [-1, 0, 0], [0, 0, 1]])
    jmetric = dft_metric(3)

    # Exact first-order expansion of Eq. (SM21), using
    # e^sigma sqrt(Pi)=L_+ and e^-sigma sqrt(Pi)=L_-.
    v0, vbar0, _ = exact_nr_vielbeins(sp.Integer(0), sp.Integer(1), sp.Integer(0))
    v1 = sp.zeros(6, 3)
    v1[1, 0] = -lp
    v1[4, 0] = -wh / 2
    v1[4, 1] = -lm
    vbar1 = sp.zeros(6, 3)
    vbar1[0, 0] = lm
    vbar1[3, 0] = -wh / 2
    vbar1[3, 1] = -lp
    dv_plus = sp.zeros(6, 3)
    dv_plus[1, 0] = -lpp * u
    dv_plus[4, 0] = -whp * u / 2
    dv_minus = sp.zeros(6, 3)
    dv_minus[4, 0] = -whm * u / 2
    dv_minus[4, 1] = -lmm * u
    dvbar_plus = sp.zeros(6, 3)
    dvbar_plus[3, 0] = -whp * u / 2
    dvbar_plus[3, 1] = -lpp * u
    dvbar_minus = sp.zeros(6, 3)
    dvbar_minus[0, 0] = lmm * u
    dvbar_minus[3, 0] = -whm * u / 2
    v1, vbar1 = rotate_nr_vielbeins_to_kim(v1, vbar1)
    dv_plus, dvbar_plus = rotate_nr_vielbeins_to_kim(dv_plus, dvbar_plus)
    dv_minus, dvbar_minus = rotate_nr_vielbeins_to_kim(dv_minus, dvbar_minus)
    v = v0 + u * v1
    vbar = vbar0 + u * vbar1
    p_lower = sp.expand(v * eta * v.T)
    pbar_lower = sp.expand(vbar * (-eta) * vbar.T)
    dv_y = -2 * u * v1 / l
    dvbar_y = -2 * u * vbar1 / l

    # Differentiate P from V so that the coset completion is retained to the
    # order needed here.
    dp_plus = sp.expand(dv_plus * eta * v.T + v * eta * dv_plus.T)
    dp_minus = sp.expand(dv_minus * eta * v.T + v * eta * dv_minus.T)
    dp_y = sp.expand(dv_y * eta * v.T + v * eta * dv_y.T)
    derivative_d = sp.zeros(6, 1)
    derivative_d[5] = -1 / l

    gamma = semi_covariant_connection(
        p_lower,
        pbar_lower,
        jmetric,
        {3: dp_plus, 4: dp_minus, 5: dp_y},
        derivative_d,
        3,
    )

    def linear(matrix: sp.Matrix) -> sp.Matrix:
        return matrix.applyfunc(
            lambda value: sp.factor(value.subs(u, 0) + u * sp.diff(value, u).subs(u, 0))
        )

    gamma = [linear(matrix) for matrix in gamma]
    phi = spin_connection(
        v,
        eta,
        jmetric,
        gamma,
        {3: dv_plus, 4: dv_minus, 5: dv_y},
    )
    phi = [linear(matrix) for matrix in phi]

    gamma_flat = [
        sp.Matrix([[0, 0], [root2, 0]]),
        sp.Matrix([[0, -root2], [0, 0]]),
        sp.diag(1, -1),
    ]
    gamma_bivectors = [
        [
            (gamma_flat[left] * gamma_flat[right]
             - gamma_flat[right] * gamma_flat[left]) / 2
            for right in range(3)
        ]
        for left in range(3)
    ]
    omega = []
    for a_index in range(6):
        matrix = sp.zeros(2)
        for left in range(3):
            for right in range(3):
                matrix += phi[a_index][left, right] * gamma_bivectors[left][right] / 4
        omega.append(linear(matrix))

    f0, f1, f2 = sp.symbols("f f_p f_pp")
    g0, g1 = sp.symbols("g0 g1")
    gp0, gp1, gm0, gm1 = sp.symbols("g0_p g1_p g0_m g1_m")
    e0 = sp.Matrix([-root2 * f0, l * f1])
    de0_plus = sp.Matrix([-root2 * f1, l * f2])
    g = sp.Matrix([g0, g1])
    spinor = e0 + u * g
    partials = {
        3: de0_plus + u * sp.Matrix([gp0, gp1]),
        4: u * sp.Matrix([gm0, gm1]),
        5: -2 * u * g / l,
    }
    covariant = [
        linear(partials.get(index, sp.zeros(2, 1)) + omega[index] * spinor)
        for index in range(6)
    ]
    v_up = linear(jmetric * v * eta)
    vbar_up = linear(jmetric * vbar * (-eta))
    equations = []
    for bar_index in range(3):
        residual = sp.zeros(2, 1)
        for a_index in range(6):
            residual += vbar_up[a_index, bar_index] * covariant[a_index]
        equations.append(linear(residual))
    dirac = sp.zeros(2, 1)
    for flat_index in range(3):
        projected = sp.zeros(2, 1)
        for a_index in range(6):
            projected += v_up[a_index, flat_index] * covariant[a_index]
        dirac += gamma_flat[flat_index] * projected
    equations.append(linear(dirac + spinor / (root2 * l)))

    print("\n=== General-profile first-subleading system ===")
    print("leading vacuum residuals vanish =", all(eq.subs(u, 0) == sp.zeros(2, 1) for eq in equations))
    for index, equation in enumerate(equations):
        coefficient = equation.applyfunc(lambda value: sp.factor(sp.diff(value, u).subs(u, 0)))
        print("order-u equation", index)
        for component in coefficient:
            print(" ", component)


def build_extremal_plus_system() -> None:
    """Exact reduced equations on the nonsingular L_-=0 branch."""
    l, u, lp, w1 = sp.symbols("l u Lp W1", positive=True, nonzero=True)
    u_y = -2 * u / l
    v, vbar, eta = extremal_plus_vielbeins(u, lp, w1)
    jmetric = dft_metric(3)
    p_lower = sp.simplify(v * eta * v.T)
    pbar_lower = sp.simplify(vbar * (-eta) * vbar.T)
    coset_checks = (
        sp.simplify(p_lower + pbar_lower - jmetric) == sp.zeros(6),
        sp.simplify((p_lower * jmetric) ** 2 - p_lower * jmetric) == sp.zeros(6),
    )
    derivative_v_y = sp.diff(v, u) * u_y
    derivative_p_y = sp.diff(p_lower, u) * u_y
    derivative_d = sp.zeros(6, 1)
    derivative_d[5] = -1 / l
    gamma = semi_covariant_connection(
        p_lower, pbar_lower, jmetric, {5: derivative_p_y}, derivative_d, 3
    )
    phi = spin_connection(v, eta, jmetric, gamma, {5: derivative_v_y})

    root2 = sp.sqrt(2)
    gamma_flat = [
        sp.Matrix([[0, 0], [root2, 0]]),
        sp.Matrix([[0, -root2], [0, 0]]),
        sp.diag(1, -1),
    ]
    gamma_bivectors = [
        [
            (gamma_flat[left] * gamma_flat[right]
             - gamma_flat[right] * gamma_flat[left]) / 2
            for right in range(3)
        ]
        for left in range(3)
    ]
    omega = []
    for a_index in range(6):
        matrix = sp.zeros(2)
        for left in range(3):
            for right in range(3):
                matrix += phi[a_index][left, right] * gamma_bivectors[left][right] / 4
        omega.append(sp.simplify(matrix))

    v_up = sp.simplify(jmetric * v * eta)
    vbar_up = sp.simplify(jmetric * vbar * (-eta))
    e0, e1 = sp.symbols("e0 e1")
    ep0, ep1 = sp.symbols("ep0 ep1")
    em0, em1 = sp.symbols("em0 em1")
    ey0, ey1 = sp.symbols("ey0 ey1")
    variables = [ep0, ep1, em0, em1, ey0, ey1, e0, e1]
    spinor = sp.Matrix([e0, e1])
    partials = {
        3: sp.Matrix([ep0, ep1]),
        4: sp.Matrix([em0, em1]),
        5: sp.Matrix([ey0, ey1]),
    }
    covariant = [partials.get(index, sp.zeros(2, 1)) + omega[index] * spinor for index in range(6)]
    equations = []
    labels = []
    for bar_index in range(3):
        residual = sp.zeros(2, 1)
        for a_index in range(6):
            residual += vbar_up[a_index, bar_index] * covariant[a_index]
        equations.append(sp.simplify(residual))
        labels.append(f"gravitino[{bar_index}]")
    dirac = sp.zeros(2, 1)
    for flat_index in range(3):
        projected = sp.zeros(2, 1)
        for a_index in range(6):
            projected += v_up[a_index, flat_index] * covariant[a_index]
        dirac += gamma_flat[flat_index] * projected
    equations.append(sp.simplify(dirac + spinor / (root2 * l)))
    labels.append("shifted dilatino")
    opposite_equations = equations[:3] + [sp.simplify(dirac - spinor / (root2 * l))]

    full = sp.Matrix.vstack(*equations)
    system = coefficient_matrix(full, variables)
    reduced_rref, pivots = system.rref()
    print("\n=== Exact L_-=0 extremal branch ===")
    print("P+Pbar=J and P^2=P =", coset_checks)
    for label, equation in zip(labels, equations):
        print(label)
        for component in equation:
            print(" ", sp.factor(component))
    print("rank =", system.rank(), "pivot columns =", pivots)
    print("independent equations (RREF) =")
    for row in reduced_rref.tolist():
        if any(value != 0 for value in row):
            print(" ", row)

    no_hair = sp.simplify(system.subs(w1, 0))
    no_hair_rref, no_hair_pivots = no_hair.rref()
    print("\nL_-=0, W1=0 rank =", no_hair.rank(), "pivot columns =", no_hair_pivots)
    print("independent equations (RREF) =")
    for row in no_hair_rref.tolist():
        if any(value != 0 for value in row):
            print(" ", row)

    opposite_full = sp.Matrix.vstack(*opposite_equations)
    opposite_system = coefficient_matrix(opposite_full, variables)
    opposite_rref, opposite_pivots = opposite_system.rref()
    print("\nopposite fermion-shift sign rank =", opposite_system.rank(), "pivot columns =", opposite_pivots)
    print("independent equations (RREF) =")
    for row in opposite_rref.tolist():
        if any(value != 0 for value in row):
            print(" ", row)

    opposite_no_hair = sp.simplify(opposite_system.subs(w1, 0))
    opposite_no_hair_rref, opposite_no_hair_pivots = opposite_no_hair.rref()
    print(
        "\nopposite sign with W1=0 rank =",
        opposite_no_hair.rank(),
        "pivot columns =",
        opposite_no_hair_pivots,
    )
    print("independent equations (RREF) =")
    for row in opposite_no_hair_rref.tolist():
        if any(value != 0 for value in row):
            print(" ", row)


def find_opposite_vacuum_projector() -> None:
    """Find the full D=10 internal projector for the opposite reduced shift."""
    l = sp.symbols("l", positive=True, nonzero=True)
    theta = sp.symbols("theta", positive=True)
    v, vbar, eta, jmetric, _, _, _, phi = build_ten_dimensional(l, theta)
    theta_point = sp.pi / 6
    substitutions = {l: 1, theta: theta_point}
    v = sp.simplify(v.subs(substitutions))
    vbar = sp.simplify(vbar.subs(substitutions))
    phi = [sp.simplify(matrix.subs(substitutions)) for matrix in phi]
    v_up = jmetric * v * eta
    vbar_up = jmetric * vbar * (-eta)
    gammas, _ = clifford_matrices()
    gamma_bivectors = [
        [
            (gammas[left] * gammas[right] - gammas[right] * gammas[left]) / 2
            for right in range(10)
        ]
        for left in range(10)
    ]
    omega = []
    for a_index in range(20):
        matrix = sp.zeros(32)
        for left in range(10):
            for right in range(10):
                if phi[a_index][left, right] != 0:
                    matrix += phi[a_index][left, right] * gamma_bivectors[left][right] / 4
        omega.append(sp.simplify(matrix))

    sigma1 = sp.Matrix([[0, 1], [1, 0]])
    sigma2 = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sigma3 = sp.diag(1, -1)
    a_sphere = sp.sin(theta_point) * sigma3 - sp.cos(theta_point) * sigma2
    e_jets = [sp.Matrix([0, 1]), sp.zeros(2, 1)]
    de_jets = [sp.zeros(2, 1), sp.Matrix([0, 1])]
    internal_identity = sp.eye(16)
    chirality_internal = sp.kronecker_product(
        sp.eye(2), sigma3, sp.kronecker_product(sigma3, sigma3)
    )

    print("\n=== Full D=10 opposite-shift vacuum projector search ===")
    for sphere_sign in (1, -1):
        sphere_derivatives = {
            9: sphere_sign * sp.I * sigma1 / 2,
            10: -sphere_sign * sp.I * a_sphere / 2,
            11: sphere_sign * sp.I * a_sphere / 2,
        }
        constraints = []
        for jet in range(2):
            derivative_maps = []
            for a_index in range(20):
                derivative = sp.zeros(32, 16)
                if a_index == 3:
                    derivative += sp.kronecker_product(de_jets[jet], internal_identity)
                if a_index in sphere_derivatives:
                    internal_derivative = sp.kronecker_product(
                        sphere_derivatives[a_index], sp.eye(2), sp.eye(4)
                    )
                    derivative += sp.kronecker_product(
                        e_jets[jet], internal_derivative
                    )
                derivative += omega[a_index] * sp.kronecker_product(
                    e_jets[jet], internal_identity
                )
                derivative_maps.append(sp.simplify(derivative))
            for bar_index in range(10):
                equation = sp.zeros(32, 16)
                for a_index in range(20):
                    if vbar_up[a_index, bar_index] != 0:
                        equation += vbar_up[a_index, bar_index] * derivative_maps[a_index]
                constraints.append(sp.simplify(equation))
            dilatino = sp.zeros(32, 16)
            for flat_index in range(10):
                projected = sp.zeros(32, 16)
                for a_index in range(20):
                    if v_up[a_index, flat_index] != 0:
                        projected += v_up[a_index, flat_index] * derivative_maps[a_index]
                dilatino += gammas[flat_index] * projected
            constraints.append(sp.simplify(dilatino))
        stacked = sp.Matrix.vstack(*constraints)
        stacked = stacked.col_join(chirality_internal - sp.eye(16))
        nullspace = stacked.nullspace()
        print("sphere torsion sign", sphere_sign, "rank/nullity =", stacked.rank(), 16 - stacked.rank())
        if nullspace:
            print("basis rows =")
            for vector in nullspace:
                print(" ", [index for index, value in enumerate(vector) if value != 0])


def full_ten_dimensional_extremal_plus() -> None:
    """Direct D=10 verification of the surviving L_-=0 chiral family."""
    l, u, lp, w0, w1 = sp.symbols("l u Lp W0 W1", positive=True, nonzero=True)
    lpp = sp.symbols("Lp_p")
    w0p, w0m = sp.symbols("W0_p W0_m")
    w1p, w1m = sp.symbols("W1_p W1_m")
    theta = sp.symbols("theta", positive=True)
    theta_point = sp.pi / 6
    v_nr, vbar_nr, eta_nr = extremal_plus_general_vielbeins(u, lp, w0, w1)
    derivative_v_nr_y = sp.diff(v_nr, u) * (-2 * u / l)
    derivative_v_nr_plus = (
        sp.diff(v_nr, lp) * lpp
        + sp.diff(v_nr, w0) * w0p
        + sp.diff(v_nr, w1) * w1p
    )
    derivative_v_nr_minus = (
        sp.diff(v_nr, w0) * w0m + sp.diff(v_nr, w1) * w1m
    )
    p_nr = sp.simplify(v_nr * eta_nr * v_nr.T)
    derivative_p_nr_y = sp.diff(p_nr, u) * (-2 * u / l)
    derivative_p_nr_plus = (
        sp.diff(p_nr, lp) * lpp
        + sp.diff(p_nr, w0) * w0p
        + sp.diff(p_nr, w1) * w1p
    )
    derivative_p_nr_minus = (
        sp.diff(p_nr, w0) * w0m + sp.diff(p_nr, w1) * w1m
    )

    metric_s = sp.diag(l**2, l**2 * sp.cos(theta) ** 2, l**2 * sp.sin(theta) ** 2)
    bfield_s = sp.zeros(3)
    bfield_s[1, 2] = l**2 * sp.cos(theta) ** 2
    bfield_s[2, 1] = -bfield_s[1, 2]
    frame_s = sp.diag(l, l * sp.cos(theta), l * sp.sin(theta))
    v_s, vbar_s, eta_s = riemannian_vielbeins(metric_s, bfield_s, frame_s)
    v_r, vbar_r, eta_r = riemannian_vielbeins(sp.eye(4), sp.zeros(4), sp.eye(4))
    v_symbolic = block_diag(v_nr, v_s, v_r)
    vbar_symbolic = block_diag(vbar_nr, vbar_s, vbar_r)
    eta = block_diag(eta_nr, eta_s, eta_r)
    jmetric = block_diag(sector_dft_metric(3), sector_dft_metric(3), sector_dft_metric(4))
    p_symbolic = sp.simplify(v_symbolic * eta * v_symbolic.T)
    pbar_symbolic = sp.simplify(vbar_symbolic * (-eta) * vbar_symbolic.T)
    derivative_v_y_symbolic = block_diag(derivative_v_nr_y, sp.zeros(6, 3), sp.zeros(8, 4))
    derivative_v_plus_symbolic = block_diag(derivative_v_nr_plus, sp.zeros(6, 3), sp.zeros(8, 4))
    derivative_v_minus_symbolic = block_diag(derivative_v_nr_minus, sp.zeros(6, 3), sp.zeros(8, 4))
    derivative_p_y_symbolic = block_diag(derivative_p_nr_y, sp.zeros(6), sp.zeros(8))
    derivative_p_plus_symbolic = block_diag(derivative_p_nr_plus, sp.zeros(6), sp.zeros(8))
    derivative_p_minus_symbolic = block_diag(derivative_p_nr_minus, sp.zeros(6), sp.zeros(8))
    derivative_v_theta_symbolic = sp.diff(v_symbolic, theta)
    derivative_p_theta_symbolic = sp.diff(p_symbolic, theta)
    substitutions = {theta: theta_point}
    v = sp.simplify(v_symbolic.subs(substitutions))
    vbar = sp.simplify(vbar_symbolic.subs(substitutions))
    p_lower = sp.simplify(p_symbolic.subs(substitutions))
    pbar_lower = sp.simplify(pbar_symbolic.subs(substitutions))
    derivative_v_y = sp.simplify(derivative_v_y_symbolic.subs(substitutions))
    derivative_v_plus = sp.simplify(derivative_v_plus_symbolic.subs(substitutions))
    derivative_v_minus = sp.simplify(derivative_v_minus_symbolic.subs(substitutions))
    derivative_p_y = sp.simplify(derivative_p_y_symbolic.subs(substitutions))
    derivative_p_plus = sp.simplify(derivative_p_plus_symbolic.subs(substitutions))
    derivative_p_minus = sp.simplify(derivative_p_minus_symbolic.subs(substitutions))
    derivative_v_theta = sp.simplify(derivative_v_theta_symbolic.subs(substitutions))
    derivative_p_theta = sp.simplify(derivative_p_theta_symbolic.subs(substitutions))
    derivative_d = sp.zeros(20, 1)
    derivative_d[5] = -1 / l
    derivative_d[9] = -sp.Rational(1, 2) * (
        sp.cot(theta_point) - sp.tan(theta_point)
    )
    gamma = semi_covariant_connection(
        p_lower,
        pbar_lower,
        jmetric,
        {
            3: derivative_p_plus,
            4: derivative_p_minus,
            5: derivative_p_y,
            9: derivative_p_theta,
        },
        derivative_d,
        10,
    )
    phi = spin_connection(
        v,
        eta,
        jmetric,
        gamma,
        {
            3: derivative_v_plus,
            4: derivative_v_minus,
            5: derivative_v_y,
            9: derivative_v_theta,
        },
    )

    gammas, _ = clifford_matrices()
    gamma_bivectors = [
        [
            (gammas[left] * gammas[right] - gammas[right] * gammas[left]) / 2
            for right in range(10)
        ]
        for left in range(10)
    ]
    omega = []
    for a_index in range(20):
        matrix = sp.zeros(32)
        for left in range(10):
            for right in range(10):
                if phi[a_index][left, right] != 0:
                    matrix += phi[a_index][left, right] * gamma_bivectors[left][right] / 4
        omega.append(sp.simplify(matrix))

    sigma1 = sp.Matrix([[0, 1], [1, 0]])
    sigma2 = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sigma3 = sp.diag(1, -1)
    candidate = sp.zeros(16, 4)
    for column, row in enumerate((0, 3, 8, 11)):
        candidate[row, column] = 1
    a_sphere = sp.sin(theta_point) * sigma3 - sp.cos(theta_point) * sigma2
    sphere_derivatives = {
        9: sp.I * sigma1 / 2,
        10: -sp.I * a_sphere / 2,
        11: sp.I * a_sphere / 2,
    }

    e0, e1 = sp.symbols("e0 e1")
    ep0, ep1 = sp.symbols("ep0 ep1")
    em0, em1 = sp.symbols("em0 em1")
    ey0, ey1 = sp.symbols("ey0 ey1")
    variables = [ep0, ep1, em0, em1, ey0, ey1, e0, e1]
    spinor = sp.kronecker_product(sp.Matrix([e0, e1]), candidate)
    partial_nr = {
        3: sp.kronecker_product(sp.Matrix([ep0, ep1]), candidate),
        4: sp.kronecker_product(sp.Matrix([em0, em1]), candidate),
        5: sp.kronecker_product(sp.Matrix([ey0, ey1]), candidate),
    }
    covariant = []
    for a_index in range(20):
        derivative = partial_nr.get(a_index, sp.zeros(32, 4))
        if a_index in sphere_derivatives:
            internal_derivative = sp.kronecker_product(
                sphere_derivatives[a_index], sp.eye(2), sp.eye(4)
            )
            derivative += sp.kronecker_product(
                sp.Matrix([e0, e1]), internal_derivative * candidate
            )
        covariant.append(derivative + omega[a_index] * spinor)

    v_up = sp.simplify(jmetric * v * eta)
    vbar_up = sp.simplify(jmetric * vbar * (-eta))
    residuals = []
    for bar_index in range(10):
        equation = sp.zeros(32, 4)
        for a_index in range(20):
            if vbar_up[a_index, bar_index] != 0:
                equation += vbar_up[a_index, bar_index] * covariant[a_index]
        residuals.append(equation)
    dilatino = sp.zeros(32, 4)
    for flat_index in range(10):
        projected = sp.zeros(32, 4)
        for a_index in range(20):
            if v_up[a_index, flat_index] != 0:
                projected += v_up[a_index, flat_index] * covariant[a_index]
        dilatino += gammas[flat_index] * projected
    residuals.append(dilatino)
    rows = []
    for residual in residuals:
        for entry in list(residual):
            if entry != 0:
                rows.append([sp.factor(sp.diff(entry, variable)) for variable in variables])
    system = sp.Matrix(rows)
    reduced_rref, pivots = system.rref()
    print("\n=== Full D=10 most-general L_-=0 extremal check ===")
    print("nonzero component equations =", system.rows)
    print("rank =", system.rank(), "pivot columns =", pivots)
    print("independent equations (RREF) =")
    for row in reduced_rref.tolist():
        if any(value != 0 for value in row):
            print(" ", row)


def build_extremal_minus_system() -> None:
    """Exact reduced equations on the mirror L_+=0 branch."""
    l, u, lm, w1 = sp.symbols("l u Lm W1", positive=True, nonzero=True)
    v, vbar, eta = extremal_minus_vielbeins(u, lm, w1)
    jmetric = dft_metric(3)
    p_lower = sp.simplify(v * eta * v.T)
    pbar_lower = sp.simplify(vbar * (-eta) * vbar.T)
    derivative_v_y = sp.diff(v, u) * (-2 * u / l)
    derivative_p_y = sp.diff(p_lower, u) * (-2 * u / l)
    derivative_d = sp.zeros(6, 1)
    derivative_d[5] = -1 / l
    gamma = semi_covariant_connection(
        p_lower, pbar_lower, jmetric, {5: derivative_p_y}, derivative_d, 3
    )
    phi = spin_connection(v, eta, jmetric, gamma, {5: derivative_v_y})
    root2 = sp.sqrt(2)
    gamma_flat = [
        sp.Matrix([[0, 0], [root2, 0]]),
        sp.Matrix([[0, -root2], [0, 0]]),
        sp.diag(1, -1),
    ]
    gamma_bivectors = [
        [
            (gamma_flat[left] * gamma_flat[right]
             - gamma_flat[right] * gamma_flat[left]) / 2
            for right in range(3)
        ]
        for left in range(3)
    ]
    omega = []
    for a_index in range(6):
        matrix = sp.zeros(2)
        for left in range(3):
            for right in range(3):
                matrix += phi[a_index][left, right] * gamma_bivectors[left][right] / 4
        omega.append(sp.simplify(matrix))
    v_up = sp.simplify(jmetric * v * eta)
    vbar_up = sp.simplify(jmetric * vbar * (-eta))
    e0, e1 = sp.symbols("e0 e1")
    ep0, ep1 = sp.symbols("ep0 ep1")
    em0, em1 = sp.symbols("em0 em1")
    ey0, ey1 = sp.symbols("ey0 ey1")
    variables = [ep0, ep1, em0, em1, ey0, ey1, e0, e1]
    spinor = sp.Matrix([e0, e1])
    partials = {
        3: sp.Matrix([ep0, ep1]),
        4: sp.Matrix([em0, em1]),
        5: sp.Matrix([ey0, ey1]),
    }
    covariant = [partials.get(index, sp.zeros(2, 1)) + omega[index] * spinor for index in range(6)]
    gravitini = []
    for bar_index in range(3):
        residual = sp.zeros(2, 1)
        for a_index in range(6):
            residual += vbar_up[a_index, bar_index] * covariant[a_index]
        gravitini.append(sp.simplify(residual))
    dirac = sp.zeros(2, 1)
    for flat_index in range(3):
        projected = sp.zeros(2, 1)
        for a_index in range(6):
            projected += v_up[a_index, flat_index] * covariant[a_index]
        dirac += gamma_flat[flat_index] * projected

    print("\n=== Exact L_+=0 mirror extremal branch ===")
    print(
        "P+Pbar=J and P^2=P =",
        sp.simplify(p_lower + pbar_lower - jmetric) == sp.zeros(6),
        sp.simplify((p_lower * jmetric) ** 2 - p_lower * jmetric) == sp.zeros(6),
    )
    for shift_sign in (1, -1):
        equations = gravitini + [sp.simplify(dirac + shift_sign * spinor / (root2 * l))]
        system = coefficient_matrix(sp.Matrix.vstack(*equations), variables)
        reduced_rref, pivots = system.rref()
        print("shift sign", shift_sign, "rank =", system.rank(), "pivot columns =", pivots)
        for row in reduced_rref.tolist():
            if any(value != 0 for value in row):
                print(" ", row)


def build_primed_sector_systems() -> None:
    """Reduced barred-Lorentz (Type-II primed) systems for the limiting branches."""
    l, u, lp, lm, w1 = sp.symbols("l u Lp Lm W1", positive=True, nonzero=True)
    root2 = sp.sqrt(2)
    eta = sp.Matrix([[0, -1, 0], [-1, 0, 0], [0, 0, 1]])
    jmetric = dft_metric(3)
    gamma_flat = [
        sp.Matrix([[0, 0], [root2, 0]]),
        sp.Matrix([[0, -root2], [0, 0]]),
        sp.diag(1, -1),
    ]
    bar_gamma = [sp.I * matrix for matrix in gamma_flat]
    bar_bivectors = [
        [
            (bar_gamma[left] * bar_gamma[right]
             - bar_gamma[right] * bar_gamma[left]) / 2
            for right in range(3)
        ]
        for left in range(3)
    ]

    e0, e1 = sp.symbols("e0 e1")
    ep0, ep1 = sp.symbols("ep0 ep1")
    em0, em1 = sp.symbols("em0 em1")
    ey0, ey1 = sp.symbols("ey0 ey1")
    variables = [ep0, ep1, em0, em1, ey0, ey1, e0, e1]
    spinor = sp.Matrix([e0, e1])
    partials = {
        3: sp.Matrix([ep0, ep1]),
        4: sp.Matrix([em0, em1]),
        5: sp.Matrix([ey0, ey1]),
    }

    backgrounds = []
    v0, vb0, eta0 = exact_nr_vielbeins(sp.Integer(0), sp.Integer(1), sp.Integer(0))
    backgrounds.append(("vacuum", v0, vb0, eta0, sp.zeros(6, 3), sp.zeros(6, 3)))
    vp, vbp, etap = extremal_plus_vielbeins(u, lp, w1)
    backgrounds.append(
        (
            "L_-=0",
            vp,
            vbp,
            etap,
            sp.diff(vp, u) * (-2 * u / l),
            sp.diff(vbp, u) * (-2 * u / l),
        )
    )
    vm, vbm, etam = extremal_minus_vielbeins(u, lm, w1)
    backgrounds.append(
        (
            "L_+=0",
            vm,
            vbm,
            etam,
            sp.diff(vm, u) * (-2 * u / l),
            sp.diff(vbm, u) * (-2 * u / l),
        )
    )

    print("\n=== Primed/barred-Lorentz reduced systems ===")
    for name, v, vbar, eta_branch, derivative_v_y, derivative_vbar_y in backgrounds:
        p_lower = sp.simplify(v * eta_branch * v.T)
        pbar_lower = sp.simplify(vbar * (-eta_branch) * vbar.T)
        derivative_p_y = sp.simplify(
            derivative_v_y * eta_branch * v.T
            + v * eta_branch * derivative_v_y.T
        )
        derivative_d = sp.zeros(6, 1)
        derivative_d[5] = -1 / l
        affine = semi_covariant_connection(
            p_lower, pbar_lower, jmetric, {5: derivative_p_y}, derivative_d, 3
        )
        phibar = spin_connection(
            vbar, -eta_branch, jmetric, affine, {5: derivative_vbar_y}
        )
        omega_bar = []
        for a_index in range(6):
            matrix = sp.zeros(2)
            for left in range(3):
                for right in range(3):
                    matrix += phibar[a_index][left, right] * bar_bivectors[left][right] / 4
            omega_bar.append(sp.simplify(matrix))
        covariant = [
            partials.get(index, sp.zeros(2, 1)) + omega_bar[index] * spinor
            for index in range(6)
        ]
        v_up = sp.simplify(jmetric * v * eta_branch)
        vbar_up = sp.simplify(jmetric * vbar * (-eta_branch))
        gravitini = []
        for flat_index in range(3):
            residual = sp.zeros(2, 1)
            for a_index in range(6):
                residual += v_up[a_index, flat_index] * covariant[a_index]
            gravitini.append(sp.simplify(residual))
        dirac_bar = sp.zeros(2, 1)
        for bar_index in range(3):
            projected = sp.zeros(2, 1)
            for a_index in range(6):
                projected += vbar_up[a_index, bar_index] * covariant[a_index]
            dirac_bar += bar_gamma[bar_index] * projected

        print("branch", name)
        for shift_sign in (1, -1):
            # The extra factor i converts the barred Clifford convention back
            # to a real reduced fermion-shift equation.
            equations = gravitini + [
                sp.simplify(dirac_bar + shift_sign * sp.I * spinor / (root2 * l))
            ]
            system = coefficient_matrix(sp.Matrix.vstack(*equations), variables)
            reduced_rref, pivots = system.rref()
            print(" shift", shift_sign, "rank", system.rank(), "pivots", pivots)
            for row in reduced_rref.tolist():
                if any(value != 0 for value in row):
                    print("  ", row)


def build_extremal_plus_variable_hair_system() -> None:
    """Exact L_-=0 equations with arbitrary W_1(x+,x-) first derivatives."""
    l, u, lp, w1 = sp.symbols("l u Lp W1", positive=True, nonzero=True)
    w1p, w1m = sp.symbols("W1_p W1_m")
    v, vbar, eta = extremal_plus_vielbeins(u, lp, w1)
    jmetric = dft_metric(3)
    p_lower = sp.simplify(v * eta * v.T)
    pbar_lower = sp.simplify(vbar * (-eta) * vbar.T)
    derivative_v_plus = sp.diff(v, w1) * w1p
    derivative_v_minus = sp.diff(v, w1) * w1m
    derivative_v_y = sp.diff(v, u) * (-2 * u / l)
    derivative_p_plus = sp.diff(p_lower, w1) * w1p
    derivative_p_minus = sp.diff(p_lower, w1) * w1m
    derivative_p_y = sp.diff(p_lower, u) * (-2 * u / l)
    derivative_d = sp.zeros(6, 1)
    derivative_d[5] = -1 / l
    affine = semi_covariant_connection(
        p_lower,
        pbar_lower,
        jmetric,
        {3: derivative_p_plus, 4: derivative_p_minus, 5: derivative_p_y},
        derivative_d,
        3,
    )
    phi = spin_connection(
        v,
        eta,
        jmetric,
        affine,
        {3: derivative_v_plus, 4: derivative_v_minus, 5: derivative_v_y},
    )
    root2 = sp.sqrt(2)
    gamma_flat = [
        sp.Matrix([[0, 0], [root2, 0]]),
        sp.Matrix([[0, -root2], [0, 0]]),
        sp.diag(1, -1),
    ]
    gamma_bivectors = [
        [
            (gamma_flat[left] * gamma_flat[right]
             - gamma_flat[right] * gamma_flat[left]) / 2
            for right in range(3)
        ]
        for left in range(3)
    ]
    omega = []
    for a_index in range(6):
        matrix = sp.zeros(2)
        for left in range(3):
            for right in range(3):
                matrix += phi[a_index][left, right] * gamma_bivectors[left][right] / 4
        omega.append(sp.simplify(matrix))

    e0, e1 = sp.symbols("e0 e1")
    ep0, ep1 = sp.symbols("ep0 ep1")
    em0, em1 = sp.symbols("em0 em1")
    ey0, ey1 = sp.symbols("ey0 ey1")
    variables = [ep0, ep1, em0, em1, ey0, ey1, e0, e1]
    spinor = sp.Matrix([e0, e1])
    partials = {
        3: sp.Matrix([ep0, ep1]),
        4: sp.Matrix([em0, em1]),
        5: sp.Matrix([ey0, ey1]),
    }
    covariant = [partials.get(index, sp.zeros(2, 1)) + omega[index] * spinor for index in range(6)]
    v_up = sp.simplify(jmetric * v * eta)
    vbar_up = sp.simplify(jmetric * vbar * (-eta))
    gravitini = []
    for bar_index in range(3):
        residual = sp.zeros(2, 1)
        for a_index in range(6):
            residual += vbar_up[a_index, bar_index] * covariant[a_index]
        gravitini.append(sp.simplify(residual))
    dirac = sp.zeros(2, 1)
    for flat_index in range(3):
        projected = sp.zeros(2, 1)
        for a_index in range(6):
            projected += v_up[a_index, flat_index] * covariant[a_index]
        dirac += gamma_flat[flat_index] * projected

    # The surviving full-D=10 projector corresponds to the opposite reduced
    # fermion-shift sign.
    equations = gravitini + [sp.simplify(dirac - spinor / (root2 * l))]
    system = coefficient_matrix(sp.Matrix.vstack(*equations), variables)
    reduced_rref, pivots = system.rref()
    print("\n=== L_-=0 with arbitrary W1(x+,x-) ===")
    print("rank =", system.rank(), "pivot columns =", pivots)
    print("independent equations (RREF) =")
    for row in reduced_rref.tolist():
        if any(value != 0 for value in row):
            print(" ", row)
    candidate = {
        e0: 0,
        ep0: 0,
        em0: 0,
        ey0: 0,
        em1: 0,
        ey1: 0,
    }
    candidate_residual = [
        sp.factor(entry.subs(candidate))
        for equation in equations
        for entry in equation
        if sp.factor(entry.subs(candidate)) != 0
    ]
    print("residuals on E=(0,F(x+)) =", candidate_residual)
    if candidate_residual:
        raise AssertionError("Kim-frame variable-hair Killing candidate failed")


def build_most_general_extremal_plus_system() -> None:
    """L_-=0 with L_+(x+), W_0(x+,x-), and W_1(x+,x-) arbitrary."""
    l, u, lp, w0, w1 = sp.symbols("l u Lp W0 W1", positive=True, nonzero=True)
    lpp = sp.symbols("Lp_p")
    w0p, w0m, w1p, w1m = sp.symbols("W0_p W0_m W1_p W1_m")
    v, vbar, eta = extremal_plus_general_vielbeins(u, lp, w0, w1)
    jmetric = dft_metric(3)
    p_lower = sp.simplify(v * eta * v.T)
    pbar_lower = sp.simplify(vbar * (-eta) * vbar.T)
    derivative_v_plus = (
        sp.diff(v, lp) * lpp
        + sp.diff(v, w0) * w0p
        + sp.diff(v, w1) * w1p
    )
    derivative_v_minus = sp.diff(v, w0) * w0m + sp.diff(v, w1) * w1m
    derivative_v_y = sp.diff(v, u) * (-2 * u / l)
    derivative_p_plus = (
        sp.diff(p_lower, lp) * lpp
        + sp.diff(p_lower, w0) * w0p
        + sp.diff(p_lower, w1) * w1p
    )
    derivative_p_minus = (
        sp.diff(p_lower, w0) * w0m + sp.diff(p_lower, w1) * w1m
    )
    derivative_p_y = sp.diff(p_lower, u) * (-2 * u / l)
    derivative_d = sp.zeros(6, 1)
    derivative_d[5] = -1 / l
    affine = semi_covariant_connection(
        p_lower,
        pbar_lower,
        jmetric,
        {3: derivative_p_plus, 4: derivative_p_minus, 5: derivative_p_y},
        derivative_d,
        3,
    )
    phi = spin_connection(
        v,
        eta,
        jmetric,
        affine,
        {3: derivative_v_plus, 4: derivative_v_minus, 5: derivative_v_y},
    )
    root2 = sp.sqrt(2)
    gamma_flat = [
        sp.Matrix([[0, 0], [root2, 0]]),
        sp.Matrix([[0, -root2], [0, 0]]),
        sp.diag(1, -1),
    ]
    gamma_bivectors = [
        [
            (gamma_flat[left] * gamma_flat[right]
             - gamma_flat[right] * gamma_flat[left]) / 2
            for right in range(3)
        ]
        for left in range(3)
    ]
    omega = []
    for a_index in range(6):
        matrix = sp.zeros(2)
        for left in range(3):
            for right in range(3):
                matrix += phi[a_index][left, right] * gamma_bivectors[left][right] / 4
        omega.append(sp.simplify(matrix))
    e0, e1 = sp.symbols("e0 e1")
    ep0, ep1 = sp.symbols("ep0 ep1")
    em0, em1 = sp.symbols("em0 em1")
    ey0, ey1 = sp.symbols("ey0 ey1")
    variables = [ep0, ep1, em0, em1, ey0, ey1, e0, e1]
    spinor = sp.Matrix([e0, e1])
    partials = {
        3: sp.Matrix([ep0, ep1]),
        4: sp.Matrix([em0, em1]),
        5: sp.Matrix([ey0, ey1]),
    }
    covariant = [partials.get(index, sp.zeros(2, 1)) + omega[index] * spinor for index in range(6)]
    v_up = sp.simplify(jmetric * v * eta)
    vbar_up = sp.simplify(jmetric * vbar * (-eta))
    gravitini = []
    for bar_index in range(3):
        residual = sp.zeros(2, 1)
        for a_index in range(6):
            residual += vbar_up[a_index, bar_index] * covariant[a_index]
        gravitini.append(sp.simplify(residual))
    dirac = sp.zeros(2, 1)
    for flat_index in range(3):
        projected = sp.zeros(2, 1)
        for a_index in range(6):
            projected += v_up[a_index, flat_index] * covariant[a_index]
        dirac += gamma_flat[flat_index] * projected
    equations = gravitini + [sp.simplify(dirac - spinor / (root2 * l))]
    system = coefficient_matrix(sp.Matrix.vstack(*equations), variables)
    reduced_rref, pivots = system.rref()
    print("\n=== Most general L_-=0 extremal branch ===")
    print("rank =", system.rank(), "pivot columns =", pivots)
    for index, equation in enumerate(equations):
        print("equation", index)
        for component in equation:
            print(" ", sp.factor(component))
    print("independent equations (RREF) =")
    for row in reduced_rref.tolist():
        if any(value != 0 for value in row):
            print(" ", row)
    candidate = {
        e0: 0,
        ep0: 0,
        em0: 0,
        ey0: 0,
        em1: 0,
        ey1: 0,
    }
    candidate_residual = [
        sp.factor(entry.subs(candidate))
        for equation in equations
        for entry in equation
        if sp.factor(entry.subs(candidate)) != 0
    ]
    print("residuals on E=(0,F(x+)) =", candidate_residual)
    if candidate_residual:
        raise AssertionError("Kim-frame general one-sided Killing candidate failed")


def check_lightcone_exchange() -> None:
    """Find the O(3,3) parity map between the two one-sided limits, if any."""
    u, lstate, wstate = sp.symbols("u L W", nonzero=True)
    vp, vbp, eta = extremal_plus_vielbeins(u, lstate, wstate)
    vm, vbm, _ = extremal_minus_vielbeins(u, lstate, wstate)
    hp = sp.simplify(vp * eta * vp.T - vbp * (-eta) * vbp.T)
    hm = sp.simplify(vm * eta * vm.T - vbm * (-eta) * vbm.T)
    jmetric = dft_metric(3)
    swap = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
    matches = []
    for sign_plus in (1, -1):
        for sign_minus in (1, -1):
            for sign_y in (1, -1):
                signs = sp.diag(sign_plus, sign_minus, sign_y)
                physical = signs * swap
                odd = sp.diag(physical, physical)
                if sp.simplify(odd.T * jmetric * odd - jmetric) != sp.zeros(6):
                    continue
                for wsign in (1, -1):
                    transformed = sp.simplify(
                        odd * hp.subs(wstate, wsign * wstate) * odd.T - hm
                    )
                    if transformed == sp.zeros(6):
                        matches.append((sign_plus, sign_minus, sign_y, wsign, odd))
    print("\n=== Lightcone-exchange O(3,3) search ===")
    print("number of exact maps =", len(matches))
    for sign_plus, sign_minus, sign_y, wsign, odd in matches:
        print(
            "coordinate signs (+,-,y), W sign =",
            (sign_plus, sign_minus, sign_y),
            wsign,
        )
        print(odd)

    signed_permutation_matches = []
    longitudinal = (0, 1, 3, 4)
    radial = (2, 5)
    for perm_long in itertools.permutations(longitudinal):
        for signs_long in itertools.product((1, -1), repeat=4):
            for perm_radial in itertools.permutations(radial):
                for signs_radial in itertools.product((1, -1), repeat=2):
                    odd = sp.zeros(6)
                    for source, target, sign in zip(longitudinal, perm_long, signs_long):
                        odd[target, source] = sign
                    for source, target, sign in zip(radial, perm_radial, signs_radial):
                        odd[target, source] = sign
                    if odd.T * jmetric * odd != jmetric:
                        continue
                    for wsign in (1, -1):
                        transformed = sp.simplify(
                            odd * hp.subs(wstate, wsign * wstate) * odd.T - hm
                        )
                        if transformed == sp.zeros(6):
                            signed_permutation_matches.append((wsign, odd))
    print("number of signed-permutation O(3,3) maps =", len(signed_permutation_matches))
    for wsign, odd in signed_permutation_matches[:8]:
        print("W sign =", wsign)
        print(odd)


if __name__ == "__main__":
    selections = set(sys.argv[1:])
    if not selections or "reduced" in selections:
        build_reduced_system()
        build_hair_only_system()
    if not selections or "full-hair" in selections:
        full_ten_dimensional_hair_only()
    if not selections or "full-nonextremal" in selections:
        full_ten_dimensional_nonextremal_point()
    if not selections or "asymptotic" in selections:
        asymptotic_general_profile()
    if not selections or "extremal" in selections:
        build_extremal_plus_system()
    if not selections or "opposite-projector" in selections:
        find_opposite_vacuum_projector()
    if not selections or "full-extremal" in selections:
        full_ten_dimensional_extremal_plus()
    if not selections or "extremal-mirror" in selections:
        build_extremal_minus_system()
    if not selections or "primed" in selections:
        build_primed_sector_systems()
    if not selections or "variable-hair" in selections:
        build_extremal_plus_variable_hair_system()
    if not selections or "most-general" in selections:
        build_most_general_extremal_plus_system()
    if not selections or "parity" in selections:
        check_lightcone_exchange()
