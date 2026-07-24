"""Verify the ten-dimensional DFT Killing-spinor data of the NR3 x S3 x R4 vacuum.

The script constructs the block-diagonal double vielbein used in
NR_Holography.tex, builds the torsionless semi-covariant DFT connection in
dimension D, and obtains the unbarred spin connection

    Phi_Apq = V^B_p (partial_A V_Bq + Gamma_AB^C V_Cq).

It first calibrates the implementation on the three-dimensional NR vacuum and
then repeats the calculation for the full ten-dimensional product.  The output
is deliberately algebraic: it reports projector identities, mixed spin-
connection components, and the Clifford operators entering the gravitino and
dilatino variations.
"""

from __future__ import annotations

import sympy as sp


def block_diag(*matrices: sp.Matrix) -> sp.Matrix:
    return sp.diag(*matrices)


def dft_metric(dimension: int) -> sp.Matrix:
    eye = sp.eye(dimension)
    zero = sp.zeros(dimension)
    return zero.row_join(eye).col_join(eye.row_join(zero))


def sector_dft_metric(dimension: int) -> sp.Matrix:
    """DFT metric for a sector ordered as (dual coordinates; coordinates)."""
    return dft_metric(dimension)


def nr_vielbeins() -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    root2 = sp.sqrt(2)
    v = sp.Matrix(
        [
            [-1, 0, 0],
            [0, 0, 0],
            [0, 0, 1 / root2],
            [0, 1, 0],
            [0, 0, 0],
            [0, 0, 1 / root2],
        ]
    )
    vbar = sp.Matrix(
        [
            [0, 0, 0],
            [1, 0, 0],
            [0, 0, 1 / root2],
            [0, 0, 0],
            [0, 1, 0],
            [0, 0, -1 / root2],
        ]
    )
    eta = sp.Matrix([[0, -1, 0], [-1, 0, 0], [0, 0, 1]])
    return v, vbar, eta


def riemannian_vielbeins(
    metric: sp.Matrix, bfield: sp.Matrix, frame: sp.Matrix
) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Canonical Riemannian double vielbeins in (dual; physical) ordering."""
    inv_transpose = frame.inv().T
    root2 = sp.sqrt(2)
    v = inv_transpose.col_join((metric + bfield) * inv_transpose) / root2
    vbar = inv_transpose.col_join((bfield - metric) * inv_transpose) / root2
    return v, vbar, sp.eye(metric.rows)


def permute_global_to_sector_order(matrix: sp.Matrix, sizes: list[int]) -> sp.Matrix:
    """Convert global (all dual; all physical) order to sector-block order."""
    total = sum(sizes)
    dual_offsets = []
    physical_offsets = []
    running = 0
    for size in sizes:
        dual_offsets.append(running)
        physical_offsets.append(total + running)
        running += size
    order: list[int] = []
    for dual, physical, size in zip(dual_offsets, physical_offsets, sizes):
        order.extend(range(dual, dual + size))
        order.extend(range(physical, physical + size))
    return matrix.extract(order, order)


def semi_covariant_connection(
    p_lower: sp.Matrix,
    pbar_lower: sp.Matrix,
    jmetric: sp.Matrix,
    derivatives_p: dict[int, sp.Matrix],
    derivative_d: sp.Matrix,
    physical_dimension: int,
) -> list[sp.Matrix]:
    """Return Gamma_CAB as a list of matrices in the last two indices.

    Antisymmetrization brackets have unit weight, e.g. X_[AB]=(X_AB-X_BA)/2.
    The formula is Eq. (13) of arXiv:1112.0069 with fermionic torsion set to
    zero.
    """
    doubled = p_lower.rows
    p_mixed = p_lower * jmetric
    pbar_mixed = pbar_lower * jmetric
    pbar_up_down = jmetric * pbar_lower

    # T_ED = (P partial^E P Pbar)_ED.  The raised derivative direction is
    # related to a lower derivative direction by the O(D,D) metric.
    tmat = sp.zeros(doubled)
    for raised_e in range(doubled):
        derivative = sp.zeros(doubled)
        for lower_f, dp in derivatives_p.items():
            if jmetric[raised_e, lower_f] != 0:
                derivative += jmetric[raised_e, lower_f] * dp
        if derivative == sp.zeros(doubled):
            continue
        row = p_mixed[raised_e, :] * derivative * pbar_up_down
        for d_index in range(doubled):
            tmat[raised_e, d_index] = row[0, d_index]
    t_antisym = (tmat - tmat.T) / 2
    trace_vector = derivative_d + sp.Matrix(
        [sum(t_antisym[e, d] for e in range(doubled)) for d in range(doubled)]
    )

    gamma: list[sp.Matrix] = []
    coefficient = -sp.Rational(4, physical_dimension - 1)
    for c_index in range(doubled):
        result = sp.zeros(doubled)

        # 2 (P partial_C P Pbar)_[AB].
        if c_index in derivatives_p:
            mixed = p_mixed * derivatives_p[c_index] * pbar_up_down
            result = sp.MutableDenseMatrix(result + mixed - mixed.T)

        # 2 (Pbar_[A^D Pbar_B]^E - P_[A^D P_B]^E) partial_D P_EC.
        for d_index, dp in derivatives_p.items():
            for a_index in range(doubled):
                for b_index in range(doubled):
                    term = 0
                    for e_index in range(doubled):
                        term += (
                            pbar_mixed[a_index, d_index]
                            * pbar_mixed[b_index, e_index]
                            - pbar_mixed[b_index, d_index]
                            * pbar_mixed[a_index, e_index]
                            - p_mixed[a_index, d_index]
                            * p_mixed[b_index, e_index]
                            + p_mixed[b_index, d_index]
                            * p_mixed[a_index, e_index]
                        ) * dp[e_index, c_index]
                    result[a_index, b_index] += term

        # -4/(D-1) (Pbar_C[A Pbar_B]^D + P_C[A P_B]^D) X_D.
        for a_index in range(doubled):
            for b_index in range(doubled):
                trace_term = 0
                for d_index in range(doubled):
                    trace_term += sp.Rational(1, 2) * (
                        pbar_lower[c_index, a_index]
                        * pbar_mixed[b_index, d_index]
                        - pbar_lower[c_index, b_index]
                        * pbar_mixed[a_index, d_index]
                        + p_lower[c_index, a_index]
                        * p_mixed[b_index, d_index]
                        - p_lower[c_index, b_index]
                        * p_mixed[a_index, d_index]
                    ) * trace_vector[d_index]
                result[a_index, b_index] += coefficient * trace_term
        gamma.append(sp.simplify(result))
    return gamma


def spin_connection(
    v: sp.Matrix,
    eta: sp.Matrix,
    jmetric: sp.Matrix,
    gamma: list[sp.Matrix],
    derivatives_v: dict[int, sp.Matrix],
) -> list[sp.Matrix]:
    """Return Phi_Apq with both Lorentz indices lowered."""
    doubled = v.rows
    v_lower_flat = v * eta
    v_up_lower_flat = jmetric * v * eta
    result: list[sp.Matrix] = []
    for a_index in range(doubled):
        derivative = derivatives_v.get(a_index, sp.zeros(*v.shape)) * eta
        covariant = derivative + gamma[a_index] * jmetric * v_lower_flat
        phi = v_up_lower_flat.T * covariant
        result.append(sp.simplify((phi - phi.T) / 2))
    return result


def nonzero_entries(matrix: sp.Matrix) -> list[tuple[int, int, sp.Expr]]:
    entries = []
    for row in range(matrix.rows):
        for column in range(matrix.cols):
            value = sp.simplify(matrix[row, column])
            if value != 0:
                entries.append((row, column, value))
    return entries


def kron(*matrices: sp.Matrix) -> sp.Matrix:
    result = sp.Matrix([[1]])
    for matrix in matrices:
        result = sp.kronecker_product(result, matrix)
    return result


def clifford_matrices() -> tuple[list[sp.Matrix], sp.Matrix]:
    """A 32-component Cl(1,9) representation adapted to 3+3+4.

    Tensor ordering is NR2 x S3_2 x auxiliary2 x R4_4.  The auxiliary
    Pauli matrices make the two odd-dimensional Clifford factors anticommute.
    """
    root2 = sp.sqrt(2)
    imaginary = sp.I
    sigma1 = sp.Matrix([[0, 1], [1, 0]])
    sigma2 = sp.Matrix([[0, -imaginary], [imaginary, 0]])
    sigma3 = sp.diag(1, -1)
    identity2 = sp.eye(2)
    identity4 = sp.eye(4)

    gamma_nr = [
        sp.Matrix([[0, 0], [root2, 0]]),
        sp.Matrix([[0, -root2], [0, 0]]),
        sigma3,
    ]
    tau = [sigma1, sigma2, sigma3]
    rho = [
        kron(sigma1, identity2),
        kron(sigma2, identity2),
        kron(sigma3, sigma1),
        kron(sigma3, sigma2),
    ]
    gammas = [kron(gamma, identity2, sigma1, identity4) for gamma in gamma_nr]
    gammas += [kron(identity2, gamma, sigma2, identity4) for gamma in tau]
    gammas += [kron(identity2, identity2, sigma3, gamma) for gamma in rho]
    chirality_r4 = kron(sigma3, sigma3)
    chirality_10 = kron(identity2, identity2, sigma3, chirality_r4)
    return gammas, chirality_10


def projected_killing_system(
    v: sp.Matrix,
    vbar: sp.Matrix,
    eta: sp.Matrix,
    jmetric: sp.Matrix,
    phi: list[sp.Matrix],
    theta: sp.Symbol,
    l: sp.Symbol,
    theta_point: sp.Expr,
) -> None:
    """Evaluate the full 10D Killing system at a regular generic S3 point.

    The S3 derivatives are those of a spinor parallel under the +H
    torsionful connection associated with B_{phi1 phi2}=l^2 cos^2(theta).
    The NR profile is E=(l f',f), and f,f',f'' are treated as independent
    jets.  We then solve all gravitino and dilatino equations simultaneously
    on a positive-chirality ten-dimensional spinor.
    """
    gammas, chirality_10 = clifford_matrices()
    identity2 = sp.eye(2)
    identity4 = sp.eye(4)
    sigma1 = sp.Matrix([[0, 1], [1, 0]])
    sigma2 = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sigma3 = sp.diag(1, -1)

    substitutions = {l: 1, theta: theta_point}
    v_point = sp.simplify(v.subs(substitutions))
    vbar_point = sp.simplify(vbar.subs(substitutions))
    eta_point = eta
    bareta_point = -eta_point
    phi_point = [sp.simplify(item.subs(substitutions)) for item in phi]
    v_up_lower = jmetric * v_point * eta_point
    vbar_up_lower = jmetric * vbar_point * bareta_point

    gamma_bivectors: list[list[sp.Matrix]] = []
    for left in range(10):
        row = []
        for right in range(10):
            row.append((gammas[left] * gammas[right] - gammas[right] * gammas[left]) / 2)
        gamma_bivectors.append(row)

    omega: list[sp.Matrix] = []
    for a_index in range(20):
        connection = sp.zeros(32)
        for p_index in range(10):
            for q_index in range(10):
                if phi_point[a_index][p_index, q_index] != 0:
                    connection += (
                        phi_point[a_index][p_index, q_index]
                        * gamma_bivectors[p_index][q_index]
                        / 4
                    )
        omega.append(sp.simplify(connection))

    # S3 torsionful-parallel derivatives at the chosen point.
    a_sphere = sp.sin(theta_point) * sigma3 - sp.cos(theta_point) * sigma2
    sphere_derivatives = {
        9: sp.I * sigma1 / 2,
        10: -sp.I * a_sphere / 2,
        11: sp.I * a_sphere / 2,
    }

    # Coefficients of the independent jets (f,f',f'') in E and partial_+ E.
    e_jets = [sp.Matrix([0, 1]), sp.Matrix([1, 0]), sp.zeros(2, 1)]
    de_jets = [sp.zeros(2, 1), sp.Matrix([0, 1]), sp.Matrix([1, 0])]
    internal_identity = sp.eye(16)

    def derivative_map(a_index: int, jet: int) -> sp.Matrix:
        result = sp.zeros(32, 16)
        if a_index == 3:
            result += sp.kronecker_product(de_jets[jet], internal_identity)
        if a_index in sphere_derivatives:
            internal_derivative = kron(
                sphere_derivatives[a_index], sp.eye(2), identity4
            )
            result += sp.kronecker_product(e_jets[jet], internal_derivative)
        result += omega[a_index] * sp.kronecker_product(
            e_jets[jet], internal_identity
        )
        return sp.simplify(result)

    derivative_maps = [
        [derivative_map(a_index, jet) for a_index in range(20)]
        for jet in range(3)
    ]

    constraints: list[sp.Matrix] = []
    for jet in range(3):
        # Ten gravitino equations D_bar{p} epsilon=0.
        for bar_index in range(10):
            equation = sp.zeros(32, 16)
            for a_index in range(20):
                coefficient = vbar_up_lower[a_index, bar_index]
                if coefficient != 0:
                    equation += coefficient * derivative_maps[jet][a_index]
            constraints.append(sp.simplify(equation))

        # Dilatino equation gamma^p D_p epsilon=0.
        equation = sp.zeros(32, 16)
        for p_index in range(10):
            projected = sp.zeros(32, 16)
            for a_index in range(20):
                coefficient = v_up_lower[a_index, p_index]
                if coefficient != 0:
                    projected += coefficient * derivative_maps[jet][a_index]
            equation += gammas[p_index] * projected
        constraints.append(sp.simplify(equation))

    stacked = constraints[0]
    for equation in constraints[1:]:
        stacked = stacked.col_join(equation)
    unconstrained_rank = stacked.rank()
    print("theta =", theta_point)
    print("Killing equations before Weyl projection: rank/nullity =", unconstrained_rank, 16 - unconstrained_rank)

    # Positive 10D Weyl chirality reduces to an internal condition because the
    # chosen chirality matrix is identity on the NR spinor factor.
    chirality_internal = kron(sp.eye(2), sigma3, kron(sigma3, sigma3))
    stacked = stacked.col_join(chirality_internal - sp.eye(16))
    rank = stacked.rank()
    nullity = 16 - rank
    print("full Killing-system rank on one NR chiral family =", rank)
    print("full Killing-system internal nullity =", nullity)
    if nullity:
        print("basis vectors:")
        for vector in stacked.nullspace():
            print(vector.T)


def prove_symbolic_candidate(
    v: sp.Matrix,
    vbar: sp.Matrix,
    eta: sp.Matrix,
    jmetric: sp.Matrix,
    phi: list[sp.Matrix],
    theta: sp.Symbol,
    l: sp.Symbol,
) -> None:
    """Prove the four-polarization Weyl candidate for symbolic theta and l."""
    gammas, _ = clifford_matrices()
    identity4 = sp.eye(4)
    sigma1 = sp.Matrix([[0, 1], [1, 0]])
    sigma2 = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sigma3 = sp.diag(1, -1)
    gamma_bivectors = [
        [
            (gammas[p] * gammas[q] - gammas[q] * gammas[p]) / 2
            for q in range(10)
        ]
        for p in range(10)
    ]
    # Rationalize all trigonometric functions with u=tan(theta).  The physical
    # coordinate patch has 0<theta<pi/2, so this is an exact parametrization
    # and turns the final identity test into rational algebra.
    u = sp.symbols("u", positive=True)
    root = sp.sqrt(1 + u**2)

    def trig_to_u(value: sp.Expr) -> sp.Expr:
        expanded = sp.expand_trig(value)
        replaced = expanded.subs(
            {
                sp.tan(theta): u,
                sp.cot(theta): 1 / u,
                sp.sin(theta): u / root,
                sp.cos(theta): 1 / root,
            },
            simultaneous=True,
        )
        return sp.factor(sp.cancel(replaced))

    v_rational = v.applyfunc(trig_to_u)
    vbar_rational = vbar.applyfunc(trig_to_u)
    phi_rational = [matrix.applyfunc(trig_to_u) for matrix in phi]
    v_up_lower = jmetric * v_rational * eta
    vbar_up_lower = jmetric * vbar_rational * (-eta)

    # Internal basis: arbitrary S3 spinor, auxiliary sigma3=-1, and negative
    # R4 chirality.  In the chosen tensor basis these are positions 5,6,13,14.
    candidate = sp.zeros(16, 4)
    for column, row in enumerate((5, 6, 13, 14)):
        candidate[row, column] = 1

    a_sphere = (u * sigma3 - sigma2) / root
    sphere_derivatives = {
        9: sp.I * sigma1 / 2,
        10: -sp.I * a_sphere / 2,
        11: sp.I * a_sphere / 2,
    }
    e_jets = [sp.Matrix([0, 1]), sp.Matrix([l, 0]), sp.zeros(2, 1)]
    de_jets = [sp.zeros(2, 1), sp.Matrix([0, 1]), sp.Matrix([l, 0])]

    def simplify_zero(matrix: sp.Matrix) -> sp.Matrix:
        return matrix.applyfunc(lambda value: sp.factor(sp.cancel(value)))

    nonzero_residuals = []
    for jet in range(3):
        spinor = sp.kronecker_product(e_jets[jet], candidate)
        derivative_maps = []
        for a_index in range(20):
            derivative = sp.zeros(32, 4)
            if a_index == 3:
                derivative += sp.kronecker_product(de_jets[jet], candidate)
            if a_index in sphere_derivatives:
                internal_derivative = kron(
                    sphere_derivatives[a_index], sp.eye(2), identity4
                )
                derivative += sp.kronecker_product(
                    e_jets[jet], internal_derivative * candidate
                )
            connection_action = sp.zeros(32, 4)
            for p_index in range(10):
                for q_index in range(10):
                    coefficient = phi_rational[a_index][p_index, q_index]
                    if coefficient != 0:
                        connection_action += (
                            coefficient
                            * gamma_bivectors[p_index][q_index]
                            * spinor
                            / 4
                        )
            derivative_maps.append(derivative + connection_action)

        for bar_index in range(10):
            equation = sp.zeros(32, 4)
            for a_index in range(20):
                coefficient = vbar_up_lower[a_index, bar_index]
                if coefficient != 0:
                    equation += coefficient * derivative_maps[a_index]
            equation = simplify_zero(equation)
            if equation != sp.zeros(32, 4):
                nonzero_residuals.append(("gravitino", jet, bar_index, nonzero_entries(equation)))

        equation = sp.zeros(32, 4)
        for p_index in range(10):
            projected = sp.zeros(32, 4)
            for a_index in range(20):
                coefficient = v_up_lower[a_index, p_index]
                if coefficient != 0:
                    projected += coefficient * derivative_maps[a_index]
            equation += gammas[p_index] * projected
        equation = simplify_zero(equation)
        if equation != sp.zeros(32, 4):
            nonzero_residuals.append(("dilatino", jet, nonzero_entries(equation)))

    print("symbolic all-theta/l candidate residual count =", len(nonzero_residuals))
    for residual in nonzero_residuals[:10]:
        print(" ", residual)


def build_three_dimensional(l: sp.Symbol):
    v, vbar, eta = nr_vielbeins()
    jmetric = sector_dft_metric(3)
    p_lower = sp.simplify(v * eta * v.T)
    pbar_lower = sp.simplify(vbar * (-eta) * vbar.T)
    derivative_d = sp.zeros(6, 1)
    derivative_d[5] = -1 / l
    gamma = semi_covariant_connection(
        p_lower, pbar_lower, jmetric, {}, derivative_d, 3
    )
    phi = spin_connection(v, eta, jmetric, gamma, {})
    return v, vbar, eta, jmetric, p_lower, pbar_lower, gamma, phi


def check_three_dimensional_profile(
    v: sp.Matrix,
    vbar: sp.Matrix,
    eta: sp.Matrix,
    jmetric: sp.Matrix,
    phi: list[sp.Matrix],
    l: sp.Symbol,
) -> None:
    """Check D_bar p E=0 and gamma^p D_p E=-E/(sqrt(2)l)."""
    root2 = sp.sqrt(2)
    gamma_flat = [
        sp.Matrix([[0, 0], [root2, 0]]),
        sp.Matrix([[0, -root2], [0, 0]]),
        sp.diag(1, -1),
    ]
    gamma_bivector = [
        [
            (gamma_flat[p] * gamma_flat[q] - gamma_flat[q] * gamma_flat[p]) / 2
            for q in range(3)
        ]
        for p in range(3)
    ]
    omega = []
    for a_index in range(6):
        connection = sp.zeros(2)
        for p_index in range(3):
            for q_index in range(3):
                connection += (
                    phi[a_index][p_index, q_index]
                    * gamma_bivector[p_index][q_index]
                    / 4
                )
        omega.append(sp.simplify(connection))
    v_up_lower = jmetric * v * eta
    vbar_up_lower = jmetric * vbar * (-eta)
    e_jets = [sp.Matrix([0, 1]), sp.Matrix([l, 0]), sp.zeros(2, 1)]
    de_jets = [sp.zeros(2, 1), sp.Matrix([0, 1]), sp.Matrix([l, 0])]
    success = True
    for jet in range(3):
        derivative_maps = []
        for a_index in range(6):
            derivative = de_jets[jet] if a_index == 3 else sp.zeros(2, 1)
            derivative_maps.append(sp.simplify(derivative + omega[a_index] * e_jets[jet]))
        for bar_index in range(3):
            residual = sp.zeros(2, 1)
            for a_index in range(6):
                residual += vbar_up_lower[a_index, bar_index] * derivative_maps[a_index]
            success = success and sp.simplify(residual) == sp.zeros(2, 1)
        dirac = sp.zeros(2, 1)
        for p_index in range(3):
            projected = sp.zeros(2, 1)
            for a_index in range(6):
                projected += v_up_lower[a_index, p_index] * derivative_maps[a_index]
            dirac += gamma_flat[p_index] * projected
        residual = sp.simplify(dirac + e_jets[jet] / (root2 * l))
        success = success and residual == sp.zeros(2, 1)
    print("E=(l f',f) satisfies the quoted D=3 equations =", success)


def build_ten_dimensional(l: sp.Symbol, theta: sp.Symbol):
    v_nr, vbar_nr, eta_nr = nr_vielbeins()

    metric_s = sp.diag(l**2, l**2 * sp.cos(theta) ** 2, l**2 * sp.sin(theta) ** 2)
    bfield_s = sp.zeros(3)
    bfield_s[1, 2] = l**2 * sp.cos(theta) ** 2
    bfield_s[2, 1] = -bfield_s[1, 2]
    frame_s = sp.diag(l, l * sp.cos(theta), l * sp.sin(theta))
    v_s, vbar_s, eta_s = riemannian_vielbeins(metric_s, bfield_s, frame_s)

    metric_r = sp.eye(4)
    bfield_r = sp.zeros(4)
    frame_r = sp.eye(4)
    v_r, vbar_r, eta_r = riemannian_vielbeins(metric_r, bfield_r, frame_r)

    v = block_diag(v_nr, v_s, v_r)
    vbar = block_diag(vbar_nr, vbar_s, vbar_r)
    eta = block_diag(eta_nr, eta_s, eta_r)
    bareta = -eta
    jmetric = block_diag(sector_dft_metric(3), sector_dft_metric(3), sector_dft_metric(4))
    p_lower = sp.simplify(v * eta * v.T)
    pbar_lower = sp.simplify(vbar * bareta * vbar.T)

    # Physical coordinate indices in sector order: y=5, theta=9.
    derivative_v_theta = sp.diff(v, theta)
    derivative_p_theta = sp.diff(p_lower, theta)
    derivative_d = sp.zeros(20, 1)
    derivative_d[5] = -1 / l
    derivative_d[9] = -sp.Rational(1, 2) * (
        sp.cot(theta) - sp.tan(theta)
    )
    gamma = semi_covariant_connection(
        p_lower,
        pbar_lower,
        jmetric,
        {9: derivative_p_theta},
        derivative_d,
        10,
    )
    phi = spin_connection(v, eta, jmetric, gamma, {9: derivative_v_theta})
    return v, vbar, eta, jmetric, p_lower, pbar_lower, gamma, phi


def main() -> None:
    l = sp.symbols("l", positive=True, nonzero=True)
    theta = sp.symbols("theta", positive=True)

    print("=== D=3 calibration ===")
    v3, vb3, eta3, j3, p3, pb3, gamma3, phi3 = build_three_dimensional(l)
    print("P+Pbar-J =", sp.simplify(p3 + pb3 - j3) == sp.zeros(6))
    for a_index, matrix in enumerate(phi3):
        entries = nonzero_entries(matrix)
        if entries:
            print(f"Phi3[{a_index}] = {entries}")
    check_three_dimensional_profile(v3, vb3, eta3, j3, phi3, l)

    print("\n=== D=10 product ===")
    v10, vb10, eta10, j10, p10, pb10, gamma10, phi10 = build_ten_dimensional(l, theta)
    print("P+Pbar-J =", sp.simplify(p10 + pb10 - j10) == sp.zeros(20))
    print("P^2-P =", sp.simplify((p10 * j10) ** 2 - p10 * j10) == sp.zeros(20))

    block_labels = ["NR", "S3", "R4"]
    tangent_ranges = [range(0, 3), range(3, 6), range(6, 10)]
    mixed = []
    for a_index, matrix in enumerate(phi10):
        for left, left_range in enumerate(tangent_ranges):
            for right, right_range in enumerate(tangent_ranges):
                if left >= right:
                    continue
                for p_index in left_range:
                    for q_index in right_range:
                        value = sp.simplify(matrix[p_index, q_index])
                        if value != 0:
                            mixed.append(
                                (a_index, block_labels[left], p_index, block_labels[right], q_index, value)
                            )
    print("mixed Phi components:")
    for entry in mixed:
        print(" ", entry)

    # Compare the NR-NR part of the ten-dimensional connection with D=3.
    print("NR-NR spin-connection blocks in D=10:")
    for a_index, matrix in enumerate(phi10):
        block = matrix[:3, :3]
        entries = nonzero_entries(block)
        if entries:
            print(f"Phi10[{a_index}]_NR = {entries}")

    print("S3-S3 spin-connection blocks in D=10:")
    for a_index, matrix in enumerate(phi10):
        block = matrix[3:6, 3:6]
        entries = nonzero_entries(block)
        if entries:
            print(f"Phi10[{a_index}]_S3 = {entries}")

    print("\n=== Full 10D Killing system at generic S3 points ===")
    projected_killing_system(v10, vb10, eta10, j10, phi10, theta, l, sp.pi / 6)
    projected_killing_system(v10, vb10, eta10, j10, phi10, theta, l, sp.pi / 3)
    print("\n=== Symbolic candidate proof ===")
    prove_symbolic_candidate(v10, vb10, eta10, j10, phi10, theta, l)


if __name__ == "__main__":
    main()

