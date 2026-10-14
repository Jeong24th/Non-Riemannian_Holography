"""Exact local SM5 checks. Python 3.12 / SymPy 1.14; no manuscript input.

Run: python checks/verify_sm5_abelian_rr.py
All functions are arbitrary smooth functions, with only the stated chiral
dependence imposed. Odd jets form a free exterior algebra: a fermion and its
derivatives are independent generators, not commuting symbols or numerical
samples. This verifies the fixed D=2 model, not a holographic identification.
"""
import json
import sys
import sympy as s

x, y = s.symbols('x_plus x_minus')
checks = []


def check(name, residual):
    if isinstance(residual, Exterior):
        ok = not residual.clean().terms
    elif isinstance(residual, (bool, s.logic.boolalg.BooleanAtom)):
        ok = bool(residual)
    else:
        values = list(residual) if isinstance(residual, s.MatrixBase) else [residual]
        ok = all(s.simplify(s.expand(v)) == 0 for v in values)
    checks.append({'check': name, 'passed': ok})
    print(('[PASS] ' if ok else '[FAIL] ') + name)
    if not ok:
        print(repr(residual))


class Exterior:
    """Sparse Grassmann polynomials, with commuting SymPy coefficients."""
    def __init__(self, terms=None):
        self.terms = terms or {}

    @staticmethod
    def scalar(value):
        return Exterior({(): s.sympify(value)})

    @staticmethod
    def odd(name, p=0, m=0):
        return Exterior({((name, p, m),): s.Integer(1)})

    def clean(self):
        return Exterior({k: s.expand(v) for k, v in self.terms.items()
                         if s.expand(v) != 0})

    def __add__(self, other):
        if not isinstance(other, Exterior):
            other = Exterior.scalar(other)
        out = dict(self.terms)
        for k, v in other.terms.items():
            out[k] = out.get(k, 0) + v
        return Exterior(out).clean()

    __radd__ = __add__

    def __neg__(self):
        return Exterior({k: -v for k, v in self.terms.items()})

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if not isinstance(other, Exterior):
            other = Exterior.scalar(other)
        out = {}
        for a, va in self.terms.items():
            for b, vb in other.terms.items():
                if set(a).intersection(b):
                    continue
                sign = (-1)**sum(i > j for i in a for j in b)
                key = tuple(sorted(a+b))
                out[key] = out.get(key, 0) + sign*va*vb
        return Exterior(out).clean()

    __rmul__ = __mul__

    def diff(self, axis):
        out = Exterior()
        for monomial, coeff in self.terms.items():
            out += Exterior({monomial: s.diff(coeff, (x, y)[axis])})
            for i, (name, p, m) in enumerate(monomial):
                # Supersymmetry parameters are arbitrary chiral functions.
                if (name.endswith('P') and axis == 1) or (name.endswith('M') and axis == 0):
                    continue
                term = Exterior.scalar(coeff)
                for j, (n, a, b) in enumerate(monomial):
                    term *= Exterior.odd(n, a + (i == j and axis == 0),
                                         b + (i == j and axis == 1))
                out += term
        return out.clean()

    def on_shell(self):
        return Exterior({k: v for k, v in self.terms.items()
                         if not any((n == 'psi_plus' and m) or
                                    (n == 'psi_minus' and p) for n, p, m in k)})

    def __repr__(self):
        return repr(self.terms)


I = s.eye(2)
J = s.Matrix([[0,0,1,0], [0,0,0,1], [1,0,0,0], [0,1,0,0]])
eta = s.Matrix([[0,-1],[-1,0]])
V = s.Matrix([[1/s.sqrt(2),0],[0,0],[0,-s.sqrt(2)],[0,0]])
Vb = s.Matrix([[0,0],[0,1/s.sqrt(2)],[0,0],[s.sqrt(2),0]])
H = V*eta*V.T - Vb*(-eta)*Vb.T
check('4b: DFT vielbein normalization and orthogonality',
      s.diag(V.T*J*V-eta, Vb.T*J*Vb+eta, V.T*J*Vb))
check('4b: generalized metric constraint', H*J*H-J)
gp = s.sqrt(2)*s.Matrix([[0,1],[0,0]])
gm = s.sqrt(2)*s.Matrix([[0,0],[-1,0]])
gbp, gbm = gp, -gm
G, Gb = s.diag(1,-1), s.diag(-1,1)
C = s.Matrix([[0,1],[1,0]])
for prefix, gammas, metric, chirality in [('unbarred',(gp,gm),eta,G), ('barred',(gbp,gbm),-eta,Gb)]:
    for a in range(2):
        for b in range(2):
            check(f'SM105 {prefix} Clifford {a}{b}', gammas[a]*gammas[b]+gammas[b]*gammas[a]-2*metric[a,b]*I)
    check(f'SM105 {prefix} chirality definition', (gammas[0]+gammas[1])*(gammas[0]-gammas[1])/2-chirality)
    for a in range(2):
        check(f'SM92 {prefix} conjugation {a}', gammas[a].T-C*gammas[a]*C)

def field(name):
    return s.Function(name)(x,y)

B = s.Matrix([[field('c'),field('b_plus')],[field('b_minus'),field('phi')]])
phi = B[1,1]

def D(z, sign=1):
    return (gp*z.diff(x)+sign*G*z.diff(y)*gbm)/s.sqrt(2)

F = D(B)
Fb = C*F.T*C
check('SM94/106 D_plus squared', D(D(B)))
check('SM94/106 D_minus squared', D(D(B,-1),-1))
check('SM107 unrestricted field strength', F-s.Matrix([[B[1,0].diff(x)+B[0,1].diff(y),phi.diff(x)],[-phi.diff(y),0]]))
check('SM108 bosonic action', -s.trace(F*Fb)/2-phi.diff(x)*phi.diff(y))
check('SM96/108 bosonic equation', D(F,-1)+2*phi.diff(x,y)*s.Matrix([[1,0],[0,0]]))
check('SM99 RR invariance and reducibility', D(D(B)))
check('SM110 RR leaves phi unchanged', D(B)[1,1])

# Recompute the d-dependent DFT Christoffel and spin connections at fixed V.
# Gamma_CAB = -4/(D-1) (P_C[A P_B]^D + Pbar_C[A Pbar_B]^D) partial_D d.
P, Pb = (J+H)/2, (J-H)/2
grad = s.Matrix([0,0,*s.symbols('d_plus d_minus')])
Gamma = []
for a in range(4):
    connection = s.zeros(4)
    for projector in (P,Pb):
        row = projector[a,:].T
        projected = projector*J*grad
        connection += -2*(row*projected.T-projected*row.T)
    Gamma.append(connection)
trace = s.Matrix([sum(J[b,c]*Gamma[c][b,a] for b in range(4) for c in range(4)) for a in range(4)])
check('SM93 Christoffel trace', trace+2*grad)
def spin_connection(vielbein, gammas):
    out = []
    for a in range(4):
        flat = (J*vielbein).T*Gamma[a]*J*vielbein
        out.append(sum((flat[b,c]*(gammas[b]*gammas[c]-gammas[c]*gammas[b])/8
                        for b in range(2) for c in range(2)), s.zeros(2)))
    return out
spin, spinb = spin_connection(V,(gp,gm)), spin_connection(Vb,(gbp,gbm))
gA = [sum(((J*V)[a,b]*(gp,gm)[b] for b in range(2)), s.zeros(2)) for a in range(4)]
gbA = [sum(((J*Vb)[a,b]*(gbp,gbm)[b] for b in range(2)), s.zeros(2)) for a in range(4)]
check('SM93-95 unbarred half-density Dirac d cancellation', sum((gA[a]*(grad[a]*I+spin[a]) for a in range(4)),s.zeros(2)))
check('SM93-95 barred half-density Dirac d cancellation', sum((gbA[a]*(grad[a]*I+spinb[a]) for a in range(4)),s.zeros(2)))
for sign in (1,-1):
    check(f'SM94 D_{sign:+d} d cancellation', sum((gA[a]*(grad[a]*B+spin[a]*B-B*spinb[a])+sign*G*(grad[a]*B+spin[a]*B-B*spinb[a])*gbA[a] for a in range(4)),s.zeros(2)))

vp, vm = s.Function('v_plus')(x), s.Function('v_minus')(y)
wp, wm = s.Function('w_plus')(x), s.Function('w_minus')(y)
def kosmann(z, p=vp, m=vm):
    return p*z.diff(x)+m*z.diff(y)+(p.diff(x)+m.diff(y))*z/2+p.diff(x)*G*z/2-m.diff(y)*z*Gb/2
def insertion(z):
    return -vp*gm*z/s.sqrt(2)+vm*G*z*gbp/s.sqrt(2)
check('SM103 Cartan identity', D(insertion(B))+insertion(D(B))-kosmann(B))
check('SM98 Witt algebra on unrestricted B', kosmann(kosmann(B,wp,wm))-kosmann(kosmann(B),wp,wm)-kosmann(B,vp*wp.diff(x)-wp*vp.diff(x),vm*wm.diff(y)-wm*vm.diff(y)))
check('SM109 phi is a scalar', kosmann(B)[1,1]-vp*phi.diff(x)-vm*phi.diff(y))
for sign in (1,-1):
    check(f'SM94/98 Kosmann compatibility D_{sign:+d}', kosmann(D(B,sign))-D(kosmann(B),sign))

ep, em = s.symbols('epsilon_plus epsilon_minus')
eps, epsb = s.Matrix([0,ep]), s.Matrix([em,0])
check('SM101/112 physical fermion SUSY signs', s.Matrix([(G*F*epsb)[1]-em*phi.diff(y),(Fb*eps)[0]-ep*phi.diff(x)]))
check('negative control: removing gamma(3) changes the minus-sector sign', (F*epsb)[1] != (G*F*epsb)[1])

pp, pm = Exterior.odd('psi_plus'), Exterior.odd('psi_minus')
eP, eM = Exterior.odd('eP'), Exterior.odd('eM')
for name, vector, operator, axis, expected in (
        ('unbarred', [Exterior.odd('rho'), pm], C*gp/s.sqrt(2), 0, pm*pm.diff(0)),
        ('barred', [pp, Exterior.odd('rho_prime')], C*gbm/s.sqrt(2), 1, pp*pp.diff(1))):
    kinetic = sum((vector[a]*vector[b].diff(axis)*operator[a,b]
                   for a in range(2) for b in range(2)), Exterior())
    check(f'SM95/108 {name} fermion kinetic reduction', kinetic-expected)
L = Exterior.scalar(phi.diff(x)*phi.diff(y))+pp*pp.diff(1)+pm*pm.diff(0)
check('Grassmann kinetic terms are not identically zero', bool((pp*pp.diff(1)).terms) and bool((pm*pm.diff(0)).terms))
def variation(dphi, dpp, dpm):
    return (dphi.diff(0)*phi.diff(y)+dphi.diff(1)*phi.diff(x)
            +dpp*pp.diff(1)+pp*dpp.diff(1)+dpm*pm.diff(0)+pm*dpm.diff(0))
dphi = pp*eP+pm*eM
dpp, dpm = eP*phi.diff(x), eM*phi.diff(y)
boundary = (pp*eP*phi.diff(y)).diff(0)+(pm*eM*phi.diff(x)).diff(1)
check('SM102/113 arbitrary chiral SUSY: off-shell total divergence', variation(dphi,dpp,dpm)-boundary)
check('negative control: wrong gamma(3) sign breaks SUSY', bool((variation(dphi,dpp,-dpm)-boundary).terms))
check('SM98 Witt invariance: weight-one action', variation(
    Exterior.scalar(vp*phi.diff(x)+vm*phi.diff(y)),
    pp.diff(0)*vp+pp.diff(1)*vm+pp*(vp.diff(x)/2),
    pm.diff(0)*vp+pm.diff(1)*vm+pm*(vm.diff(y)/2))
    -(L*vp).diff(0)-(L*vm).diff(1))
a0, ap, am = field('alpha0'),field('alpha_plus'),field('alpha_minus')
dp = pm.diff(0)*a0+pp.diff(1)*ap
dm = pp.diff(1)*a0+pm.diff(0)*am
check('SM115 trivial variations: off-shell total divergence', variation(Exterior(),dp,dm)-(pp*dp).diff(1)-(pm*dm).diff(0))
check('SM115 trivial variations vanish on shell', dp.on_shell()+dm.on_shell())
e1P,e2P,e1M,e2M = [Exterior.odd(n) for n in ('e1P','e2P','e1M','e2M')]
vP,vM = 2*e1P*e2P,2*e1M*e2M
comm_phi = e1P*phi.diff(x)*e2P+e1M*phi.diff(y)*e2M-e2P*phi.diff(x)*e1P-e2M*phi.diff(y)*e1M
check('SM114 scalar SUSY closure off shell',comm_phi-vP*phi.diff(x)-vM*phi.diff(y))
for name,fermion,axis,e1,e2 in [('plus',pp,0,e1P,e2P),('minus',pm,1,e1M,e2M)]:
    comm = e2*(pp*e1P+pm*e1M).diff(axis)-e1*(pp*e2P+pm*e2M).diff(axis)
    v = (vP,vM)[axis]
    residual = comm-vP*fermion.diff(0)-vM*fermion.diff(1)-v.diff(axis)*fermion*s.Rational(1,2)
    check(f'SM114 {name} SUSY closure modulo fermion equations', residual.on_shell())

if __name__ == '__main__':
    passed = sum(c['passed'] for c in checks)
    print(json.dumps({'passed':passed,'total':len(checks),'python':sys.version.split()[0], 'sympy':s.__version__}))
    sys.exit(0 if passed == len(checks) else 1)
