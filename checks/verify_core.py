"""Independent SymPy counterparts of the five short Mathematica files.
Run: python checks/verify_core.py. No manuscript or other file input.
These are scoped algebraic checks, not a full tensor-EDFE derivation.
"""
import json
import sys
import sympy as s

results = []
def check(name, expression):
    terms = list(expression) if isinstance(expression, s.MatrixBase) else [expression]
    passed = all(s.simplify(s.expand(t)) == 0 for t in terms)
    results.append({'check':name,'passed':passed})
    print(('[PASS] ' if passed else '[FAIL] ')+name)
    if not passed:
        print(expression)

# 01: generalized metrics and shared leading boundary geometry.
u, Lp, Lm, chi, sigma, W = s.symbols('u Lp Lm chi sigma W')
j = s.BlockMatrix([[s.zeros(3),s.eye(3)],[s.eye(3),s.zeros(3)]]).as_explicit()
g = s.Matrix([[2*Lp,-u-Lp*Lm/u,0],[-u-Lp*Lm/u,2*Lm,0],[0,0,1]])
b = (u+Lp*Lm/u)*s.Matrix([[0,-1,0],[1,0,0],[0,0,0]])
hR = s.BlockMatrix([[g.inv(),-g.inv()*b],[b*g.inv(),g-b*g.inv()*b]]).as_explicit().applyfunc(s.factor)
check('01 R: H J H = J',hR*j*hR-j)
check('01 R: measure squared = -det(g)',(u-Lp*Lm/u)**2+g.det())
c, h, e = s.cosh(chi),s.sinh(chi),s.exp(sigma)
hNR = s.Matrix([[0,0,0,c,-h/e,0],[0,0,0,e*h,-c,0],[0,0,1,0,0,0],
                [c,e*h,0,-W*e*h,W*c,0],[-h/e,-c,0,W*c,-W*h/e,0],[0,0,0,0,0,1]])
check('01 NR: H J H = J',hNR*j*hNR-j)
check('01 NR: inverse-metric block',hNR[:3,:3]-s.diag(0,0,1))
hinf = s.Matrix([[0,0,0,1,0,0],[0,0,0,0,-1,0],[0,0,1,0,0,0],
                 [1,0,0,0,0,0],[0,-1,0,0,0,0],[0,0,0,0,0,1]])
check('01 R boundary',hR.applyfunc(lambda z:s.limit(z,u,s.oo))-hinf)
check('01 NR boundary at W0=0',hNR.subs({chi:0,W:0})-hinf)

# 02: reduced radial equation, boundary data and finite action.
l, GN, product, A, P, M, W0, W1, mu = s.symbols('l GN product A P M W0 W1 mu',nonzero=True)
radial = s.Function('radialIntegral')(chi)
radial_prime = 4*s.sqrt(2)*(s.sinh(chi)/s.sinh(chi/s.sqrt(2))-s.sqrt(2))
hair = W0+W1*chi/(2*s.sqrt(product))+l**2*s.sqrt(product)/16*(
    A*radial-2*P**2*(s.exp(chi)-1-chi)+2*M**2*(s.exp(-chi)-1+chi))
source = l**2*s.sqrt(product)/16*(A*s.diff(radial_prime,chi)-2*P**2*s.exp(chi)+2*M**2*s.exp(-chi))
check('02 radial hair equation',s.diff(hair,chi,2).subs(s.diff(radial,chi,2),s.diff(radial_prime,chi))-source)
check('02 W0 boundary value',hair.subs(chi,0).subs(radial.subs(chi,0),0)-W0)
check('02 W1 radial slope',s.diff(hair,chi).subs(s.diff(radial,chi),0).subs(chi,0)-W1/(2*s.sqrt(product)))
check('02 integral leading order',s.limit(radial_prime/chi**2,chi,0)-s.Rational(2,3))
measure = u-mu/u
flux = -2*(2*u/l)*s.diff(measure,u)
check('02 radial boundary flux',flux+4*(u+mu/u)/l)
renormalized = (-flux-8*s.sqrt(mu)/l-4*measure/l)/(16*s.pi*GN)
check('02 finite saddle action',s.limit(renormalized,u,s.oo)+8*s.sqrt(mu)/(16*s.pi*GN*l))

# 03: source operators and prescribed separated-point kernels.
d, dp, dm, LL, Lprime, WW, Wprime, Lopposite = s.symbols('delta deltaPlus deltaMinus LL Lprime WW Wprime Lopposite')
inverse = 1/(2*s.pi*s.I*d)
normalization = 1/(16*s.pi*GN*l)
ward_nr = 2*LL*s.diff(inverse,d)+Lprime*inverse
ward_r = ward_nr-l**2*s.diff(inverse,d,3)/4
hair_response = 2*WW*s.diff(inverse,d)+Wprime*inverse-l**2*Lopposite*s.diff(inverse,d,3)
diagonal_nr, diagonal_r = [-normalization*z/(4*s.I) for z in (ward_nr,ward_r)]
mixed_nr = -normalization*hair_response/(16*s.I)
check('03 NR diagonal kernel',diagonal_nr-(Lprime/d-2*LL/d**2)/(128*s.pi**2*GN*l))
check('03 R central term',diagonal_r-diagonal_nr-3*l/(256*s.pi**2*GN*d**4))
check('03 NR mixed kernel',mixed_nr-(Wprime/d-2*WW/d**2+6*l**2*Lopposite/d**4)/(512*s.pi**2*GN*l))
hair_self = (l**6/32)*s.diff(s.I/(2*s.pi*dp*dm),dp,3,dm,3)/(64*s.pi*s.I*GN*l)
check('03 R-vacuum hair self-correlator',hair_self-9*l**5/(1024*s.pi**2*GN*dp**4*dm**4))

# 04: auxiliary fields, clocks and specified reduced spinor channels.
plus, minus, bp, bm, beta, betab, F, radial_term = s.symbols('plus minus barredPlus barredMinus beta betaBar F radial')
first = radial_term+2*Lp*plus*bp+2*Lm*minus*bm-beta*bm-betab*plus+beta*betab/(2*F)
aux = s.solve([s.diff(first,beta),s.diff(first,betab)],(beta,betab))
check('04 auxiliary worldsheet reduction',first.subs(aux)-(radial_term+2*Lp*plus*bp+2*Lm*minus*bm-2*F*plus*bm))
clock = s.Matrix([s.cosh(chi/2),-s.exp(-sigma)*s.sinh(chi/2)])
bar_clock = s.Matrix([-s.exp(sigma)*s.sinh(chi/2),s.cosh(chi/2)])
check('04 clock wedge',s.Matrix.hstack(clock,bar_clock).det()-1)
check('04 hair metric block',W*(clock*bar_clock.T+bar_clock*clock.T)-W*s.Matrix([[-e*h,c],[c,-h/e]]))
x,y = s.symbols('x y')
f,gg = s.Function('f')(x),s.Function('g')(x)
chirality = s.diag(1,-1)
dirac = lambda z:chirality*z/(s.sqrt(2)*l)+s.Matrix([[0,0],[1,0]])*z.diff(x)
positive,negative = s.Matrix([s.sqrt(2)*f,l*s.diff(f,x)]),s.Matrix([0,gg])
check('04 arbitrary chiral f channel',dirac(positive)-positive/(s.sqrt(2)*l))
check('04 arbitrary chiral g channel',dirac(negative)+negative/(s.sqrt(2)*l))

# 05: unrestricted RR bispinor, kinetic reduction and Cartan identity.
gamma_plus = s.sqrt(2)*s.Matrix([[0,1],[0,0]])
gamma_minus = s.sqrt(2)*s.Matrix([[0,0],[-1,0]])
bar_gamma_minus = -gamma_minus
C = s.Matrix([[0,1],[1,0]])
B = s.Matrix([[s.Function('c')(x,y),s.Function('bplus')(x,y)],
              [s.Function('bminus')(x,y),s.Function('phi')(x,y)]])
rr = lambda z,sign=1:(gamma_plus*z.diff(x)+sign*chirality*z.diff(y)*bar_gamma_minus)/s.sqrt(2)
check('05 D_plus nilpotency and RR invariance',rr(rr(B)))
check('05 D_minus nilpotency',rr(rr(B,-1),-1))
strength = rr(B)
bar_strength = C*strength.T*C
check('05 scalar action',-s.trace(strength*bar_strength)/2-B[1,1].diff(x)*B[1,1].diff(y))
check('05 unbarred fermion kinetic term',C*gamma_plus/s.sqrt(2)-s.diag(0,1))
check('05 barred fermion kinetic term',C*bar_gamma_minus/s.sqrt(2)-s.diag(1,0))
vp,vm = s.Function('vp')(x),s.Function('vm')(y)
insertion = lambda z:(-vp*gamma_minus*z+vm*chirality*z*gamma_plus)/s.sqrt(2)
kosmann = lambda z:vp*z.diff(x)+vm*z.diff(y)+(vp.diff(x)+vm.diff(y))*z/2+vp.diff(x)*chirality*z/2+vm.diff(y)*z*chirality/2
check('05 Cartan anticommutator equals Kosmann',rr(insertion(B))+insertion(rr(B))-kosmann(B))

if __name__ == '__main__':
    passed = sum(r['passed'] for r in results)
    print(json.dumps({'passed':passed,'total':len(results),'python':sys.version.split()[0],'sympy':s.__version__}))
    sys.exit(0 if passed==len(results) else 1)
