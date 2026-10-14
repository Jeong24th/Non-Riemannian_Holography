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
y, l = s.symbols('y l', positive=True)
Lp, Lm, chi, sigma, W = s.symbols('Lplus Lminus chi sigma calW')
# Directly use exp(2y/l); u_pm below retain their manuscript meaning.
radial_exponential = s.exp(2*y/l)
j = s.BlockMatrix([[s.zeros(3),s.eye(3)],[s.eye(3),s.zeros(3)]]).as_explicit()
g = s.Matrix([[2*Lp,-radial_exponential-Lp*Lm/radial_exponential,0],[-radial_exponential-Lp*Lm/radial_exponential,2*Lm,0],[0,0,1]])
b = (radial_exponential+Lp*Lm/radial_exponential)*s.Matrix([[0,-1,0],[1,0,0],[0,0,0]])
hR = s.BlockMatrix([[g.inv(),-g.inv()*b],[b*g.inv(),g-b*g.inv()*b]]).as_explicit().applyfunc(s.factor)
check('01 R: H J H = J',hR*j*hR-j)
check('01 R: measure squared = -det(g)',(radial_exponential-Lp*Lm/radial_exponential)**2+g.det())
c, h, e = s.cosh(chi),s.sinh(chi),s.exp(sigma)
hNR = s.Matrix([[0,0,0,c,-h/e,0],[0,0,0,e*h,-c,0],[0,0,1,0,0,0],
                [c,e*h,0,-W*e*h,W*c,0],[-h/e,-c,0,W*c,-W*h/e,0],[0,0,0,0,0,1]])
check('01 NR: H J H = J',hNR*j*hNR-j)
check('01 NR: inverse-metric block',hNR[:3,:3]-s.diag(0,0,1))
hinf = s.Matrix([[0,0,0,1,0,0],[0,0,0,0,-1,0],[0,0,1,0,0,0],
                 [1,0,0,0,0,0],[0,-1,0,0,0,0],[0,0,0,0,0,1]])
check('01 R boundary',hR.applyfunc(lambda z:s.limit(z,y,s.oo))-hinf)
check('01 NR boundary at W0=0',hNR.subs({chi:0,W:0})-hinf)

# 02: reduced radial equation, boundary data and finite action.
xPlus,xMinus,xPrimePlus,xPrimeMinus = s.symbols('xPlus xMinus xPrimePlus xPrimeMinus')
Lplus,Lminus = s.Function('Lplus')(xPlus),s.Function('Lminus')(xMinus)
G,mu,Y = s.symbols('G mu Y',positive=True)
W0,W1 = s.Function('W0')(xPlus,xMinus),s.Function('W1')(xPlus,xMinus)
calPi = Lplus*Lminus
u = {1:Lplus**(-s.Rational(1,2)), -1:Lminus**(-s.Rational(1,2))}
calA = {1:s.diff(u[1],xPlus,2)/u[1], -1:s.diff(u[-1],xMinus,2)/u[-1]}
calI = s.Function('calI')(chi)
calIPrime = 4*s.sqrt(2)*(s.sinh(chi)/s.sinh(chi/s.sqrt(2))-s.sqrt(2))
calW = W0+W1*chi/(2*s.sqrt(calPi))+l**2*s.sqrt(calPi)/16*sum(
    calA[sign]*u[sign]**2*calI-2*sign*(s.diff(u[1],xPlus)+sign*s.diff(u[-1],xMinus))**2
    *(s.exp(sign*chi)-1-sign*chi) for sign in (-1,1))
# Input: reduced radial EDFE; its full tensor derivation is outside this file.
radialSource = l**2*s.sqrt(calPi)/16*sum(
    calA[sign]*u[sign]**2*s.diff(calIPrime,chi)-2*sign
    *(s.diff(u[1],xPlus)+sign*s.diff(u[-1],xMinus))**2*s.exp(sign*chi) for sign in (-1,1))
check('02 radial hair equation',s.diff(calW,chi,2).subs(s.diff(calI,chi,2),s.diff(calIPrime,chi))-radialSource)
check('02 W0 boundary value',calW.subs(chi,0).subs(calI.subs(chi,0),0)-W0)
check('02 W1 radial slope',s.diff(calW,chi).subs(s.diff(calI,chi),0).subs(chi,0)-W1/(2*s.sqrt(calPi)))
check('02 integral leading order',s.limit(calIPrime/chi**2,chi,0)-s.Rational(2,3))
# SM2-SM5: constant positive L_pm, prescribed endpoint, no extra interior action.
measure = s.exp(2*y/l)-mu*s.exp(-2*y/l)
flux = -2*s.diff(measure,y)
check('02 radial boundary flux',flux+4*(s.exp(2*y/l)+mu*s.exp(-2*y/l))/l)
Sren = (-flux.subs(y,Y)-8*s.sqrt(mu)/l-4*measure.subs(y,Y)/l)/(16*s.pi*G)
check('02 finite saddle action',s.limit(Sren,Y,s.oo)+8*s.sqrt(mu)/(16*s.pi*G*l))

# 03: source operators and prescribed separated-point kernels.
# SM27 factor 1/(4i), SM43 response, SM47 prescription, SM53-SM55 NR response.
# The - sector follows by interchanging + and -.
deltaPlus,deltaMinus = xPlus-xPrimePlus,xMinus-xPrimeMinus
inverseDerivative = 1/(2*s.pi*s.I*deltaPlus)
calKplusNR = -2*Lplus*s.diff(inverseDerivative,xPlus)-s.diff(Lplus,xPlus)*inverseDerivative
calKplusR = calKplusNR+l**2*s.diff(inverseDerivative,xPlus,3)/4
AplusNR,AplusR = [z/(64*s.pi*s.I*G*l) for z in (calKplusNR,calKplusR)]
MplusNR = (-2*W1*s.diff(inverseDerivative,xPlus)-s.diff(W1,xPlus)*inverseDerivative
           +l**2*Lminus*s.diff(inverseDerivative,xPlus,3))/(256*s.pi*s.I*G*l)
check('03 NR diagonal kernel',AplusNR-(s.diff(Lplus,xPlus)/deltaPlus-2*Lplus/deltaPlus**2)/(128*s.pi**2*G*l))
check('03 R central term',AplusR-AplusNR-3*l/(256*s.pi**2*G*deltaPlus**4))
check('03 NR mixed kernel',MplusNR-(s.diff(W1,xPlus)/deltaPlus-2*W1/deltaPlus**2
                                   +6*l**2*Lminus/deltaPlus**4)/(512*s.pi**2*G*l))
hair_self = (l**6/32)*s.diff(s.I/(2*s.pi*deltaPlus*deltaMinus),xPlus,3,xMinus,3)/(64*s.pi*s.I*G*l)
check('03 R-vacuum hair self-correlator',hair_self-9*l**5/(1024*s.pi**2*G*deltaPlus**4*deltaMinus**4))

# 04: SM58 auxiliary fields, SM60-SM62 clock forms, SM82, SM85 -> SM87.
plus, minus, bp, bm, beta, betab, F, radial_term = s.symbols('plus minus barredPlus barredMinus beta betaBar F radial')
# plus/minus=partial x^+/x^-, bp/bm=bar-partial x^+/x^-; radial=partial y bar-partial y.
first = radial_term+2*Lp*plus*bp+2*Lm*minus*bm+beta*bp+betab*minus+beta*betab/(2*F)
aux = s.solve([s.diff(first,beta),s.diff(first,betab)],(beta,betab))
check('04 auxiliary worldsheet reduction',first.subs(aux)-(radial_term+2*Lp*plus*bp+2*Lm*minus*bm-2*F*minus*bp))
clock = s.Matrix([s.cosh(chi/2),-s.exp(-sigma)*s.sinh(chi/2)])
bar_clock = s.Matrix([-s.exp(sigma)*s.sinh(chi/2),s.cosh(chi/2)])
check('04 clock wedge',s.Matrix.hstack(clock,bar_clock).det()-1)
check('04 hair metric block',W*(clock*bar_clock.T+bar_clock*clock.T)-W*s.Matrix([[-e*h,c],[c,-h/e]]))
lambdaMatrices = [s.sqrt(2)*s.Matrix([[0,0],[1,0]]),-s.sqrt(2)*s.Matrix([[0,1],[0,0]]),s.diag(1,-1)]
eta = s.Matrix([[0,-1,0],[-1,0,0],[0,0,1]])
check('04 three-dimensional Clifford algebra',s.Matrix.vstack(*[
    lambdaMatrices[a]*lambdaMatrices[b]+lambdaMatrices[b]*lambdaMatrices[a]-2*eta[a,b]*s.eye(2)
    for a in range(3) for b in range(3)]))
sPlusI,sPlusJ = s.Function('sPlusI')(xPlus),s.Function('sPlusJ')(xPlus)
kPlus = sPlusI*sPlusJ
hillRules = {}
for spinor in (sPlusI,sPlusJ):
    hillRules[s.diff(spinor,xPlus,2)] = 2*Lplus*spinor/l**2
    hillRules[s.diff(spinor,xPlus,3)] = s.diff(2*Lplus*spinor/l**2,xPlus)
check('04 Hill bilinear is a local isometry',(kPlus*s.diff(Lplus,xPlus)
      +2*Lplus*s.diff(kPlus,xPlus)-l**2*s.diff(kPlus,xPlus,3)/4).subs(hillRules,simultaneous=True))

# 05: SM105-SM108. Flat oplus/ominus indices differ from coordinate +/-.
gammaOplus = s.sqrt(2)*s.Matrix([[0,1],[0,0]])
gammaOminus = s.sqrt(2)*s.Matrix([[0,0],[-1,0]])
barGammaOplus,barGammaOminus = gammaOplus,-gammaOminus
gamma3 = s.diag(1,-1)
barGamma3 = -gamma3
chargeC = barChargeC = s.Matrix([[0,1],[1,0]])
calB = s.Matrix([[s.Function('c')(xPlus,xMinus),s.Function('bPlus')(xPlus,xMinus)],
                 [s.Function('bMinus')(xPlus,xMinus),s.Function('phi')(xPlus,xMinus)]])
Dplus = lambda z:(gammaOplus*z.diff(xPlus)+gamma3*z.diff(xMinus)*barGammaOminus)/s.sqrt(2)
Dminus = lambda z:(gammaOplus*z.diff(xPlus)-gamma3*z.diff(xMinus)*barGammaOminus)/s.sqrt(2)
check('05 D_+ nilpotency and RR invariance',Dplus(Dplus(calB)))
check('05 D_- nilpotency',Dminus(Dminus(calB)))
calF = Dplus(calB)
barCalF = barChargeC.inv()*calF.T*chargeC
check('05 scalar action',-s.trace(calF*barCalF)/2-calB[1,1].diff(xPlus)*calB[1,1].diff(xMinus))
check('05 unbarred fermion kinetic term',chargeC*gammaOplus/s.sqrt(2)-s.diag(0,1))
check('05 barred fermion kinetic term',barChargeC*barGammaOminus/s.sqrt(2)-s.diag(1,0))
# SM97, SM103: weight-1/2 Kosmann derivative, no field chirality constraint.
vPlus,vMinus = s.Function('vPlus')(xPlus),s.Function('vMinus')(xMinus)
insertion = lambda z:(-vPlus*gammaOminus*z+vMinus*gamma3*z*barGammaOplus)/s.sqrt(2)
kosmann = lambda z:(vPlus*z.diff(xPlus)+vMinus*z.diff(xMinus)+(vPlus.diff(xPlus)+vMinus.diff(xMinus))*z/2
                    +vPlus.diff(xPlus)*gamma3*z/2-vMinus.diff(xMinus)*z*barGamma3/2)
check('05 Cartan anticommutator equals Kosmann',Dplus(insertion(calB))+insertion(Dplus(calB))-kosmann(calB))

if __name__ == '__main__':
    passed = sum(r['passed'] for r in results)
    print(json.dumps({'passed':passed,'total':len(results),'python':sys.version.split()[0],'sympy':s.__version__}))
    sys.exit(0 if passed==len(results) else 1)
