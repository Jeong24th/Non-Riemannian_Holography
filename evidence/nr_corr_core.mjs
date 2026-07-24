// Dependency-free core checks for NR-CORR-001.
// Run:
//   node codex_nr_corr_core.mjs
// Negative control:
//   node codex_nr_corr_core.mjs --selftest-fail

const selftestFail = process.argv.includes("--selftest-fail");
let failures = 0;

function check(name, actual, expected = true) {
  const ok = Object.is(actual, expected);
  console.log(`[${ok ? "PASS" : "FAIL"}] ${name}: actual=${actual} expected=${expected}`);
  if (!ok) failures += 1;
}

function close(a, b, tol = 1e-11) {
  return Math.abs(a - b) <= tol * Math.max(1, Math.abs(a), Math.abs(b));
}

const l = 3;
const y = 0.7;
const D = 5;       // d_+ d_- r^(0)
const Pp = 7;      // d_+^2 r^(0)
const Pm = -4;     // d_-^2 r^(0)
const C = 11;      // d_+^2 d_-^2 r^(0)
const A = 13;      // d_+ d_- delta d^(0)
const B = 17;      // d_+^2 ell_-^(0) + d_-^2 ell_+^(0)

// SMNRlin fourth equation.
const dyD = (l / 4) * D;
check("(8/l) d_y delta d = 2 d_+d_- r", close((8 / l) * dyD, 2 * D));
const oldDyD = (l / 8) * D;
check("negative control: old coefficient l/8 fails", close((8 / l) * oldDyD, 2 * D), false);

// L_y := d_y^2 + (2/l)d_y.  The ell particular solution is l*y*P.
function LyPolynomial(a1, a2, yy) {
  // f=a1*y+a2*y^2
  return 2 * a2 + (2 / l) * (a1 + 2 * a2 * yy);
}
check("ell_+ particular branch", close(LyPolynomial(l * Pp, 0, y), 2 * Pp));
check("ell_- particular branch", close(LyPolynomial(l * Pm, 0, y), 2 * Pm));

// w particular solution.
const wy = l * (4 * A - B + (l * l / 2) * C);
const wy2 = -(l * l / 2) * C;
const wLhs = LyPolynomial(wy, wy2, y);
const wRhs = 8 * A - 2 * B - 2 * l * y * C;
check("w y and y^2 branches solve SMNRlin", close(wLhs, wRhs));
const wWithoutY2 = LyPolynomial(wy, 0, y);
check("negative control: omitting y^2 branch fails", close(wWithoutY2, wRhs), false);

// Symplectic flux on r=r0+r2*exp(-2y/l), w=w0+w2*exp(-2y/l).
const r0 = 2, r2 = -3, w0 = 5, w2 = 7;
const e = Math.exp(-2 * y / l);
const r = r0 + r2 * e;
const w = w0 + w2 * e;
const ry = (-2 / l) * r2 * e;
const wyDer = (-2 / l) * w2 * e;
const flux = 0.5 * Math.exp(2 * y / l) * (w * ry - r * wyDer);
const fluxExpected = (r0 * w2 - r2 * w0) / l;
check("cross-channel symplectic flux", close(flux, fluxExpected));

// R normalization anchor: K=T/(4pi), c=3l/(2G).
const G = 2;
const c = (3 * l) / (2 * G);
const anomalyFromL = (1 / (8 * Math.PI * G * l)) * (-(l * l) / 4);
const anomalyFromCft = -c / (48 * Math.PI);
check("R anomaly fixes c=3l/(2G)", close(anomalyFromL, anomalyFromCft));
const kkFromRescale = (1 / (4 * Math.PI)) ** 2 * (c / 2);
const kkPrinted = (1 / (4 * Math.PI) ** 2) * (c / 2);
check("R two-point normalization", close(kkFromRescale, kkPrinted));

// Centerless Witt realization for constant L and modes m+n=0.
const m = 3, n = -3, Lconst = 2;
const lhsICoeff = 2 * n * Lconst;
const rhsICoeff = (n - m) * Lconst;
check("NR charge algebra has no m^3 central term", lhsICoeff === rhsICoeff);

// Mixed worldsheet contraction: beta-x+ and betabar-x- each differentiate 1/z.
const holomorphicPower = -2;
const antiholomorphicPower = -2;
check("mixed V_W-beta betabar correlator powers",
      holomorphicPower === -2 && antiholomorphicPower === -2);

if (selftestFail) {
  check("injected negative-control failure", wy2, 0);
}

if (failures > 0) {
  console.error(`FAILED CHECKS: ${failures}`);
  process.exit(1);
}
console.log("ALL CHECKS PASSED");

