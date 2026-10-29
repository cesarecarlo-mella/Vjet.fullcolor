"""evaluate the 84 finite remainders of ../../Zgqq_finite_remainders at points of the decay triangle:
alpha1-3 = HA1-3 and gamma1-3 = HA4-6, each with 14 colour layers "Nc^c*Nf^f*NfZ^z", as
    sum over rows {r, g, m, q} of  q * RationalFunctionBasis[r] * T[g] * TranscendentalMultipliers[m],
with the generators T[g] from gen_eval (tridents_bs + out/generators_bs.bin).

    from remainder_eval import evaluate_remainders, KEYS
    R, info = evaluate_remainders([(0.3, 0.2), (0.1, 0.7)])    # R.shape == (2, 84), column j is KEYS[j]
    KEYS[j] == ('alpha1', 'Nc^0*Nf^0*NfZ^1')                   # LABELS[j] == 'C0F0Z1_HA1' (old file names)
    python3 remainder_eval.py x y                              # prints the 84 values

Variables: (x, y) of the alphabet (letter 1 = x, letter 2 = y, letter 3 = 1 - x - y).  The remainder files use
t = y, u = 1 - x - y.
"""
import os, re, sys
import numpy as np
import gen_eval

HERE = os.path.dirname(os.path.abspath(__file__))
ZG = os.path.join(HERE, '..', 'amplitudes')                     # in the shipped package
if not os.path.isdir(ZG):
    ZG = os.path.join(HERE, '..', '..', 'Zgqq_finite_remainders')   # in the repository
NAMES = ['alpha1', 'alpha2', 'alpha3', 'gamma1', 'gamma2', 'gamma3']
HA = {n: i + 1 for i, n in enumerate(NAMES)}


def _items(fn):
    """the elements of a Mathematica list {a, b, ...} in the file (comments removed)"""
    src = re.sub(r'\(\*.*?\*\)', '', open(os.path.join(ZG, fn)).read(), flags=re.S).strip()
    return [s.strip() for s in src[1:-1].split(',')]


# rational functions of (t, u); the multipliers are constants times a power of r = +Sqrt[x y z] (= Sqrt[t u (1-t-u)])
_RFB = [eval('lambda t, u: ' + s.replace('^', '**')) for s in _items('RationalFunctionBasis.wl')]
_CST = dict(Pi=np.pi, I=1j, z3=1.2020569031595942854, z5=1.0369277551433699263)


def _tm(s, r):
    return complex(eval(s.replace('Zeta[3]', 'z3').replace('Zeta[5]', 'z5').replace('^', '**'), dict(_CST, r=r)))


TM_STR = _items('TranscendentalMultipliers.wl')
TM = np.array([_tm(s, 1.0) for s in TM_STR])                          # the constant
RPOW = [int(round(np.log2(abs(_tm(s, 2.0) / _tm(s, 1.0))))) for s in TM_STR]    # the power of r
assert all(abs(_tm(s, 3.0) - c * 3.0 ** k) < 1e-12 * abs(c) for s, c, k in zip(TM_STR, TM, RPOW)), 'multiplier not c * r^k'

# rows, merged per (r, g, power of r): coefficient sum_m q * TM[m]
KEYS, LABELS, ROWS = [], [], []
for n in NAMES:
    src = open(os.path.join(ZG, 'FiniteRemainders', n + '.wl')).read()
    body = src[src.index('*)') + 2:]
    for key, rows in re.findall(r'"([^"]+)"\s*->\s*\{(.*?)\}\s*(?=,\s*"|\s*\}\s*$)', body, re.S):
        c, f, z = re.fullmatch(r'Nc\^(-?\d+)\*Nf\^(\d+)\*NfZ\^(\d+)', key).groups()
        acc = {}
        for r, g, m, q in re.findall(r'\{(\d+),(\d+),(\d+),(-?\d+(?:/\d+)?)\}', rows):
            a, b = q.split('/') if '/' in q else (q, '1')
            k = (int(r) - 1, int(g) - 1, RPOW[int(m) - 1])
            acc[k] = acc.get(k, 0) + int(a) / int(b) * TM[int(m) - 1]
        KEYS.append((n, key))
        LABELS.append('C%dF%sZ%s_HA%d' % (int(c) % 10, f, z, HA[n]))
        ROWS.append(tuple(np.array([k[i] for k in acc]) for i in range(3)) + (np.array(list(acc.values())),))


# the error estimate: ERR_FACTOR times the difference of the two evaluations (different paths: the path-dependent
# integration error), plus ERR_EPS times sum |terms| of the sum over rows (the systematic error of the generators, the
# same on every path -- the double system's matrices are rounded -- amplified by the cancellation in the sum, which
# reaches 1e9 near the edges).  Against full-quad values (2026-10-29, 8 points x 84): the actual error divided by
# sum |terms| was at most 5e-14 (at the corner point (1/30, 1/30)); ERR_EPS is provisional: it overestimates ~500x
# in the median.  To be calibrated / replaced after the 1000-point study (accuracy_study/).
ERR_FACTOR = 10
ERR_EPS = 1e-13


def remainders_from(points, T, rsign=1, absolute=False):
    """the 84 remainders from generator values T (npts, 579) at the (physical) points; r = rsign Sqrt[x y z]
    (rsign = -1 in the production regions: the continued root).  absolute=True: sum |terms| instead (real)"""
    P = np.asarray(points, float)
    t, u = P[:, 1], 1 - P[:, 0] - P[:, 1]
    RF = np.array([f(t, u) * np.ones_like(t) for f in _RFB]).T          # (npts, 129)
    rr = rsign * np.sqrt(P[:, 0] * P[:, 1] * (1 - P[:, 0] - P[:, 1]))
    out = np.zeros((len(P), len(KEYS)), float if absolute else complex)
    for j, (r, g, k, c) in enumerate(ROWS):
        if absolute:
            out[:, j] = np.abs(RF[:, r] * T[:, g] * rr[:, None] ** k) @ np.abs(c)
        else:
            out[:, j] = (RF[:, r] * T[:, g] * rr[:, None] ** k) @ c
    return out


def evaluate_remainders(points, tol=1e-10, error=False, region='decay'):
    """the 84 remainders at the points: (x, y) of the decay triangle, or (X, Y) of the triangle of a production region
    ('x', 'y', 'z': A > 1, see REGION_MAP; out/region_<A>_bs.bin, rational functions at the physical point, continued
    root r = -Sqrt[x y z]).  A point closer than NEAR to a spurious curve is replaced by the symmetric average across
    the curve (gen_eval.plan, with the letters of the region): this covers both the poles of the DE and those of the
    rational functions, which cancel in the remainder.
    error=True: every point is evaluated again on a different complexified path (and, for an averaged point, with half
    the stencil step); the error estimate ERR_FACTOR |difference| + ERR_EPS sum |terms| is returned: (R, info, err).
    Near the edges of the triangle the rational functions (up to 1/x^3, 1/y^3, 1/z^3, and (x+y)^-5 etc. at the corners)
    cancel against large logarithms, and double precision loses digits there; the estimate shows how many."""
    if region == 'decay':
        data, data2, rsign, phys = DATA_DECAY, gen_eval.second_data(), 1, (lambda P: P)
    else:
        data = data2 = os.path.join(HERE, 'out', 'region_%s_bs.bin' % region)
        rsign, phys = -1, (lambda P: region_to_phys(region, P))
    flat, pl = gen_eval.plan([tuple(p) for p in points], region=region)
    T, infos = gen_eval.evaluate_direct(flat, tol, data=data, region=region)
    R, info = gen_eval.collapse(remainders_from(phys(flat), T, rsign), infos, pl)
    if not error:
        return R, info
    flat2, pl2 = gen_eval.plan([tuple(p) for p in points], hfac=0.5, region=region)
    T2, infos2 = gen_eval.evaluate_direct(flat2, tol, data=data2, second=True, region=region)
    R2, _ = gen_eval.collapse(remainders_from(phys(flat2), T2, rsign), infos2, pl2)
    S, _ = gen_eval.collapse(remainders_from(phys(flat), T, rsign, absolute=True), infos, pl)
    return R, info, ERR_FACTOR * np.abs(R - R2) + ERR_EPS * np.abs(S)


DATA_DECAY = gen_eval.DATA
REGION_MAP = {'z': lambda X, Y, Z: (-X / Z, -Y / Z), 'x': lambda X, Y, Z: (1 / Z, -X / Z), 'y': lambda X, Y, Z: (-X / Z, 1 / Z)}


def region_to_phys(region, points):
    """points (X, Y) of the triangle of the production region A > 1 (A = 1/Z) -> physical (x, y)"""
    return [REGION_MAP[region](X, Y, 1 - X - Y) for X, Y in points]


def evaluate_remainders_region(region, points, tol=1e-12, error=False):
    """the 84 remainders in the production region A > 1 at points (X, Y) of its triangle: evaluate_remainders"""
    return evaluate_remainders(points, tol, error, region=region)


if __name__ == '__main__':
    R, info = evaluate_remainders([(float(sys.argv[1]), float(sys.argv[2]))])
    for (n, key), lab, v in zip(KEYS, LABELS, R[0]):
        print('%-7s %-17s (%-11s) = %.15e %+.15e i' % (n, key, lab, v.real, v.imag))
    print(info[0])
