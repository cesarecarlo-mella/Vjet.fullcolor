#!/usr/bin/env python3
"""zgqq_eval -- full-colour three-loop finite remainders of the helicity amplitudes for V -> qbar q g, at any point
of the decay region.  See README.md.

Python:
    import zgqq_eval as Z
    res = Z.evaluate([(0.2, 0.3), (0.1, 0.5)])
    res.values[p, i, l]    complex value: point p, coefficient Z.COEFFICIENTS[i], colour layer Z.LAYERS[l]
    res.errors[p, i, l]    absolute error estimate (same shape; None with error=False)
    v, e = Z.colour_sum(res, Nc=3, Nf=5, NfZ=...)        -> the six coefficients, shape (npts, 6), and their errors

Command line:
    python3 zgqq_eval.py X Y                              all 84 values (and errors) at one point
    python3 zgqq_eval.py --points FILE [--out OUT.csv | OUT.json]
    python3 zgqq_eval.py X Y --Nc 3 --Nf 5 --NfZ 1        the six coefficients summed over the colour layers
"""
import argparse, collections, json, os, re, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'core'))
import gen_eval          # noqa: E402
import remainder_eval    # noqa: E402

__version__ = open(os.path.join(HERE, 'VERSION')).read().strip() if os.path.isfile(os.path.join(HERE, 'VERSION')) else 'dev'

COEFFICIENTS = list(remainder_eval.NAMES)                         # alpha1, alpha2, alpha3, gamma1, gamma2, gamma3
LAYERS = [k for n, k in remainder_eval.KEYS if n == COEFFICIENTS[0]]    # 'Nc^c*Nf^f*NfZ^z', the same for all six
POWERS = [tuple(int(v) for v in re.fullmatch(r'Nc\^(-?\d+)\*Nf\^(\d+)\*NfZ\^(\d+)', k).groups()) for k in LAYERS]
_COL = np.array([[remainder_eval.KEYS.index((n, k)) for k in LAYERS] for n in COEFFICIENTS])    # (6, 14) -> column

Result = collections.namedtuple('Result', 'points values errors info')


def _check(points):
    P = np.atleast_2d(np.asarray(points, float))
    if P.ndim != 2 or P.shape[1] != 2:
        raise ValueError('points must be pairs (x, y)')
    bad = ~((P[:, 0] > 0) & (P[:, 1] > 0) & (P[:, 0] + P[:, 1] < 1))
    if bad.any():
        raise ValueError('points outside the decay region x > 0, y > 0, x + y < 1: %s' % P[bad][:5].tolist())
    return [tuple(p) for p in P]


def evaluate(points, tol=1e-12, error=True):
    """the 84 finite remainders at the points (x, y) = (s12, s13) / q^2, with z = s23 / q^2 = 1 - x - y.
    Returns Result(points, values, errors, info): values[p, i, l] for COEFFICIENTS[i] and LAYERS[l]; errors the same
    shape (a conservative absolute error estimate, from a second, independent evaluation), or None if error=False."""
    pts = _check(points)
    if error:
        R, info, err = remainder_eval.evaluate_remainders(pts, tol, error=True)
        return Result(pts, R[:, _COL], err[:, _COL], info)
    R, info = remainder_eval.evaluate_remainders(pts, tol)
    return Result(pts, R[:, _COL], None, info)


def colour_sum(res, Nc=3.0, Nf=5.0, NfZ=None):
    """sum over the colour layers, sum_l Nc^c Nf^f NfZ^z values[:, :, l]: (values (npts, 6), errors (npts, 6) or None).
    NfZ must be given whenever a layer contains it (see README)."""
    if NfZ is None and any(z for _, _, z in POWERS):
        raise ValueError('give NfZ (the layers with NfZ^1 contribute)')
    w = np.array([Nc ** c * Nf ** f * (NfZ ** z if z else 1.0) for c, f, z in POWERS])
    v = res.values @ w
    e = None if res.errors is None else res.errors @ np.abs(w)
    return v, e


def _read_points(fn):
    pts = []
    for line in open(fn):
        line = line.split('#')[0].replace(',', ' ').split()
        if line:
            pts.append((float(line[0]), float(line[1])))
    return pts


def main(argv=None):
    ap = argparse.ArgumentParser(description='full-colour three-loop finite remainders for V -> qbar q g (decay region)')
    ap.add_argument('xy', nargs='*', type=float, help='x y (= s12/q^2, s13/q^2) of one point')
    ap.add_argument('--points', help='file with one point "x y" per line')
    ap.add_argument('--out', help='write the results to OUT.csv or OUT.json instead of printing them')
    ap.add_argument('--tol', type=float, default=1e-12, help='integration tolerance (default 1e-12)')
    ap.add_argument('--no-error', action='store_true', help='skip the error estimate (twice as fast)')
    ap.add_argument('--Nc', type=float, help='also give the six coefficients summed over the colour layers')
    ap.add_argument('--Nf', type=float, default=5.0)
    ap.add_argument('--NfZ', type=float)
    ap.add_argument('--version', action='version', version='zgqq_eval ' + __version__)
    a = ap.parse_args(argv)
    if a.points:
        pts = _read_points(a.points)
    elif len(a.xy) == 2:
        pts = [tuple(a.xy)]
    else:
        ap.error('give x y, or --points FILE')
    res = evaluate(pts, a.tol, error=not a.no_error)
    err = res.errors if res.errors is not None else np.full(res.values.shape, np.nan)
    rows = [(p[0], p[1], c, l, res.values[ip, i, j], err[ip, i, j]) for ip, p in enumerate(res.points)
            for i, c in enumerate(COEFFICIENTS) for j, l in enumerate(LAYERS)]
    summed = None
    if a.Nc is not None:
        summed = colour_sum(res, a.Nc, a.Nf, a.NfZ)
    if a.out and a.out.endswith('.json'):
        out = []
        for ip, p in enumerate(res.points):
            d = dict(x=p[0], y=p[1], values={c: {l: [res.values[ip, i, j].real, res.values[ip, i, j].imag, float(err[ip, i, j])]
                                                 for j, l in enumerate(LAYERS)} for i, c in enumerate(COEFFICIENTS)})
            if summed is not None:
                d['colour_sum'] = dict(Nc=a.Nc, Nf=a.Nf, NfZ=a.NfZ, values={c: [summed[0][ip, i].real, summed[0][ip, i].imag,
                                       float(summed[1][ip, i]) if summed[1] is not None else None] for i, c in enumerate(COEFFICIENTS)})
            out.append(d)
        json.dump(dict(version=__version__, layers='[Re, Im, error]', points=out), open(a.out, 'w'), indent=1)
    elif a.out:
        with open(a.out, 'w') as f:
            f.write('x,y,coefficient,layer,re,im,error\n')
            for x, y, c, l, v, e in rows:
                f.write('%.17g,%.17g,%s,%s,%.16e,%.16e,%.3e\n' % (x, y, c, l, v.real, v.imag, e))
            if summed is not None:
                for ip, p in enumerate(res.points):
                    for i, c in enumerate(COEFFICIENTS):
                        e = summed[1][ip, i] if summed[1] is not None else np.nan
                        f.write('%.17g,%.17g,%s,sum(Nc=%g;Nf=%g;NfZ=%s),%.16e,%.16e,%.3e\n'
                                % (p[0], p[1], c, a.Nc, a.Nf, a.NfZ, summed[0][ip, i].real, summed[0][ip, i].imag, e))
    else:
        for ip, p in enumerate(res.points):
            print('# x = %.17g, y = %.17g, z = %.17g' % (p[0], p[1], 1 - p[0] - p[1]))
            print('# %-7s %-18s %24s %24s %10s' % ('coeff', 'colour layer', 'Re', 'Im', 'error'))
            for i, c in enumerate(COEFFICIENTS):
                for j, l in enumerate(LAYERS):
                    v = res.values[ip, i, j]
                    print('  %-7s %-18s %24.16e %24.16e %10.2e' % (c, l, v.real, v.imag, err[ip, i, j]))
            if summed is not None:
                print('# summed over the colour layers with Nc = %g, Nf = %g, NfZ = %s' % (a.Nc, a.Nf, a.NfZ))
                for i, c in enumerate(COEFFICIENTS):
                    v = summed[0][ip, i]; e = summed[1][ip, i] if summed[1] is not None else np.nan
                    print('  %-7s %-18s %24.16e %24.16e %10.2e' % (c, 'sum', v.real, v.imag, e))


if __name__ == '__main__':
    main()
