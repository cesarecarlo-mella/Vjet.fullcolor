#!/usr/bin/env python3
"""zgqq_eval -- full-colour three-loop finite remainders of the helicity amplitudes for V -> qbar q g, in the decay
region and, analytically continued, in the three production regions.  See README.md.

Kinematics: x = s12/q^2, y = s13/q^2, z = s23/q^2 = 1 - x - y.  Decay: 0 < x, y, z < 1.  Production region A (A = x, y
or z): A > 1 and the other two negative, parametrised by its own triangle (X, Y), X, Y > 0, X + Y < 1, Z = 1 - X - Y:
    region 'z':  (x, y, z) = (-X/Z, -Y/Z, 1/Z)        region 'x':  (x, y, z) = (1/Z, -X/Z, -Y/Z)
    region 'y':  (x, y, z) = (-X/Z, 1/Z, -Y/Z)        (Z = 1/A: the inverse of the large invariant)

Python:
    import zgqq_eval as Z
    res = Z.evaluate([(0.2, 0.3), (0.1, 0.5)])             decay, points (x, y)
    res = Z.evaluate([(0.2, 0.3)], region='z')             production region z > 1, points (X, Y) of its triangle
    res = Z.evaluate_physical([(-0.5, -0.25), (0.2, 0.3)]) physical (x, y) anywhere: the region is chosen per point
    res.values[p, i, l]    complex value: point p, coefficient Z.COEFFICIENTS[i], colour layer Z.LAYERS[l]
    res.errors[p, i, l]    absolute error estimate (same shape; None with error=False)
    res = Z.evaluate(points, precision='quad')            about 32 digits (27 guaranteed by the start values), ~5 s
    res.values + res.values_lo                            per point and path: values rounded to double, values_lo the rest
    v, e = Z.colour_sum(res, Nc=3, Nf=5, NfZ=...)        -> the six coefficients, shape (npts, 6), and their errors

Command line:
    python3 zgqq_eval.py X Y                              all 84 values (and errors) at one point
    python3 zgqq_eval.py --points FILE [--out OUT.csv | OUT.json]
    python3 zgqq_eval.py X Y --Nc 3 --Nf 5 --NfZ 1        the six coefficients summed over the colour layers
    python3 zgqq_eval.py X Y --quad                       full double-double precision (needs core/cpp/iter_quad: make)
    python3 zgqq_eval.py X Y --region z                   production region z > 1, (X, Y) of its triangle
    python3 zgqq_eval.py x y --physical                   physical (x, y), decay or production
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

Result = collections.namedtuple('Result', 'points values errors info values_lo region', defaults=(None, 'decay'))
REGIONS = ('decay', 'x', 'y', 'z')


def _check(points, region='decay'):
    if region not in REGIONS:
        raise ValueError('region must be one of %s' % (REGIONS,))
    P = np.atleast_2d(np.asarray(points, float))
    if P.ndim != 2 or P.shape[1] != 2:
        raise ValueError('points must be pairs (x, y)')
    bad = ~((P[:, 0] > 0) & (P[:, 1] > 0) & (P[:, 0] + P[:, 1] < 1))
    if bad.any():
        raise ValueError('points outside the triangle X > 0, Y > 0, X + Y < 1 (%s): %s' % (region, P[bad][:5].tolist()))
    return [tuple(p) for p in P]


def physical_point(region, X, Y):
    """(X, Y) of the triangle of a region -> physical (x, y, z)"""
    if region == 'decay':
        return X, Y, 1 - X - Y
    x, y = remainder_eval.REGION_MAP[region](X, Y, 1 - X - Y)
    return x, y, 1 - x - y


def region_point(x, y):
    """physical (x, y) -> (region, X, Y): decay if 0 < x, y, z < 1; region A if A > 1 and the other two are negative"""
    z = 1 - x - y
    if x > 0 and y > 0 and z > 0:
        return 'decay', x, y
    for A, a, b, c in (('z', z, x, y), ('x', x, y, z), ('y', y, x, z)):
        if a > 1 and b < 0 and c < 0:
            return A, -b / a, -c / a                    # Z = 1/a, X = -b Z, Y = -c Z (region map above)
    raise ValueError('(x, y, z) = (%g, %g, %g) is in no physical region (decay, or one invariant > 1 and two < 0)' % (x, y, z))


def evaluate(points, tol=None, error=True, precision='double', region='decay'):
    """the 84 finite remainders at the points: (x, y) = (s12, s13) / q^2 of the decay (z = 1 - x - y), or (X, Y) of the
    triangle of a production region ('x', 'y', 'z'; see the module docstring and physical_point).
    Returns Result(points, values, errors, info, values_lo, region): values[p, i, l] for COEFFICIENTS[i] and
    LAYERS[l]; errors the same shape (a conservative absolute error estimate, from a second evaluation on another path),
    or None if error=False.
    precision='quad': everything in double-double (about 32 digits; the start values hold 27): the exact integer
    system integrated by core/cpp/iter_quad, then the generators and the remainders.  values is the double-double value
    rounded to double, values_lo the remainder (the full value is values + values_lo).  Default tol: 1e-12 (double),
    1e-28 (quad)."""
    pts = _check(points, region)
    if precision == 'quad':
        import quad_eval
        out = quad_eval.evaluate_remainders_quad(pts, tol or 1e-28, error=error, region=region)
        R = out[0]
        hi = (R.re.hi + 1j * R.im.hi)[:, _COL]
        lo = (R.re.lo + 1j * R.im.lo)[:, _COL]
        return Result(pts, hi, out[2][:, _COL] if error else None, out[1], lo, region)
    if precision != 'double':
        raise ValueError("precision must be 'double' or 'quad'")
    tol = tol or 1e-12
    if error:
        R, info, err = remainder_eval.evaluate_remainders(pts, tol, error=True, region=region)
        return Result(pts, R[:, _COL], err[:, _COL], info, None, region)
    R, info = remainder_eval.evaluate_remainders(pts, tol, region=region)
    return Result(pts, R[:, _COL], None, info, None, region)


def evaluate_physical(points, tol=None, error=True, precision='double'):
    """the 84 remainders at physical points (x, y) in any region (region_point): Result with points the physical
    points, region a list (one per point), info[p]['region_point'] = (region, X, Y)"""
    P = [tuple(map(float, p)) for p in np.atleast_2d(np.asarray(points, float))]
    where = [region_point(*p) for p in P]
    shape = (len(P), len(COEFFICIENTS), len(LAYERS))
    vals, lo = np.zeros(shape, complex), None
    errs = np.zeros(shape) if error else None
    info = [None] * len(P)
    for reg in REGIONS:
        idx = [i for i, w in enumerate(where) if w[0] == reg]
        if not idx:
            continue
        r = evaluate([where[i][1:] for i in idx], tol, error, precision, reg)
        vals[idx] = r.values
        if error:
            errs[idx] = r.errors
        if r.values_lo is not None:
            lo = np.zeros(shape, complex) if lo is None else lo
            lo[idx] = r.values_lo
        for i, inf in zip(idx, r.info):
            info[i] = dict(inf, region_point=where[i])
    return Result(P, vals, errs, info, lo, [w[0] for w in where])


def colour_sum(res, Nc=3.0, Nf=5.0, NfZ=None):
    """sum over the colour layers, sum_l Nc^c Nf^f NfZ^z values[:, :, l]: (values (npts, 6), errors (npts, 6) or None).
    NfZ must be given whenever a layer contains it (see README)."""
    if NfZ is None and any(z for _, _, z in POWERS):
        raise ValueError('give NfZ (the layers with NfZ^1 contribute)')
    w = np.array([Nc ** c * Nf ** f * (NfZ ** z if z else 1.0) for c, f, z in POWERS])
    v = res.values @ w
    if res.values_lo is not None:                     # quad: the sum in double-double, rounded
        from ddnp import DD, dd_sum
        re_ = dd_sum(DD(res.values.real, res.values_lo.real) * w, axis=-1)
        im_ = dd_sum(DD(res.values.imag, res.values_lo.imag) * w, axis=-1)
        v = re_.to_float() + 1j * im_.to_float()
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
    ap.add_argument('--quad', action='store_true', help='full double-double precision (about 32 digits, ~5 s per point and path)')
    ap.add_argument('--region', choices=REGIONS, default='decay', help='the points are (X, Y) of this region (default decay)')
    ap.add_argument('--physical', action='store_true', help='the points are physical (x, y) in any region')
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
    prec = 'quad' if a.quad else 'double'
    if a.physical:
        res = evaluate_physical(pts, a.tol, error=not a.no_error, precision=prec)
    else:
        res = evaluate(pts, a.tol, error=not a.no_error, precision=prec, region=a.region)
    regs = res.region if isinstance(res.region, list) else [res.region] * len(res.points)
    err = res.errors if res.errors is not None else np.full(res.values.shape, np.nan)
    if res.values_lo is not None:                     # print all digits: the values as decimal strings
        import mpmath as mp
        full = lambda h, l: mp.nstr(mp.mpf(h) + l, 32, min_fixed=-mp.inf, max_fixed=-mp.inf) if np.isfinite(h) else repr(h)
        vstr = lambda ip, i, j: (full(res.values[ip, i, j].real, res.values_lo[ip, i, j].real),
                                 full(res.values[ip, i, j].imag, res.values_lo[ip, i, j].imag))
    else:
        vstr = lambda ip, i, j: ('%.16e' % res.values[ip, i, j].real, '%.16e' % res.values[ip, i, j].imag)
    summed = None
    if a.Nc is not None:
        summed = colour_sum(res, a.Nc, a.Nf, a.NfZ)
    phys = [p + (1 - p[0] - p[1],) if (a.physical or rg == 'decay') else physical_point(rg, *p) for p, rg in zip(res.points, regs)]
    if a.out and a.out.endswith('.json'):
        out = []
        for ip, p in enumerate(res.points):
            d = dict(region=regs[ip], x=phys[ip][0], y=phys[ip][1], z=phys[ip][2], values={c: {l: [res.values[ip, i, j].real, res.values[ip, i, j].imag, float(err[ip, i, j])]
                                                 for j, l in enumerate(LAYERS)} for i, c in enumerate(COEFFICIENTS)})
            if regs[ip] != 'decay':
                d['X'], d['Y'] = (res.info[ip]['region_point'][1:] if a.physical else p)
            if res.values_lo is not None:
                d['values_quad'] = {c: {l: list(vstr(ip, i, j)) for j, l in enumerate(LAYERS)} for i, c in enumerate(COEFFICIENTS)}
            if summed is not None:
                d['colour_sum'] = dict(Nc=a.Nc, Nf=a.Nf, NfZ=a.NfZ, values={c: [summed[0][ip, i].real, summed[0][ip, i].imag,
                                       float(summed[1][ip, i]) if summed[1] is not None else None] for i, c in enumerate(COEFFICIENTS)})
            out.append(d)
        json.dump(dict(version=__version__, layers='[Re, Im, error]', points=out), open(a.out, 'w'), indent=1)
    elif a.out:
        with open(a.out, 'w') as f:
            f.write('region,x,y,coefficient,layer,re,im,error\n')               # x, y: the physical point
            for ip, p in enumerate(phys):
                for i, c in enumerate(COEFFICIENTS):
                    for j, l in enumerate(LAYERS):
                        f.write('%s,%.17g,%.17g,%s,%s,%s,%s,%.3e\n' % ((regs[ip], p[0], p[1], c, l) + vstr(ip, i, j) + (err[ip, i, j],)))
            if summed is not None:
                for ip, p in enumerate(phys):
                    for i, c in enumerate(COEFFICIENTS):
                        e = summed[1][ip, i] if summed[1] is not None else np.nan
                        f.write('%s,%.17g,%.17g,%s,sum(Nc=%g;Nf=%g;NfZ=%s),%.16e,%.16e,%.3e\n'
                                % (regs[ip], p[0], p[1], c, a.Nc, a.Nf, a.NfZ, summed[0][ip, i].real, summed[0][ip, i].imag, e))
    else:
        for ip, p in enumerate(res.points):
            if regs[ip] == 'decay':
                print('# decay: x = %.17g, y = %.17g, z = %.17g' % phys[ip])
            else:
                XY = res.info[ip]['region_point'][1:] if a.physical else p
                print('# region %s: (X, Y) = (%.17g, %.17g) -> x = %.17g, y = %.17g, z = %.17g' % ((regs[ip],) + tuple(XY) + tuple(phys[ip])))
            wd = 40 if res.values_lo is not None else 24
            print('# %-7s %-18s %*s %*s %10s' % ('coeff', 'colour layer', wd, 'Re', wd, 'Im', 'error'))
            for i, c in enumerate(COEFFICIENTS):
                for j, l in enumerate(LAYERS):
                    vr, vi = vstr(ip, i, j)
                    print('  %-7s %-18s %*s %*s %10.2e' % (c, l, wd, vr, wd, vi, err[ip, i, j]))
            if summed is not None:
                print('# summed over the colour layers with Nc = %g, Nf = %g, NfZ = %s' % (a.Nc, a.Nf, a.NfZ))
                for i, c in enumerate(COEFFICIENTS):
                    v = summed[0][ip, i]; e = summed[1][ip, i] if summed[1] is not None else np.nan
                    print('  %-7s %-18s %24.16e %24.16e %10.2e' % (c, 'sum', v.real, v.imag, e))


if __name__ == '__main__':
    main()
