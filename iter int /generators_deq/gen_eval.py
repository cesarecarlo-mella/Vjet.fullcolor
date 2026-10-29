"""evaluate the generators T[1..579] (Zgqq_finite_remainders/Generators.wl) at points of the Euclidean triangle with
tridents_bs + generators_bs.bin.

    from gen_eval import evaluate
    T, info = evaluate([(0.3, 0.2), (0.1, 0.7)])     # T.shape == (2, 579), T[p, i] = T[i+1] at point p
    python3 gen_eval.py x y                           # prints the values
"""
import json, os, subprocess, sys, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'out/generators_bs.bin')
HEAD_Y = os.path.join(HERE, 'out/start_y.head')     # start values at P0 from the y-first regularisation
META = json.load(open(os.path.join(HERE, 'out/meta.json')))
NAMES, WT = META['names'], np.array(META['wt'])


# Path z(t) = P0 + (t + 4 i d t (1 - t)) (P - P0), per coordinate.  The forms are singular on the zeros of the letters
# in use (1-6, 11-16) and, for the odd letters 7, 8, of r (x y z) and of x z + y, x y + z.  The quadratic ones vanish on
# curves crossing the triangle: their complex zeros can lie exactly on the deformed path (e.g. letter 11 at t = 3/4
# for P = (7/30, 19/30) with the default deformation), and then the stepper fails.  Along the path each polynomial has
# degree <= 4 in t, so its zeros are roots; a deformation is used only if they stay away from 0 < t < 1.
P0 = (0.3, 0.3)                                     # start point of out/generators_bs.bin
SING = {1: lambda x, y: x, 2: lambda x, y: y, 3: lambda x, y: 1 - x - y, 4: lambda x, y: 1 - x, 5: lambda x, y: 1 - y,
        6: lambda x, y: x + y, '7': lambda x, y: x * (1 - x - y) + y, '8': lambda x, y: x * y + 1 - x - y,
        11: lambda x, y: x - y + y * y, 12: lambda x, y: -x + x * x + y, 13: lambda x, y: (x - 1) ** 2 - y,
        14: lambda x, y: -y + (x + y) ** 2, 15: lambda x, y: -x + (x + y) ** 2, 16: lambda x, y: -x + (1 - y) ** 2}
# letters that enter only in the production regions (in the region's own variables X, Y; see PROFILES)
SING_REGION = {9: lambda x, y: 1 - y - x * y, 10: lambda x, y: 1 - x - x * y, 17: lambda x, y: x + x * y + y * y,
               18: lambda x, y: x * x + y + x * y, 19: lambda x, y: (y - 1) ** 2 + x * y, 20: lambda x, y: (x - 1) ** 2 + x * y}
DEFORM = [(0.1, 0.2), (0.15, 0.25), (0.25, 0.1), (-0.15, -0.25), (0.2, 0.3), (-0.25, -0.1)]
MARGIN = 1e-3                                       # minimal distance (in t) of a zero from the open path


def path_margin(P, d, region='decay'):
    """smallest distance, in the path parameter, of a complex zero of a letter from the segment 0 < t < 1, and its key;
    zeros at the end point t = 1 itself (P on a curve) are not counted: no deformation can avoid them"""
    ts = np.linspace(0, 1, 5)
    X = P0[0] + (ts + 4j * d[0] * ts * (1 - ts)) * (P[0] - P0[0])
    Y = P0[1] + (ts + 4j * d[1] * ts * (1 - ts)) * (P[1] - P0[1])
    best = (np.inf, None)
    for k, f in profile(region)['sing'].items():
        c = np.polyfit(ts, f(X, Y), 4)
        c = c[np.argmax(np.abs(c) > 1e-14 * np.abs(c).max()):]
        for r in (np.roots(c) if len(c) > 1 else []):
            if abs(r - 1) < 1e-9:
                continue
            dist = abs(r.imag) if 0 <= r.real <= 1 else min(abs(r), abs(r - 1))
            best = min(best, (dist, k), key=lambda b: b[0])
    return best


def choose_deformation(P, second=False, region='decay'):
    """the default deformation if the path keeps MARGIN from every zero, else the one with the largest margin.
    second: a different deformation (the best of the others), for an independent second evaluation"""
    m = [path_margin(P, d, region)[0] for d in DEFORM]
    first = 0 if m[0] >= MARGIN else int(np.argmax(m))
    if not second:
        return DEFORM[first]
    return DEFORM[max((i for i in range(len(DEFORM)) if i != first), key=lambda i: m[i])]


def second_data(data=DATA):
    """the data file with independent start values (y-first regularisation instead of x-first): the same system, only the
    head differs.  Built once from start_y.head and the system part of data, next to data or in the temp directory.
    The start values are the error shared by all evaluations from one file, and near the edges of the triangle the
    rational functions amplify it, so an error estimate must compare two files."""
    head = open(HEAD_Y, 'rb').read()
    for target in (data[:-4] + '_y.bin', os.path.join(tempfile.gettempdir(), 'zgqq_generators_bs_y.bin')):
        try:
            if os.path.isfile(target) and os.path.getsize(target) == os.path.getsize(data) and open(target, 'rb').read(len(head)) == head:
                return target
            with open(data, 'rb') as f:
                main = f.read()
            assert main[:8] == head[:8] and len(main) > len(head)
            with open(target + '.tmp', 'wb') as g:
                g.write(head + main[len(head):])
            os.replace(target + '.tmp', target)
            return target
        except OSError:
            continue
    raise OSError('cannot write the second data file')


# fallback if the stepper still fails at a point
RETRY = [(0.15, 0.25), (0.25, 0.1), (-0.15, -0.25)]

# End points on (or very close to) a spurious curve: letters 11-16 are spurious, the functions (and the amplitude) are
# smooth across their curves, but the DE has poles there (and the rational functions of the remainders have poles on
# 15, 16).  Such a point is replaced by the symmetric average across the curve, along its normal n:
#     f(P) = (15 S1 - 6 S2 + S3) / 10 + O(h^6),   S_k = (f(P + k h n) + f(P - k h n)) / 2.
# All lengths are relative to the local scale rho = min(x, y, z), the distance to the edges (where the functions are
# singular); all spurious curves run into the corners.
SPURIOUS = (11, 12, 13, 14, 15, 16)
NEAR = {11: 1e-8, 12: 1e-8, 13: 1e-8, 14: 1e-8,    # poles of the DE only: the direct integration works until ~exact hits
        15: 1e-5, 16: 1e-5}                         # also poles of the rational functions: digits lost closer than this
H = 1e-3                                            # stencil step, at most rho / 100 (O(h^6) error ~ 36 (h/rho)^6)

# Production regions (region_build.py): the region's triangle in its own variables (X, Y) carries the same 20 letter
# functions, and the region's system is A'_j = sum_k c_kj A_k (out/region_<A>.json: letter_map[k][j] = c_kj).  Its
# letters are the j reached from the decay letters in use; its spurious curves are those among 11-16, and the rational
# functions (decay letters 15, 16 in their denominators) have poles on the images of 15 and 16.
_PROFILES = {}


def profile(region='decay'):
    """letters (for the path check), spurious letters and their NEAR thresholds of the decay triangle or a region"""
    if region not in _PROFILES:
        if region == 'decay':
            _PROFILES[region] = dict(sing=SING, spurious=SPURIOUS, near=NEAR)
        else:
            cmap = json.load(open(os.path.join(HERE, 'out', 'region_%s.json' % region)))['letter_map']
            used = [1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16]
            img = lambda k: {int(j) for j, c in cmap[str(k)].items() if c}
            lets = set().union(*[img(k) for k in used])
            pole = img(15) | img(16)
            allf = {**SING, **SING_REGION}
            sing = {(str(j) if j in (7, 8) else j): allf[str(j) if j in (7, 8) else j] for j in sorted(lets)}
            spur = tuple(j for j in sorted(lets) if 11 <= j <= 16)
            _PROFILES[region] = dict(sing=sing, spurious=spur, near={j: 1e-5 if j in pole else 1e-8 for j in spur})
    return _PROFILES[region]


def spurious_near(P, region='decay'):
    """(distance, unit normal, letter) of the closest spurious curve if it is closer than NEAR[letter] * rho, else None"""
    x, y = P
    rho = min(x, y, 1 - x - y)
    best, e = None, 1e-7 * rho
    pr = profile(region)
    for k in pr['spurious']:
        f = pr['sing'][k]
        g = np.array([(f(x + e, y) - f(x - e, y)) / (2 * e), (f(x, y + e) - f(x, y - e)) / (2 * e)])
        d = abs(f(x, y)) / np.linalg.norm(g)
        if d < pr['near'][k] * rho and (best is None or d < best[0]):
            best = (d, g / np.linalg.norm(g), k)
    return best


def plan(points, hfac=1.0, region='decay'):
    """the points to integrate: each point itself, or the 6 stencil points if it is near a spurious curve
    (step hfac * min(H, rho / 100)).  Returns (flat list of points, [(start in flat, count, averaging info or None)])"""
    flat, pl = [], []
    for P in points:
        near = spurious_near(P, region)
        if near is None:
            pl.append((len(flat), 1, None)); flat.append(P)
            continue
        d, n, k = near
        h = hfac * min(H, min(P[0], P[1], 1 - P[0] - P[1]) / 100)
        pl.append((len(flat), 6, dict(letter=k, distance=d, h=h)))
        flat += [(P[0] + s * j * h * n[0], P[1] + s * j * h * n[1]) for j in (1, 2, 3) for s in (1, -1)]
    return flat, pl


def collapse(values, infos, pl):
    """values on the flat points -> values on the original points (stencils combined)"""
    out, inf = [], []
    for i0, c, av in pl:
        if c == 1:
            out.append(values[i0]); inf.append(infos[i0])
            continue
        V = values[i0:i0 + 6]
        S = [(V[2 * j] + V[2 * j + 1]) / 2 for j in range(3)]
        out.append((15 * S[0] - 6 * S[1] + S[2]) / 10)
        inf.append(dict(infos[i0], time=sum(i['time'] for i in infos[i0:i0 + 6]), averaged=av))
    return np.array(out), inf


def _run(points, tol, binary, data, delta=None):
    """tridents_bs on the points: ([(info, values)] of the points done, whether it stopped on a step-size failure)"""
    if not os.path.isfile(binary):
        raise FileNotFoundError('%s not found: build it first (make)' % binary)
    inp = '\n'.join('%.17g %.17g' % tuple(p) for p in points)
    args = [binary, data, str(tol), str(tol)] + (['%.17g' % delta[0], '%.17g' % delta[1]] if delta else [])
    r = subprocess.run(args, input=inp, capture_output=True, text=True)
    failed = r.returncode != 0
    if failed and 'Max number of iterations' not in r.stderr:
        raise RuntimeError('tridents_bs: ' + r.stderr)
    done, cur = [], None
    for l in r.stdout.split('\n'):
        if l.startswith('#'):
            h = l.split(); cur = []
            done.append((dict(x=float(h[1]), y=float(h[2]), steps=int(h[3]), evals=int(h[4]), time=float(h[5]),
                              delta=delta or (0.1, 0.2)), cur))
        elif l:
            a, b = l.split(); cur.append(complex(float(a), float(b)))
    return done, failed


def _run_all(points, tol, binary, data, delta):
    """all points with one deformation; a point where the stepper fails is redone with the RETRY deformations"""
    res, i = [], 0
    while i < len(points):
        done, failed = _run(points[i:], tol, binary, data, delta)
        res += done; i += len(done)
        if failed:
            for d in RETRY:
                one, f1 = _run([points[i]], tol, binary, data, d)
                if not f1:
                    break
            else:
                P = points[i]
                on = [k for k, f in SING.items() if abs(f(*P)) < 1e-12]
                raise RuntimeError('no path deformation works at %s%s' % (P, ' (end point on the zero curve of %s)' % on if on else ''))
            res += one; i += 1
    return res


def evaluate(points, tol=1e-10, binary=os.path.join(HERE, 'cpp/tridents_bs'), data=DATA):
    """generators at the points (points near a spurious curve: averaged across it, see plan)"""
    flat, pl = plan([tuple(p) for p in points])
    T, infos = evaluate_direct(flat, tol, binary, data)
    return collapse(T, infos, pl)


def evaluate_direct(points, tol=1e-10, binary=os.path.join(HERE, 'cpp/tridents_bs'), data=DATA,
                    second=False, region='decay'):
    """generators at the points themselves (path deformation chosen per point; second: on a second, different path)"""
    points = [tuple(p) for p in points]
    groups = {}
    for j, P in enumerate(points):
        groups.setdefault(choose_deformation(P, second, region), []).append(j)
    res = [None] * len(points)
    for d, idx in groups.items():
        for j, r in zip(idx, _run_all([points[j] for j in idx], tol, binary, data, None if d == DEFORM[0] else d)):
            res[j] = r
    NT = len(WT)
    return (np.array([np.array(v).reshape(7, NT)[WT, np.arange(NT)] for _, v in res]),    # generator i in the slot of its weight
            [inf for inf, _ in res])


if __name__ == '__main__':
    T, info = evaluate([(float(sys.argv[1]), float(sys.argv[2]))])
    for i, v in enumerate(T[0]):
        print('%s (weight %d) = %.15e %+.15e i' % (NAMES[i], WT[i], v.real, v.imag))
    print(info[0])
