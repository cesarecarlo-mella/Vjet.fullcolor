"""
charts.py -- the expansion charts of the blow-up recursion: one registry read by every module
(letters, recursion, file names, map back to physical variables).

Every chart is SMALL = t * BIG^2 around a corner where SMALL, BIG -> 0; the series are in t, BIG,
Log[t], Log[BIG], and phys.py maps them back to SMALL, BIG.

Decay corners (Euclidean triangle, x, y, z > 0, z = 1 - x - y); id 1..6, file tag c<id>:

    1: x = t y^2   2: z = t y^2   3: x = t z^2          4-6: the mirrors of 1-3 (SMALL <-> BIG),
                                                         for the crossed families; same constants

Continued regions A > 1 (A one of x, y, z; the other two are negative).  Triangle coordinates

    v = 1/A ,   bh = -b/A  for each of the other two variables b ,     v + bh + ch = 1 ,

(symbols v, xh, yh, zh).  At the high-energy corner bh, v -> 0 (ch -> 1) there are two charts,
with DIFFERENT constants (from weight 2 the masters depend on the ratio bh/v there):

    R<A>_<b>v : bh = t v^2   (BIG v)    bh << v     the 's = 0' chart  (CC2, CC3 of cc_series)
    R<A>_v<b> : v = t bh^2   (BIG bh)   v << bh     the 's = 1' chart  (CC2e, CC3e)

e.g. RZ_yv is the old CC2 (z > 1: yh = u, xh = w), RZ_vy the old CC2e.  The threshold corner of
region A (bh, ch -> 0) is the decay corner b, c -> 0 continued: it needs no chart of its own.

In the continued region x y z > 0 and Sqrt[x] Sqrt[y] Sqrt[z] = (I Sqrt|b|)(I Sqrt|c|) Sqrt[A]
= -Sqrt|xyz| with the continuation Log b -> log|b| + I Pi: ODD_SIGN = -1 for the odd letters l7, l8.
"""

DECAY = {
    1: ('x', 'y'), 2: ('z', 'y'), 3: ('x', 'z'),        # (SMALL, BIG)
    4: ('y', 'x'), 5: ('y', 'z'), 6: ('z', 'x'),
}
CONT = {}
for _A in 'xyz':
    for _b in 'xyz':
        if _b != _A:
            CONT['R%s_%sv' % (_A.upper(), _b)] = dict(region=_A, small=_b + 'h', big='v')
            CONT['R%s_v%s' % (_A.upper(), _b)] = dict(region=_A, small='v', big=_b + 'h')
CHART_IDS = list(DECAY) + sorted(CONT)


def parse(s):
    """'1' -> 1, 'RZ_yv' -> 'RZ_yv' (command line); raises on unknown charts."""
    s = str(s).strip()
    c = int(s) if s.isdigit() else s
    if c not in DECAY and c not in CONT:
        raise ValueError("unknown chart %r (known: %s)" % (s, CHART_IDS))
    return c


def parse_list(s):
    return [parse(p) for p in str(s).split(',') if p.strip()]


def tag(c):
    """file tag: c1..c6, or the chart name."""
    return 'c%d' % c if c in DECAY else c


def small_big(c):
    """(SMALL, BIG) symbol names."""
    if c in DECAY:
        return DECAY[c]
    return CONT[c]['small'], CONT[c]['big']


def big(c):
    return small_big(c)[1]


def is_continued(c):
    return c in CONT


def bdy_tag(c):
    """tag of the boundary-constant files: a decay mirror uses its standard corner."""
    return tag({4: 1, 5: 2, 6: 3}.get(c, c))


def odd_sign(c):
    return -1 if c in CONT else 1


def blowup(c):
    sm, bg = small_big(c)
    return '%s = t*%s^2' % (sm, bg)


def param(c, t, B):
    """(x, y) of the chart as expressions in t, B (any objects with + - * /, e.g. sympy symbols)."""
    if c in DECAY:
        sm, bg = DECAY[c]
        val = {sm: t * B**2, bg: B}
        rest = ({'x', 'y', 'z'} - {sm, bg}).pop()
        val[rest] = 1 - B - t * B**2
        return val['x'], val['y']
    R = CONT[c]
    A = R['region']
    others = [b for b in 'xyz' if b != A]
    coord = {R['small']: t * B**2, R['big']: B}
    rest = ({'v'} | {b + 'h' for b in others}) - set(coord)
    coord[rest.pop()] = 1 - B - t * B**2
    val = {A: 1 / coord['v']}
    for b in others:
        val[b] = -coord[b + 'h'] / coord['v']
    return val['x'], val['y']
