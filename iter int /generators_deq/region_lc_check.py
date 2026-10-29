"""the remainders in region y > 1 against the published planar results continued to production kinematics
(2024 ancillary, FiniteRemainders/lmL.m: alpha[i,3]; lpL.m: gamma[i,3]; q bar g -> V q bar; variables v = 1/y,
u = -z/y = our triangle of region y with Z = v, Y = u).  Their three layers Nc^3, Nc^2 Nf, Nc Nf^2 are separated by
evaluating at several (Nc, Nf); every G is evaluated by GiNaC.
usage: region_lc_check.py [points u,v;u,v;...] [--json FILE]   (FILE: the published values, as a test reference)"""
import os, re, subprocess, sys
import mpmath as mp
import numpy as np
import remainder_eval as R

HERE = os.path.dirname(os.path.abspath(__file__))
ANC = os.path.join(HERE, '..', '..', 'reproduceLC', '2024_Vjet_3L_LC_ancillary', 'FiniteRemainders')
mp.mp.dps = 30
GRE = re.compile(r'(?<![A-Za-z])G\[[^\[\]]*\]')


def split_top(s):
    out, d, cur = [], 0, ''
    for ch in s:
        if ch in '[({': d += 1
        elif ch in '])}': d -= 1
        if ch == ',' and d == 0:
            out.append(cur); cur = ''
        else:
            cur += ch
    out.append(cur)
    return out


def pyexpr(e):
    e = e.replace('\n', ' ').replace('\\', '').replace('^', '**')
    e = e.replace('Pi', 'PI').replace('Zeta[3]', 'Z3').replace('Zeta[5]', 'Z5')
    e = re.sub(r'\bI\b', 'IU', e)
    return re.sub(r'(?<![\w.])(\d+)(?![\w.])', r'Q(\1)', e)


def parse_c(v):
    v = v.replace(' ', '')
    if not v.endswith('I'):
        return mp.mpc(mp.mpf(v), 0)
    body = v[:-1].rstrip('*'); cut = -1
    for k in range(len(body) - 1, 0, -1):
        if body[k] in '+-' and body[k - 1] not in 'eE':
            cut = k; break
    re_, im_ = ('0', body) if cut < 0 else (body[:cut], body[cut:])
    im_ = {'': '1', '+': '1', '-': '-1'}.get(im_, im_)
    return mp.mpc(mp.mpf(re_), mp.mpf(im_))


def load(fn):
    s = re.sub(r'\(\*.*?\*\)', '', open(os.path.join(ANC, fn)).read(), flags=re.S).strip()[1:-1]
    pub = {}
    for r in split_top(s):
        lhs, rhs = r.split('->', 1)
        pub[lhs.strip()] = rhs.strip()
    atoms = sorted({g for e in pub.values() for g in GRE.findall(e)}, key=len, reverse=True)
    code = {}
    for k, e in pub.items():
        for n, g in enumerate(atoms):
            e = e.replace(g, ' GG%d ' % n)
        code[k] = compile(pyexpr(e), k, 'eval')
    return code, atoms


def ginsh_values(atoms, u, v, digits=25):
    lines = ['Digits=%d;' % digits, 'u=%s:' % u, 'v=%s:' % v]
    for g in atoms:
        a = split_top(g[2:-1])
        lines.append('evalf(G({%s},%s));' % (','.join(t.strip() for t in a[:-1]), a[-1].strip()))
    r = subprocess.run(['ginsh'], input='\n'.join(lines) + '\n', capture_output=True, text=True)
    out = [o.strip() for o in r.stdout.split('\n') if o.strip()][-len(atoms):]
    assert len(out) == len(atoms) and not r.stderr.strip(), r.stderr[:300]
    return [parse_c(x) for x in out]


JSON = sys.argv[sys.argv.index('--json') + 1] if '--json' in sys.argv else None
args = [a for i, a in enumerate(sys.argv[1:], 1) if a != '--json' and sys.argv[i - 1] != '--json']
pts = [tuple(float(a) for a in p.split(',')) for p in args[0].split(';')] if args else \
      [(0.2, 0.3), (0.3, 0.2), (0.15, 0.6), (0.4, 0.4)]
LAY = {3: 'Nc^3*Nf^0*NfZ^0', 2: 'Nc^2*Nf^1*NfZ^0', 1: 'Nc^1*Nf^2*NfZ^0'}
ours, info = R.evaluate_remainders_region('y', [(1 - u - v, u) for u, v in pts])
worst = 0
published = []
for fn, coef, names in (('lmL.m', 'alpha', ('alpha1', 'alpha2', 'alpha3')), ('lpL.m', 'gamma', ('gamma1', 'gamma2', 'gamma3'))):
    code, atoms = load(fn)
    print('%s: %d G atoms' % (fn, len(atoms)), flush=True)
    for ip, (u, v) in enumerate(pts):
        env = {'Q': mp.mpf, 'PI': mp.pi, 'Z3': mp.zeta(3), 'Z5': mp.zeta(5), 'IU': mp.mpc(0, 1), 'u': mp.mpf(u), 'v': mp.mpf(v), 'q2': mp.mpf(1)}
        env.update({'GG%d' % n: g for n, g in enumerate(ginsh_values(atoms, u, v))})
        for i, nm in enumerate(names, 1):
            E = lambda Nc, Nf: complex(eval(code['%s[%d, 3]' % (coef, i)], dict(env, Nc=mp.mpf(Nc), Nf=mp.mpf(Nf))))
            a = E(1, 0); p1 = E(1, 1); m1 = E(1, -1)
            b = (p1 - m1) / 2; c = (p1 + m1) / 2 - a
            chk = abs(E(2, 1) - (8 * a + 4 * b + 2 * c)) / max(abs(E(2, 1)), 1)
            for nc, val in ((3, a), (2, b), (1, c)):
                j = R.KEYS.index((nm, LAY[nc]))
                o = ours[ip, j]
                rel = abs(o - val) / abs(val)
                worst = max(worst, rel)
                published.append(dict(point=ip, coefficient=nm, layer=LAY[nc], value=[val.real, val.imag]))
                print('   (u, v) = (%g, %g) %s %-17s ours %+.10e %+.10ei  published %+.10e %+.10ei  rel diff %.1e%s'
                      % (u, v, nm, LAY[nc], o.real, o.imag, val.real, val.imag, rel, '' if chk < 1e-12 else '  (!! layer split %.0e)' % chk), flush=True)
print('worst relative difference: %.1e' % worst)
if JSON:
    import json
    json.dump(dict(source='reproduceLC/2024_Vjet_3L_LC_ancillary/FiniteRemainders lmL.m (alpha[i,3]), lpL.m (gamma[i,3]), '
                          'evaluated with GiNaC at 25 digits; layers from evaluations at several (Nc, Nf)',
                   region='y', points=[[1 - u - v, u] for u, v in pts], uv=[list(p) for p in pts], values=published),
              open(JSON, 'w'), indent=1)
    print('written', JSON)
