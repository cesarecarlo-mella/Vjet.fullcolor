#!/usr/bin/env python3
"""Compare the GPL form of the leading-colour HA coefficients (HA_GPL_leading/*_G.m, variables x,y)
with the published finite remainders 2024_Vjet_3L_LC_ancillary/FiniteRemainders/decay.m
(variables y, z = 1-x-y, q2 -> 1).  All G functions are evaluated with GiNaC (ginsh, 30 digits);
the rational prefactors are evaluated exactly-ish in mpmath (40 digits)."""
import re, sys, subprocess, itertools
import mpmath as mp
mp.mp.dps = 40
DEC = sys.argv[1]; HADIR = sys.argv[2]
PTS = [(mp.mpf(3)/10, mp.mpf(2)/11), (mp.mpf(1)/7, mp.mpf(1)/5), (mp.mpf(1)/2, mp.mpf(1)/4)]
PTS_Q = ['3/10,2/11', '1/7,1/5', '1/2,1/4']
LAYERS = {'C3F0Z0': (3, 0), 'C2F1Z0': (2, 1), 'C1F2Z0': (1, 2)}
HANAME = {1: 'alpha[1, 3]', 2: 'alpha[2, 3]', 3: 'alpha[3, 3]', 4: 'gamma[1, 3]', 5: 'gamma[2, 3]', 6: 'gamma[3, 3]'}
GRE = re.compile(r'G\[[^\[\]]*\]')

def split_top(s):
    out = []; d = 0; cur = ''
    for ch in s:
        if ch in '[({': d += 1
        elif ch in '])}': d -= 1
        if ch == ',' and d == 0: out.append(cur); cur = ''
        else: cur += ch
    out.append(cur); return out

def pyexpr(e):
    e = e.replace('\n', ' ').replace('\\', '')
    e = e.replace('^', '**').replace('Pi', 'PI').replace('Zeta[3]', 'Z3').replace('Zeta[5]', 'Z5')
    e = re.sub(r'\bI\b', 'IU', e)
    e = re.sub(r'(?<![\w.])(\d+)(?![\w.])', r'Q(\1)', e)
    return e

def ginsh_values(gatoms, xq, yq):
    """gatoms: list of Mathematica G strings in x,y,z (z = 1-x-y). Returns dict -> mpc."""
    lines = ['Digits=32;', 'x=%s;' % xq, 'y=%s;' % yq, 'z=1-x-y;']
    for g in gatoms:
        args = split_top(g[2:-1])
        idx, var = args[:-1], args[-1].strip()
        lines.append('evalf(G({%s},%s));' % (','.join(a.strip() for a in idx), var))
    out = subprocess.run(['ginsh'], input='\n'.join(lines) + '\n', capture_output=True, text=True).stdout.split('\n')
    out = [o.strip() for o in out if o.strip()]
    vals = out[-len(gatoms):] if gatoms else []
    assert len(vals) == len(gatoms), (len(vals), len(gatoms), out[:5])
    res = {}
    for g, v in zip(gatoms, vals):
        res[g] = parse_c(v.replace(' ', ''))
    return res

def parse_c(v):
    if not v.endswith('I'):
        return mp.mpc(mp.mpf(v), 0)
    body = v[:-1].rstrip('*')
    cut = -1
    for k in range(len(body) - 1, 0, -1):
        if body[k] in '+-' and body[k - 1] not in 'eE':
            cut = k; break
    if cut < 0:
        re_, im_ = '0', body
    else:
        re_, im_ = body[:cut], body[cut:]
    if im_ in ('', '+'): im_ = '1'
    if im_ == '-': im_ = '-1'
    return mp.mpc(mp.mpf(re_), mp.mpf(im_))

def evaluate(expr_py, gmap, env):
    e = expr_py
    for i, g in enumerate(gmap):
        pass
    return eval(e, env)

def prep(expr):
    atoms = sorted(set(GRE.findall(expr)), key=len, reverse=True)
    names = {}
    for i, g in enumerate(atoms):
        names[g] = 'GG%d' % i
    e = expr
    # replace atoms (longest first)
    for g in atoms:
        e = e.replace(g, ' ' + names[g] + ' ')
    return atoms, names, compile(pyexpr(e), 'x', 'eval')

# ---- load published
s = open(DEC).read()
s = re.sub(r'\(\*.*?\*\)', '', s, flags=re.S).strip()[1:-1]
pub = {}
for r in split_top(s):
    lhs, rhs = r.split('->', 1); pub[lhs.strip()] = rhs.strip()
# ---- load mine
mine = {}
for L in LAYERS:
    for i in range(1, 7):
        t = open('%s/%s_HA%d_G.m' % (HADIR, L, i)).read()
        t = re.sub(r'\(\*.*?\*\)', '', t, flags=re.S)
        mine[(L, i)] = t.split('=', 1)[1].strip().rstrip(';')
prepped = {k: prep(v) for k, v in list(pub.items()) + [((k[0], k[1]), v) for k, v in mine.items()]}
allatoms = sorted(set(a for k, (atoms, _, _) in prepped.items() for a in atoms))
print('distinct G atoms: %d' % len(allatoms), flush=True)
base = {'Q': mp.mpf, 'PI': mp.pi, 'Z3': mp.zeta(3), 'Z5': mp.zeta(5), 'IU': mp.mpc(0, 1)}
worst = {}
for (xv, yv), q in zip(PTS, PTS_Q):
    xq, yq = q.split(',')
    gv = ginsh_values(allatoms, xq, yq)
    print('\n=== point x=%s, y=%s (z = 1-x-y)' % (xq, yq), flush=True)
    for i in range(1, 7):
        name = HANAME[i]
        atoms, names, code = prepped[name]
        def pubval(Nc, Nf):
            env = dict(base, x=xv, y=yv, z=1 - xv - yv, q2=mp.mpf(1), Nc=mp.mpf(Nc), Nf=mp.mpf(Nf))
            env.update({names[g]: gv[g] for g in atoms})
            return eval(code, env)
        # colour coefficients: the published expression is a polynomial in Nc, Nf of degree 3
        mons = [(a, b) for a in range(0, 4) for b in range(0, 4) if a + b <= 3]
        ptsC = [(2,1),(3,4),(5,2),(7,9),(11,3),(13,5),(4,7),(6,11),(9,1),(8,6),(10,13),(12,2)]
        A = mp.matrix([[mp.mpf(Nc) ** a * mp.mpf(Nf) ** b for (a, b) in mons] for (Nc, Nf) in ptsC])
        vals = [pubval(Nc, Nf) for (Nc, Nf) in ptsC]
        # least squares on complex: solve normal equations separately for re/im
        coef = {}
        for part in (mp.re, mp.im):
            sol = mp.lu_solve(A.T * A, A.T * mp.matrix([part(v) for v in vals]))
            for k, m in enumerate(mons):
                coef.setdefault(m, [0, 0])
                coef[m][0 if part is mp.re else 1] = sol[k]
        resid = max(abs(sum(complex(coef[m][0] + 1j * coef[m][1]) * Nc ** m[0] * Nf ** m[1] for m in mons) - complex(v)) for (Nc, Nf), v in zip(ptsC, vals))
        present = [m for m in mons if abs(complex(coef[m][0] + 1j * coef[m][1])) > 1e-20]
        for L, (a, b) in LAYERS.items():
            atoms2, names2, code2 = prepped[(L, i)]
            env = dict(base, x=xv, y=yv)
            env.update({names2[g]: gv[g] for g in atoms2})
            mv = eval(code2, env)
            pv = mp.mpc(coef[(a, b)][0], coef[(a, b)][1])
            d = abs(mv - pv); rel = d / max(abs(pv), mp.mpf(10) ** -30)
            worst[(L, i)] = max(worst.get((L, i), 0), float(rel))
            print('  HA%d=%-12s %s  mine %s  published %s  |diff| %.1e  rel %.1e' % (
                i, name, L, mp.nstr(mv, 15), mp.nstr(pv, 15), float(d), float(rel)), flush=True)
        print('     (colour fit residual %.1e; monomials present: %s)' % (resid, ['Nc^%d Nf^%d' % m for m in present]))
print('\nSUMMARY  worst relative difference over the 3 points:')
for k in sorted(worst): print('  %s HA%d  %.1e' % (k[0], k[1], worst[k]))
