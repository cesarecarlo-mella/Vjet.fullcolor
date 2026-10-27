#!/usr/bin/env python3
"""check_series.py <G-file> <series-file> : value of the series vs the exact G form (GiNaC) at small points"""
import re, sys, subprocess
import mpmath as mp
sys.path.insert(0, '/home/claude/cmpLC')
mp.mp.dps = 40
gfile, sfile = sys.argv[1], sys.argv[2]
PTS = sys.argv[3:] or ['1/1000,1/50', '1/100,1/20', '1/50,1/10']
def split_top(s):
    out = []; d = 0; cur = ''
    for ch in s:
        if ch in '[({': d += 1
        elif ch in '])}': d -= 1
        if ch == ',' and d == 0: out.append(cur); cur = ''
        else: cur += ch
    out.append(cur); return out
def parse_c(v):
    if not v.endswith('I'): return mp.mpc(mp.mpf(v), 0)
    body = v[:-1].rstrip('*'); cut = -1
    for k in range(len(body) - 1, 0, -1):
        if body[k] in '+-' and body[k - 1] not in 'eE': cut = k; break
    re_, im_ = ('0', body) if cut < 0 else (body[:cut], body[cut:])
    if im_ in ('', '+'): im_ = '1'
    if im_ == '-': im_ = '-1'
    return mp.mpc(mp.mpf(re_), mp.mpf(im_))
def pyexpr(e):
    e = e.replace('\n', ' ').replace('^', '**').replace('Pi', 'PI').replace('Zeta[3]', 'Z3').replace('Zeta[5]', 'Z5')
    e = e.replace('Log[x]', 'LX').replace('Log[y]', 'LY')
    e = re.sub(r'\bI\b', 'IU', e)
    return re.sub(r'(?<![\w.])(\d+)(?![\w.])', r'Q(\1)', e)
def body(f):
    b = re.sub(r'\(\*.*?\*\)', '', open(f).read(), flags=re.S).split('=', 1)[1].strip().rstrip(';').strip()
    if b.startswith('(') and b.endswith(')'): b = b[1:-1].strip()
    return b
gexpr = body(gfile); sexpr = body(sfile)
atoms = sorted(set(re.findall(r'G\[[^\[\]]*\]', gexpr)), key=len, reverse=True)
ge = gexpr
for i, g in enumerate(atoms): ge = ge.replace(g, ' GG%d ' % i)
gcode = compile(pyexpr(ge), 'g', 'eval'); scodes = [compile(pyexpr(t), 's', 'eval') for t in sexpr.split('\n + ')]
base = {'Q': mp.mpf, 'PI': mp.pi, 'Z3': mp.zeta(3), 'Z5': mp.zeta(5), 'IU': mp.mpc(0, 1)}
for pt in PTS:
    xq, yq = pt.split(',')
    lines = ['Digits=40;', 'x=%s;' % xq, 'y=%s;' % yq]
    for g in atoms:
        a = split_top(g[2:-1]); lines.append('evalf(G({%s},%s));' % (','.join(t.strip() for t in a[:-1]), a[-1].strip()))
    out = [o.strip() for o in subprocess.run(['ginsh'], input='\n'.join(lines) + '\n', capture_output=True, text=True).stdout.split('\n') if o.strip()]
    vals = out[-len(atoms):]
    xv = mp.mpf(int(xq.split('/')[0])) / int(xq.split('/')[1]); yv = mp.mpf(int(yq.split('/')[0])) / int(yq.split('/')[1])
    env = dict(base, x=xv, y=yv, LX=mp.log(xv), LY=mp.log(yv)); env.update({'GG%d' % i: parse_c(v) for i, v in enumerate(vals)})
    ex = eval(gcode, env); se = mp.fsum([eval(c, env) for c in scodes])
    print('  (x,y)=(%s,%s)  exact %s  series %s  rel.diff %.1e' % (xq, yq, mp.nstr(ex, 18), mp.nstr(se, 18), float(abs(ex - se) / abs(ex))), flush=True)
