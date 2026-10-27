import re, glob, os, sys, random, collections, json, itertools, time, pickle
from functools import lru_cache
P=1048573
def inv(x): return pow(x,P-2,P)
DATA='/mnt/user-data/uploads/denominator_bracket'
CACHE='/home/claude/fit/rows.pkl'
if os.path.exists(CACHE): by=pickle.load(open(CACHE,'rb'))
else:
    by=collections.defaultdict(list)
    for f in sorted(glob.glob(DATA+'/*.wl')):
        amp=os.path.basename(f)[:-3]; L=open(f).read().splitlines()
        i=next(k for k,l in enumerate(L) if l.startswith('+eps^(0)*')); buf=[]
        for l in L[i+1:]:
            m=re.match(r'^\+\(([-0-9/]+)\)\*\((.*)\)\*\((.*)\)$',l)
            if m:
                c=m.group(3); w=re.search(r'II\[([0-9,]*)\]',c).group(1)
                a,_,b=m.group(1).partition('/')
                buf.append((c[:c.index('II[')].rstrip('*'),tuple(int(x) for x in w.split(',')) if w else (),m.group(2),int(a)*inv(int(b) if b else 1)%P))
            elif l.startswith(')/('):
                by[amp]+=[(c,w,p,l[3:-1],co) for c,w,p,co in buf]; buf=[]
    by=dict(by); pickle.dump(by,open(CACHE,'wb'))
@lru_cache(None)
def rho(w):
    n=len(w)
    if n==1: return ((w,1),)
    out=collections.Counter(); f=(n-1)*inv(n)%P
    for x,c in rho(w[:-1]): out[x+(w[-1],)]+=c*f
    for x,c in rho(w[1:]): out[x+(w[0],)]-=c*f
    return tuple((k,v%P) for k,v in out.items() if v%P)
def rhoS(S):
    out=collections.Counter()
    for w,c in S.items():
        for y,d in rho(w): out[y]+=c*d
    return {k:v%P for k,v in out.items() if v%P}
def delta(F,a,b):
    Z=collections.defaultdict(collections.Counter)
    for w,c in F.items(): Z[w[:a]][w[a:]]+=c; Z[w[b:]][w[:b]]-=c
    out=collections.Counter()
    for Lw,R in Z.items():
        RR=collections.Counter()
        for r,c in R.items():
            for y,d in rho(r): RR[y]+=c*d
        RR={y:v%P for y,v in RR.items() if v%P}
        if not RR: continue
        for x,c in rho(Lw):
            for y,d in RR.items(): out[(x,y)]+=c*d
    return {k:v%P for k,v in out.items() if v%P}
class Span:
    def __init__(s): s.b={}
    def reduce(s,v):
        v=dict(v); b=s.b
        while v:
            piv=min(v)
            if piv not in b: return v
            r=b[piv]; f=v[piv]
            for k,c in r.items():
                x=(v.get(k,0)-f*c)%P
                if x: v[k]=x
                else: v.pop(k,None)
        return v
    def add(s,v):
        v=s.reduce(v)
        if not v: return 0
        piv=min(v); f=inv(v[piv]); s.b[piv]={k:c*f%P for k,c in v.items()}; return 1
    def contains(s,v): return not s.reduce(v)
    def copy(s): n=Span(); n.b={k:dict(v) for k,v in s.b.items()}; return n
    def __len__(s): return len(s.b)
def data_vectors(const,W,amps=None,even_only=True,maxpts=12,seed=3):
    rnd=random.Random(seed); out=[]
    for amp in sorted(by):
        if amps and amp not in amps: continue
        g=collections.defaultdict(collections.Counter)
        for c,wd,p,d,co in by[amp]:
            if c==const and len(wd)==W and (not even_only or not {7,8}&set(wd)): g[(p,d)][wd]+=co
        if not g: continue
        G=[(compile(p.replace('^','**'),'p','eval'),compile(d.replace('^','**'),'d','eval'),{w:c%P for w,c in v.items()}) for (p,d),v in g.items()]
        loc=Span(); miss=0
        while miss<3 and len(loc)<maxpts:
            env={'t':rnd.randrange(P),'u':rnd.randrange(P)}; F=collections.Counter()
            for pc,dc,vec in G:
                r=eval(pc,env)%P*inv(eval(dc,env)%P)%P
                for w,c in vec.items(): F[w]+=r*c
            v=rhoS({w:c%P for w,c in F.items() if c%P})
            if loc.add(v): out.append((amp,v)); miss=0
            else: miss+=1
    return out
# ---- candidate functions
C=json.load(open('/home/claude/fit/pairs.json')); LET=C['letters']
def rhoA(S):
    out=collections.Counter()
    for w,c in S:
        for y,d in rho(tuple(w)): out[y]+=c*d
    return {k:v%P for k,v in out.items() if v%P}
def subst(rA,vecs):
    out=collections.Counter()
    sp={a:[(LET[i],e) for i,e in enumerate(v) if e] for a,v in vecs.items()}
    for w,c in rA.items():
        for combo in itertools.product(*[sp[a] for a in w]):
            m=c
            for _,e in combo: m*=e
            out[tuple(k for k,_ in combo)]+=m
    return {k:v%P for k,v in out.items() if v%P}
def classical(W):
    rA=rhoA(C['symbols'][f'Li{W}'])
    return [(('Li%d'%W,a),subst(rA,{0:a[1],2:b[1]})) for a,b in C['args']]
def depth2(fn,pairs):
    rA=rhoA(C['symbols'][fn]); out=[]
    for p in pairs:
        vecs={0:p['x1'][1],1:p['x2'][1],2:p['1-x1'][1],3:p['1-x2'][1],4:p['1-x1x2'][1]}
        out.append(((fn,p['x1'],p['x2']),subst(rA,vecs)))
    return out
TAG=10**9
def support(S_cl, cand, data):
    """S_cl: Span of classical; cand: list of (key,vec) independent mod S_cl; returns set of cand indices used and uncovered count"""
    B=S_cl.copy()
    for j,(k,v) in enumerate(cand):
        vv=dict(v); vv[(TAG,j)]=1; B.add(vv)
    used=set(); unc=0
    for a,v in data:
        r=B.reduce(v)
        if any(k[0]!=TAG for k in r): unc+=1
        used|={k[1] for k in r if k[0]==TAG}
    return used,unc
