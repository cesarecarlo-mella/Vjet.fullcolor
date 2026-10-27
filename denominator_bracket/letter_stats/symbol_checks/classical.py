import pickle,collections,os,sys,random,json,time
from functools import lru_cache
P=(1<<61)-1
inv=lambda x:pow(x,P-2,P)
rows=pickle.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'rows.pkl'),'rb'))
Q=set(range(9,17)); ALG={7,8}
by=collections.defaultdict(list)
for r in rows: by[r[0]].append(r)
@lru_cache(None)
def rho(w):
    n=len(w)
    if n==1: return {w:1}
    out=collections.Counter(); f=(n-1)*inv(n)%P
    for x,c in rho(w[:-1]).items(): out[x+(w[-1],)]+=c*f
    for x,c in rho(w[1:]).items(): out[x+(w[0],)]-=c*f
    return {k:v%P for k,v in out.items() if v%P}
def ev(expr,t,u,cache={}):
    k=(expr,t,u)
    if k not in cache: cache[k]=eval(expr.replace('^','**'),{'t':t,'u':u})%P
    return cache[k]
def evc(s):
    a,_,b=s.partition('/'); return int(a)*inv(int(b) if b else 1)%P
def pure(amp,t,u):
    """pure functions F_P: dict (const, weight) -> {word: coef mod P}"""
    F=collections.defaultdict(collections.Counter)
    for a,c,w,d,co,p in by[amp]:
        F[(c,len(w))][w]+=evc(co)*ev(p,t,u)*inv(ev(d,t,u))
    return {k:{w:x%P for w,x in v.items() if x%P} for k,v in F.items()}
def delta(Fw,a,b):
    Z=collections.defaultdict(collections.Counter)
    for w,c in Fw.items():
        Z[w[:a]][w[a:]]+=c; Z[w[b:]][w[:b]]-=c
    out=collections.Counter()
    for L,R in Z.items():
        RR=collections.Counter()
        for r,c in R.items():
            if c%P:
                for y,d in rho(r).items(): RR[y]+=c*d
        RR={y:v%P for y,v in RR.items() if v%P}
        if not RR: continue
        for x,c in rho(L).items():
            for y,d in RR.items(): out[(x,y)]+=c*d
    return {k:v%P for k,v in out.items() if v%P}
def cls(k):
    s=set(k[0])|set(k[1])
    return ('Q' if s&Q else '')+('A' if s&ALG else '') or 'lin'
TESTS=[((('',4)),[(2,2)]),((('',5)),[(3,2)]),((('',6)),[(4,2),(3,3)]),((('Pi^2',4)),[(2,2)])]
amps=sys.argv[1:] or sorted(by)
random.seed(7)
pts=[(random.randrange(P),random.randrange(P)) for _ in range(2)]
res={}
for amp in amps:
    t0=time.time(); out={}
    for i,(t,u) in enumerate(pts):
        F=pure(amp,t,u)
        for key,comps in TESTS:
            if key not in F: continue
            for (a,b) in comps:
                D=delta(F[key],a,b)
                out.setdefault(f"{key[0] or '1'}*w{key[1]} d{a},{b}",[]).append((len(D),dict(collections.Counter(cls(k) for k in D))))
    res[amp]=out
    print(amp,round(time.time()-t0,1),json.dumps(out),flush=True)
