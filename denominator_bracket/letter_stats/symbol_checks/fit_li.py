"""Is the weight-W part of the eps^0 remainders (even letters, mod products) in
   span{ rho Li_W(x) } + span{ rho Li_{a,b}(x1,x2), a+b=W } with admissible arguments?
Usage: python3 fit_li.py CONST W FUNCS        e.g.  fit_li.py 1 4 Li2_2
       FUNCS: comma list of depth-2 functions from li_candidates.json (Li2_2,Li3_1,...), or 'none'
"""
import sys,os,random,collections,json,itertools,time
HERE=os.path.dirname(os.path.abspath(__file__))
_argv=sys.argv[:]; sys.argv=['x','NONE']
exec(open(os.path.join(HERE,'classical.py')).read().split("TESTS=")[0]); sys.argv=_argv
const='' if sys.argv[1]=='1' else sys.argv[1]; W=int(sys.argv[2]); FUN=[] if sys.argv[3]=='none' else sys.argv[3].split(',')
C=json.load(open(os.path.join(HERE,'li_candidates.json'))); LET=C['letters']
def rhoA(S):   # rho on abstract symbol (exact ints scaled mod P)
    out=collections.Counter()
    for w,c in S:
        for y,d in rho(tuple(w)).items(): out[y]+=c*d
    return {k:v%P for k,v in out.items() if v%P}
def subst(rA, vecs):
    out=collections.Counter()
    for w,c in rA.items():
        parts=[[(LET[i],e) for i,e in enumerate(vecs[a]) if e] for a in w]
        for combo in itertools.product(*parts):
            m=c
            for _,e in combo: m*=e
            out[tuple(k for k,_ in combo)]+=m
    return {k:v%P for k,v in out.items() if v%P}
class Span:
    def __init__(s): s.b={}
    def add(s,v):
        v=dict(v); b=s.b
        while v:
            piv=min(v)
            if piv not in b: f=inv(v[piv]); b[piv]={k:c*f%P for k,c in v.items()}; return 1
            r=b[piv]; f=v[piv]
            for k,c in r.items():
                x_=(v.get(k,0)-f*c)%P
                if x_: v[k]=x_
                else: v.pop(k,None)
        return 0
    def copy(s): n=Span(); n.b={k:dict(v) for k,v in s.b.items()}; return n
t0=time.time()
# classical Li_W(x) for all 216 arguments: abstract letters 0 -> x, 2 -> 1-x
rLi=rhoA(C['symbols'][f'Li{W}'])
S=Span()
for a,b in C['args']:
    vecs={0:a[1],2:b[1]}
    S.add(subst(rLi,vecs))
r_cl=len(S.b); print(f"Li_{W} span: {r_cl}",flush=True)
# data
rnd=random.Random(3); data=[]
for amp in sorted(by):
    g=collections.defaultdict(collections.Counter)
    for _,c,wd,d,co,p in by[amp]:
        if c==const and len(wd)==W and not {7,8}&set(wd): g[(p,d)][wd]+=evc(co)
    if not g: continue
    for _ in range(6):
        env={'t':rnd.randrange(P),'u':rnd.randrange(P)}; F=collections.Counter()
        for (p,d),vec in g.items():
            r=eval(p.replace('^','**'),env)%P*inv(eval(d.replace('^','**'),env)%P)%P
            for x_,c in vec.items(): F[x_]+=r*c
        out=collections.Counter()
        for w,c in F.items():
            for y,dd in rho(w).items(): out[y]+=c*dd
        data.append({k:v%P for k,v in out.items() if v%P})
T=S.copy(); grow=sum(T.add(v) for v in data)
print(f"{sys.argv[1]}*w{W} (even letters): {len(data)} pure functions; outside Li_{W} span: {grow} dims",flush=True)
if FUN:
    D=S.copy(); cand=[]
    for fn in FUN:
        rA=rhoA(C['symbols'][fn])
        for p in C['pairs']:
            vecs={0:p['x1'][1],1:p['x2'][1],2:p['1-x1'][1],3:p['1-x2'][1],4:p['1-x1x2'][1]}
            cand.append((fn,p,subst(rA,vecs)))
    for _,_,v in cand: D.add(v)
    print(f"Li_{W} + {FUN} on {len(C['pairs'])} pairs: span {len(D.b)} (+{len(D.b)-r_cl})",flush=True)
    E=D.copy(); g2=sum(E.add(v) for v in data)
    print(">>> data inside span:", g2==0, f"(rank growth {g2})",flush=True)
    # greedy minimal choice: depth-2 functions needed modulo Li_W, restricted to span(Li_W + data)
    SD=S.copy(); [SD.add(v) for v in data]
    chosen=[]; cur=S.copy()
    for fn,p,v in cand:
        if SD.copy().add(v)==0 and cur.copy().add(v)==1: cur.add(v); chosen.append((fn,p))
    full=cur.copy(); g3=sum(full.add(v) for v in data)
    print(f"greedy: {len(chosen)} single depth-2 functions lying in span(Li_{W}+data); with them data covered: {g3==0} (missing {g3})",flush=True)
    json.dump([[fn,p['x1'],p['x2']] for fn,p in chosen],open(os.path.join(HERE,f'chosen_{sys.argv[1]}_w{W}.json'),'w'))
print(f"time {time.time()-t0:.0f}s")
