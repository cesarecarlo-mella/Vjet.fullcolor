# Is the rho-projected symbol of each eps^0 piece in span{rho(Li_n(x))}, x over the 36 admissible orbits?
import sys,os,random,collections,json,itertools
_argv=sys.argv[:]; sys.argv=['x','NONE']
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'classical.py')).read().split("TESTS=")[0]); sys.argv=_argv
A=json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'args.json'))); LET=A['letters']
def tensor(vecs):
    out=collections.Counter({():1})
    for v in vecs:
        new=collections.Counter()
        for w,c in out.items():
            for k,e in zip(LET,v):
                if e: new[w+(k,)]+=c*e
        out=new
    return out
def rhoS(S):
    out=collections.Counter()
    for w,c in S.items():
        for y,d in rho(w).items(): out[y]+=c*d
    return {k:v%P for k,v in out.items() if v%P}
def lin_vecs(n, linear_only):
    V=[]
    for a,b,c in A['triples']:
        for (X,Y,Z) in [(a,b,c),(b,a,c)]:          # Li_n(x) and Li_n(1-x)
            ex=[p-q for p,q in zip(X,Z)]; e1x=[p-q for p,q in zip(Y,Z)]
            if linear_only and any(ex[6:]+e1x[6:]): continue
            V.append(rhoS({w:-c_ for w,c_ in tensor([e1x]+[ex]*(n-1)).items()}))
    return V
basis={}
def add(v):
    v=dict(v)
    while v:
        piv=min(v)
        if piv not in basis: f=inv(v[piv]); basis[piv]={k:c*f%P for k,c in v.items()}; return 1
        r=basis[piv]; f=v[piv]
        for k,c in r.items():
            x_=(v.get(k,0)-f*c)%P
            if x_: v[k]=x_
            else: v.pop(k,None)
    return 0
const=sys.argv[1] if sys.argv[1]!='1' else ''; n=int(sys.argv[2]); linear_only=sys.argv[3:4]==['lin']
for v in lin_vecs(n,linear_only): add(v)
r0=len(basis); print(f"Li_{n} span ({'linear args' if linear_only else 'all 36 orbits'}): dim {r0}",flush=True)
rnd=random.Random(3); extra=collections.Counter(); tested=0
for amp in sorted(by):
    g=collections.defaultdict(collections.Counter)
    for _,c,wd,d,co,p in by[amp]:
        if c==const and len(wd)==n: g[(p,d)][wd]+=evc(co)
    if not g: continue
    for _ in range(6):
        env={'t':rnd.randrange(P),'u':rnd.randrange(P)}; F=collections.Counter()
        for (p,d),vec in g.items():
            r=eval(p.replace('^','**'),env)%P*inv(eval(d.replace('^','**'),env)%P)%P
            for x_,c in vec.items(): F[x_]+=r*c
        tested+=1
        if add(rhoS(F)): extra[amp]+=1
print(f"{sys.argv[1]}*w{n}: tested {tested} pure functions; rank grew by {len(basis)-r0} ->", "ALL IN SPAN" if len(basis)==r0 else dict(extra))
