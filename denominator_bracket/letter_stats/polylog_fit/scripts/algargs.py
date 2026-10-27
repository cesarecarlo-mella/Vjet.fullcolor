# S-unit search including the square-root sector: rho = sqrt(-s t u) (= i*Sqrt[x y (1-x-y)])
import random, itertools, collections, json
p=(1<<61)-1            # p = 3 mod 4
def sqrt_mod(a): r=pow(a,(p+1)//4,p); return r if r*r%p==a%p else None
EV={1:lambda s,t,u:s,2:lambda s,t,u:t,3:lambda s,t,u:u,4:lambda s,t,u:1-s,5:lambda s,t,u:1-t,6:lambda s,t,u:1-u,
    11:lambda s,t,u:s-t*(1-t),12:lambda s,t,u:t-s*(1-s),13:lambda s,t,u:(1-s)**2-t,14:lambda s,t,u:(1-u)**2-t,
    15:lambda s,t,u:(1-u)**2-s,16:lambda s,t,u:(1-t)**2-s}
pts=[]
while len(pts)<3:
    s,t=random.randrange(2,p),random.randrange(2,p); u=(1-s-t)%p
    r=sqrt_mod((-s*t*u)%p)
    if r: pts.append((s,t,u,r))
# generators: even letters, rho, alpha1=su+rho, alpha1b=su-rho, alpha2=st+rho, alpha2b, alpha3=tu+rho, alpha3b
GEN=list(EV)+['r','a1','a1b','a2','a2b','a3','a3b']
def gval(g,s,t,u,r):
    if g in EV: return EV[g](s,t,u)%p
    return {'r':r,'a1':s*u+r,'a1b':s*u-r,'a2':s*t+r,'a2b':s*t-r,'a3':t*u+r,'a3b':t*u-r}[g]%p
V={g:[gval(g,*q) for q in pts] for g in GEN}
# vector representation over basis [1..6,11..16, 7, 8, 9'] (9' = (tu-rho)/(tu+rho)), half-integers allowed
BAS=[1,2,3,4,5,6,11,12,13,14,15,16,7,8,'9p']
def gvec(g):
    v=collections.Counter()
    if g in EV: v[g]=1
    elif g=='r': v[1]=v[2]=v[3]=0.5
    else:
        N={'1':[1,3,4,6],'2':[1,2,4,5],'3':[2,3,5,6]}[g[1]]; odd={'1':7,'2':8,'3':'9p'}[g[1]]
        for k in N: v[k]+=0.5
        v[odd]+= -0.5 if not g.endswith('b') else 0.5
    return v
# check multiplicative dependence of the three odd phases: l7^a l8^b l9'^c = +-1 ?
def ph(i,q):
    s,t,u,r=q; a={1:s*u,2:s*t,3:t*u}[i]; return (a-r)*pow(a+r,p-2,p)%p
for a,b,c in itertools.product(range(-2,3),repeat=3):
    if (a,b,c)!=(0,0,0) and all(pow(ph(1,q),a%(p-1),p)*pow(ph(2,q),b%(p-1),p)*pow(ph(3,q),c%(p-1),p)%p in (1,p-1) for q in pts):
        print("odd-phase relation: l7^%d l8^%d l9'^%d = +-1"%(a,b,c))
mons=[]
for deg in range(0,4):
    for combo in itertools.combinations_with_replacement(GEN,deg):
        val=[1]*len(pts)
        for g in combo: val=[x*y%p for x,y in zip(val,V[g])]
        mons.append((collections.Counter(combo),tuple(val)))
H={}
for m,v in mons:
    H.setdefault(v,(m,1)); H.setdefault(tuple((-x)%p for x in v),(m,-1))
sols=set(); out=[]
alg=lambda m: any(not isinstance(g,int) for g in m)
for (m1,v1),(m2,v2) in itertools.combinations(mons,2):
    if set(m1)&set(m2): continue
    for sg in (1,-1):
        sv=tuple((a+sg*b)%p for a,b in zip(v1,v2))
        if sv in H:
            m3,s3=H[sv]
            if set(m3)&set(m1) or set(m3)&set(m2): continue
            if not (alg(m1) or alg(m2) or alg(m3)): continue
            key=frozenset([tuple(sorted(m1.items(),key=str)),tuple(sorted(m2.items(),key=str)),tuple(sorted(m3.items(),key=str))])
            if key in sols: continue
            sols.add(key)
            def vec(m):
                v=collections.Counter()
                for g,e in m.items():
                    for k,c in gvec(g).items(): v[k]+=c*e
                return v
            X={k:vec(m1)[k]-vec(m3)[k] for k in BAS}; Y={k:vec(m2)[k]-vec(m3)[k] for k in BAS}
            out.append({"z":[X[k] for k in BAS],"1-z":[Y[k] for k in BAS],
                        "eq":f"{dict(m1)} {'+' if sg>0 else '-'} {dict(m2)} = {'+' if s3>0 else '-'}{dict(m3)}"})
print("algebraic S-unit equations (orbits):",len(out))
for o in out[:40]: print("  ",o["eq"])
json.dump({"basis":[str(b) for b in BAS],"sols":out},open('algargs.json','w'))
