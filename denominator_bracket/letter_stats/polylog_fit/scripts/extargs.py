# build extended argument list (rational 216 + algebraic), with exact values at random points for product checks
import json, random, collections, itertools
from fractions import Fraction as Fr
exec(open('algargs.py').read().split("mons=[]")[0].replace("for a,b,c in itertools.product","for a,b,c in []and itertools.product"))
C=json.load(open('pairs.json')); LETR=C['letters']
BASX=[1,2,3,4,5,6,11,12,13,14,15,16,7,8]
def clean(vec):  # vec over BAS (incl 9p) -> over BASX, with 9p = -(7+8)
    d=dict(zip(BAS,vec)); n9=d.pop('9p',0); d[7]=d.get(7,0)-n9; d[8]=d.get(8,0)-n9
    return [Fr(d.get(k,0)).limit_denominator(4) for k in BASX]
def value_from_vec(sign,vec,q):
    # only for rational args: vec over LETR
    s,t,u,r=q; v=sign%p
    for k,e in zip(LETR,vec):
        lv=EV[k](s,t,u)%p
        v=v*pow(lv,e%(p-1),p)%p
    return v
ARGS={}; NAME={}
import sympy as sp
S_,T_,U_,RHO=sp.symbols('s t u rho')
LN={1:S_,2:T_,3:U_,4:1-S_,5:1-T_,6:1-U_,11:S_-T_*(1-T_),12:T_-S_*(1-S_),13:(1-S_)**2-T_,14:(1-U_)**2-T_,15:(1-U_)**2-S_,16:(1-T_)**2-S_}
GN={'r':RHO,'a1':S_*U_+RHO,'a1b':S_*U_-RHO,'a2':S_*T_+RHO,'a2b':S_*T_-RHO,'a3':T_*U_+RHO,'a3b':T_*U_-RHO}
GN.update(LN)
def mexpr(m): return sp.Mul(*[GN[g]**e for g,e in m.items()])
for a,b in C['args']:
    val=tuple(value_from_vec(a[0],a[1],q) for q in pts)
    zv=[Fr(x) for x in a[1]]+[Fr(0),Fr(0)]; ov=[Fr(x) for x in b[1]]+[Fr(0),Fr(0)]
    ARGS[val]=(zv,ov); NAME[val]=str(sp.factor(a[0]*sp.Mul(*[LN[k]**e for k,e in zip(LETR,a[1])])))
nrat=len(ARGS)
# algebraic: recompute equations with values
alg_out=json.load(open('algargs.json'))['sols']
import re
def mval(m,q):
    v=1
    for g,e in m.items(): v=v*pow(gval(g,*q),e,p)%p
    return v
# redo search quickly to get monomials (reuse code path)
mons=[]
for deg in range(0,4):
    for combo in itertools.combinations_with_replacement(GEN,deg):
        mons.append(collections.Counter(combo))
Hv={}
for m in mons:
    v=tuple(mval(m,q) for q in pts); Hv.setdefault(v,(m,1)); Hv.setdefault(tuple((-x)%p for x in v),(m,-1))
def vec(m):
    v=collections.Counter()
    for g,e in m.items():
        for k,c in gvec(g).items(): v[k]+=c*e
    return [v[k] for k in BAS]
for m1,m2 in itertools.combinations(mons,2):
    if set(m1)&set(m2): continue
    if not any(not isinstance(g,int) for g in list(m1)+list(m2)): continue
    v1=tuple(mval(m1,q) for q in pts); v2=tuple(mval(m2,q) for q in pts)
    for sg in (1,-1):
        sv=tuple((a+sg*b)%p for a,b in zip(v1,v2))
        if sv in Hv:
            m3,s3=Hv[sv]
            if set(m3)&set(m1) or set(m3)&set(m2): continue
            # z = m1/(s3 m3), 1-z = -sg m2/(s3 m3)
            z=tuple(a*pow(s3*c%p,p-2,p)%p for a,c in zip(v1,[mval(m3,q) for q in pts]))
            X=clean([a-b for a,b in zip(vec(m1),vec(m3))]); Y=clean([a-b for a,b in zip(vec(m2),vec(m3))])
            one=lambda zz: tuple((1-x)%p for x in zz)
            invv=lambda zz: tuple(pow(x,p-2,p) for x in zz)
            members=[(z,X,Y)]
            members.append((one(z),Y,X))
            members.append((invv(z),[-x for x in X],[y-x for x,y in zip(X,Y)]))
            members.append((invv(one(z)),[-y for y in Y],[x-y for x,y in zip(X,Y)]))
            zz=tuple(a*pow((a-1)%p,p-2,p)%p for a in z); members.append((zz,[x-y for x,y in zip(X,Y)],[-y for y in Y]))
            zz=tuple((a-1)*pow(a,p-2,p)%p for a in z); members.append((zz,[y-x for x,y in zip(X,Y)],[-x for x in X]))
            z0=mexpr(m1)/(s3*mexpr(m3))
            exprs=[z0,1-z0,1/z0,1/(1-z0),z0/(z0-1),(z0-1)/z0]
            for (val,zv,ov),ex in zip(members,exprs):
                if val not in ARGS: ARGS[val]=(zv,ov); NAME[val]=str(sp.simplify(ex))
print("arguments: rational",nrat,"algebraic",len(ARGS)-nrat)
keys=list(ARGS)
pairs=[]
for a,b in itertools.product(keys,repeat=2):
    if a==b: continue
    prod=tuple(x*y%p for x,y in zip(a,b))
    if prod in ARGS:
        pairs.append((a,b,prod))
print("pairs with product admissible:",len(pairs),"(involving algebraic:",sum(1 for a,b,c in pairs if any(ARGS[x][0][12] or ARGS[x][0][13] or ARGS[x][1][12] or ARGS[x][1][13] for x in (a,b))),")")
ser=lambda v:[str(x) for x in v]
json.dump({"basis":BASX,"args":[[ser(ARGS[k][0]),ser(ARGS[k][1])] for k in keys],
           "pairs":[[ser(ARGS[a][0]),ser(ARGS[b][0]),ser(ARGS[a][1]),ser(ARGS[b][1]),ser(ARGS[c][1])] for a,b,c in pairs],
           "names":[NAME[k] for k in keys],"alg_flag":[bool(ARGS[k][0][12] or ARGS[k][0][13] or ARGS[k][1][12] or ARGS[k][1][13]) for k in keys]},open('ext_args.json','w'))
