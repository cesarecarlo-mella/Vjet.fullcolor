"""Symbol of G(a1..an; z) via Goncharov's I-recursion, entries factored over abstract letters."""
import sympy as sp, collections, itertools
from functools import lru_cache
X1,X2=sp.symbols('X1 X2')
ABS=[X1,X2,1-X1,1-X2,1-X1*X2]          # abstract letters 0..4
def factor_entry(e):
    e=sp.factor(sp.cancel(e))
    if e==0: return None
    num,den=sp.fraction(sp.together(e))
    vec=[0]*len(ABS)
    fl=[(f,m) for f,m in sp.factor_list(num)[1]]+[(f,-m) for f,m in sp.factor_list(den)[1]]
    for f,m in fl:
        f=sp.expand(f); hit=None
        for k,a in enumerate(ABS):
            if sp.expand(f-a)==0 or sp.expand(f+a)==0: hit=k
        if hit is None: raise ValueError(f"factor {f} not in alphabet")
        vec[hit]+=m
    return vec
@lru_cache(None)
def symI(args):
    """args = (a0, a1..an, a_{n+1}) as sympy exprs; returns Counter{tuple(letter idx): coef}"""
    n=len(args)-2
    if n==0: return collections.Counter({():1})
    out=collections.Counter()
    for i in range(1,n+1):
        sub=symI(args[:i]+args[i+1:])
        for nb,sg in ((args[i+1],1),(args[i-1],-1)):
            v=factor_entry(args[i]-nb)
            if v is None: continue
            for w,c in sub.items():
                for k,e in enumerate(v):
                    if e: out[w+(k,)]+=sg*c*e
    return collections.Counter({k:v for k,v in out.items() if v})
def symG(a,z): return symI(tuple([sp.Integer(0)]+list(reversed(a))+[z]))
def Li(ns,depth2=True):
    """Li_{n1,n2}(X1,X2) := G(0^{n2-1},1/X2,0^{n1-1},1/(X1 X2);1)   (or depth1: -G(0^{n-1},1/X1;1))"""
    if len(ns)==1: return collections.Counter({w:-c for w,c in symG([0]*(ns[0]-1)+[1/X1],sp.Integer(1)).items()})
    n1,n2=ns
    return symG([0]*(n2-1)+[1/X2]+[0]*(n1-1)+[1/(X1*X2)],sp.Integer(1))
