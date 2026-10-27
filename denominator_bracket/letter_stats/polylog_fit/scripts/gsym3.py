import sympy as sp, collections
from functools import lru_cache
X1,X2,X3=sp.symbols('X1 X2 X3')
ABS=[X1,X2,X3,1-X1,1-X2,1-X3,1-X1*X2,1-X2*X3,1-X1*X2*X3]
def factor_entry(e):
    e=sp.cancel(e)
    if e==0: return None
    num,den=sp.fraction(sp.together(e)); vec=[0]*len(ABS)
    fl=[(f,m) for f,m in sp.factor_list(num)[1]]+[(f,-m) for f,m in sp.factor_list(den)[1]]
    for f,m in fl:
        hit=[k for k,a in enumerate(ABS) if sp.expand(f-a)==0 or sp.expand(f+a)==0]
        assert hit, f
        vec[hit[0]]+=m
    return vec
@lru_cache(None)
def symI(args):
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
S=symG([0,1/X3,0,1/(X2*X3),0,1/(X1*X2*X3)],sp.Integer(1))
import json; json.dump([[list(w),int(c)] for w,c in S.items()],open('Li2_2_2.json','w'))
print("Li_{2,2,2} symbol terms:",len(S))
