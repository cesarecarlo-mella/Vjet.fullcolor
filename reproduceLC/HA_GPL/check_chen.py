# numerical check of ii2g.chen(): path-ordered integral computed directly vs the G-form
import numpy as np, itertools, random, sys
import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import ii2g
from scipy.integrate import cumulative_trapezoid as ct
L={1:lambda X,Y:np.log(X),2:lambda X,Y:np.log(Y),3:lambda X,Y:np.log(1-X-Y),4:lambda X,Y:np.log(1-X),5:lambda X,Y:np.log(1-Y),6:lambda X,Y:np.log(X+Y)}
def dlog(k,X,Y,var,h=1e-7):   # analytic derivatives
    if var=='y': d={1:0*Y,2:1/Y,3:-1/(1-X-Y),4:0*Y,5:-1/(1-Y),6:1/(X+Y)}
    else: d={1:1/X,2:0*X,3:-1/(1-X-Y),4:-1/(1-X),5:0*X,6:1/(X+Y)}
    return d[k]
N=200001
def direct(word,x0,y0):
    s=np.linspace(0,y0,N); F=np.ones(N)
    Fs=[F]
    for a in word:
        with np.errstate(all='ignore'): integ=dlog(a,np.zeros(N),s,'y')*Fs[-1]
        integ[0]=0 if Fs[-1][0]==0 else integ[0]
        Fs.append(ct(integ,s,initial=0))
    vals=[f[-1] for f in Fs]       # prefix values at (0,y0)
    # x-leg: each prefix continues
    s=np.linspace(0,x0,N); cur=[np.full(N,v) for v in vals]
    out=[cur[0]]
    for j,a in enumerate(word):
        with np.errstate(all='ignore'): integ=dlog(a,s,np.full(N,y0),'x')*out[-1]
        if not np.isfinite(integ[0]): integ[0]=0
        out.append(vals[j+1]+ct(integ,s,initial=0))
    return out[-1][-1]
def Gnum(idx,z,y0):
    idx=[complex(eval(i.replace('y','(%r)'%y0))) for i in idx]
    s=np.linspace(0,z,N); G=np.ones(N,complex)
    for b in reversed(idx):
        with np.errstate(all='ignore'): integ=G/(s-b)
        if not np.isfinite(integ[0]): integ[0]=0
        G=ct(integ,s,initial=0)
    return G[-1]
def chen_num(word,x0,y0):
    tot=0
    for gx,gy in ii2g.chen(word):
        v=1
        for g,var,z in ((gx,'x',x0),(gy,'y',y0)):
            if g:
                idx=[a.strip() for a in g[2:-1].rsplit(',',1)[0].split(',')]
                v*=Gnum(idx,z,y0)
        tot+=v
    return tot
random.seed(3)
x0,y0=0.31,0.22
worst=0
for n in (1,2,3):
    for _ in range(6):
        w=(random.choice([3,4,5]),)+tuple(random.choice([2,3,4,5,6]) for _ in range(n-1))
        a=direct(w,x0,y0); b=chen_num(w,x0,y0)
        worst=max(worst,abs(a-b)); print(w,'direct %.9f  G-form %.9f%+.1ei'%(a,b.real,b.imag))
print('max |diff|',worst)
