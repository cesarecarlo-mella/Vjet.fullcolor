import numpy as np, scipy.sparse as sps, collections, itertools, json, time, random
from fractions import Fraction as Fr
exec(open('core.py').read()); exec(open('sk2.py').read())
E=json.load(open('ext_args.json')); BX=E['basis']; LI={l:i for i,l in enumerate(BX)}
def m(x): f=Fr(x); return f.numerator%P*inv(f.denominator%P)%P
class R1:
    """rank-one sketch phi_j(w)=prod_i R[j,i,w_i] (mod p)"""
    def __init__(s,K,W,seed=0):
        rng=np.random.default_rng(seed); s.K,s.W=K,W
        s.R=rng.integers(1,p,size=(W,len(BX),K),dtype=np.int64)   # [pos, letter, j]
        h=W//2; s.h=h
        s.pre={}; s.suf={}
        for w in itertools.product(range(len(BX)),repeat=h):
            v=np.ones(K,dtype=np.int64)
            for i,l in enumerate(w): v=v*s.R[i,l]%p
            s.pre[w]=v
        for w in itertools.product(range(len(BX)),repeat=W-h):
            v=np.ones(K,dtype=np.int64)
            for i,l in enumerate(w): v=v*s.R[h+i,l]%p
            s.suf[w]=v
        s.prekeys=list(s.pre); s.sufkeys=list(s.suf)
        s.preidx={k:i for i,k in enumerate(s.prekeys)}; s.sufidx={k:i for i,k in enumerate(s.sufkeys)}
        s.A=np.array([s.pre[k] for k in s.prekeys]); s.B=np.array([s.suf[k] for k in s.sufkeys])
    def data(s,v):
        rows=[];cols=[];vals=[]
        for w,c in v.items():
            ww=tuple(LI[l] for l in w); rows.append(s.preidx[ww[:s.h]]); cols.append(s.sufidx[ww[s.h:]]); vals.append(c%p)
        M=sps.csr_matrix((np.array(vals,dtype=np.int64),(rows,cols)),shape=(len(s.prekeys),len(s.sufkeys)))
        T=(M@s.B)%p                      # (npre x K)
        return np.asarray((T*s.A%p).sum(axis=0)%p).ravel()
    def abstract(s,rA_words,rA_coef,vecs):
        """rA_words: (n,W) abstract letter idx; vecs: dict a -> list over BX (mod p ints)"""
        G={}
        for a,vec in vecs.items():
            va=np.array(vec,dtype=np.int64)%p
            G[a]=np.stack([(va[:,None]*s.R[i]%p).sum(axis=0)%p for i in range(s.W)])  # (W,K)
        Ga=np.stack([G[a] for a in sorted(G)]) ; order={a:k for k,a in enumerate(sorted(G))}
        idx=np.vectorize(order.get)(rA_words)
        prod=Ga[idx[:,0],0]
        for i in range(1,s.W): prod=prod*Ga[idx[:,i],i]%p
        return (prod*rA_coef[:,None]%p).sum(axis=0)%p
def abstract_rho(S):
    r=rhoA(S); ws=list(r); return np.array(ws,dtype=np.int64),np.array([r[w] for w in ws],dtype=np.int64)
