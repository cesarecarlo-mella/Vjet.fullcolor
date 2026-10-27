"""Move the single quadratic letter of each eps^0 word to the LAST integration by shuffle:
   II[u,q,v1..vm] = sum_k (-1)^k II[(u sh v_k..v_1), q] * II[v_{k+1}..v_m]
Verifies the identity exactly on every quadratic word and prints what the rewriting produces."""
import pickle,collections,os,itertools,random
from functools import lru_cache
import sys, importlib.util
D=sys.argv[1] if len(sys.argv)>1 else os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import re, glob
rows=[]
for f in sorted(glob.glob(os.path.join(D,'*.wl'))):
    L=open(f).read().splitlines(); i=next(k for k,l in enumerate(L) if l.startswith('+eps^(0)*'))
    for l in L[i+1:]:
        m=re.search(r'^\+\(.*\)\*\((.*?)II\[([0-9,]*)\]\)$',l)
        if m: rows.append((os.path.basename(f)[:-3], m.group(1).rstrip('*'), tuple(int(x) for x in m.group(2).split(',')) if m.group(2) else ()))
Q=set(range(9,17))
@lru_cache(None)
def sh(a,b):
    if not a: return {b:1}
    if not b: return {a:1}
    out=collections.Counter()
    for w,c in sh(a[:-1],b).items(): out[w+(a[-1],)]+=c
    for w,c in sh(a,b[:-1]).items(): out[w+(b[-1],)]+=c
    return dict(out)
def decompose(w):
    """w = u q v  ->  list of (coef, A, B) with A ending in q, B linear; w = sum coef * A sh B"""
    i=[k for k,l in enumerate(w) if l in Q][0]; u,q,v=w[:i],w[i],w[i+1:]
    res=collections.Counter()
    for k in range(len(v)+1):
        rev=tuple(reversed(v[:k]))
        for x,c in sh(u,rev).items(): res[(x+(q,),v[k:])]+=(-1)**k*c
    return {k:c for k,c in res.items() if c}
def expand(dec):
    out=collections.Counter()
    for (A,B),c in dec.items():
        for x,m in sh(A,B).items(): out[x]+=c*m
    return {k:v for k,v in out.items() if v}
# 1) identity check on random words and on all quadratic words of the data
random.seed(1)
for _ in range(2000):
    n=random.randint(1,7); w=[random.choice([1,2,3,4,5,6]) for _ in range(n)]; w.insert(random.randint(0,n),11); w=tuple(w)
    assert expand(decompose(w))=={w:1},w
qw=sorted({r[2] for r in rows if Q&set(r[2])})
for w in qw: assert expand(decompose(w))=={w:1},w
print("identity verified on 2000 random words and on all",len(qw),"quadratic words of the data")
# 2) what the rewriting produces
A=collections.Counter(); B=collections.Counter(); pairs=set(); nterms=collections.Counter()
for w in qw:
    d=decompose(w); nterms[len(d)]+=1
    for (a,b) in d: A[a]+=1; B[b]+=1; pairs.add((a,b))
print("terms per decomposed word:",dict(sorted(nterms.items())))
print("distinct q-final words II[..,q]:",len(A),"by length",dict(sorted(collections.Counter(len(a) for a in A).items())))
print("distinct linear cofactors II[v]:",len(B),"by length",dict(sorted(collections.Counter(len(b) for b in B).items())))
print("distinct (q-final, cofactor) pairs:",len(pairs))
print("q-final words per quadratic letter:",dict(sorted(collections.Counter(a[-1] for a in A).items())))
# effective depth: prefix length before q in q-final words, with the constant
P=collections.Counter()
for a,c,w,*_ in rows:
    if Q&set(w):
        for (x,y) in decompose(w): P[(c,len(x)-1,len(y))]+=1
print("(const, prefix length before q, cofactor length):",dict(sorted(P.items())))
