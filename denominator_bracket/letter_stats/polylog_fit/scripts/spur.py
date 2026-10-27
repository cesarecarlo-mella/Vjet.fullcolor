import sympy as sp, json, itertools, collections, time, pickle
exec(open('letters.py').read())
C=json.load(open('pairs.json')); LETR=C['letters']
args=[(a[0],tuple(a[1])) for a,b in C['args']]
adm=set(args)
mul=lambda a,b:(a[0]*b[0],tuple(x+y for x,y in zip(a[1],b[1])))
W=collections.defaultdict(list)
for i,a in enumerate(args):
    for j,b in enumerate(args):
        if i==j: continue
        w=mul(a,b)
        if w not in adm: W[w].append((i,j))
print("non-admissible products:",len(W),"pairs:",sum(len(v) for v in W.values()),flush=True)
FULL=[1,2,3,4,5,6,9,10,11,12,13,14,15,16,17,18,19,20]
t0=time.time(); res={}; newf={}
for n,(w,prs) in enumerate(W.items()):
    e=w[0]*sp.Mul(*[L[k]**ex for k,ex in zip(LETR,w[1])])
    num,den=sp.fraction(sp.together(1-e))
    fl=sp.factor_list(sp.expand(num))[1]+[(f,-mm) for f,mm in sp.factor_list(sp.expand(den))[1]]
    vec=collections.Counter()
    for f,mm in fl:
        f=sp.expand(f); hit=None
        for k in FULL:
            if sp.expand(f-L[k])==0 or sp.expand(f+L[k])==0: hit=k;break
        if hit is None:
            key=str(f) if str(f) in newf or str(sp.expand(-f)) not in newf else str(sp.expand(-f))
            newf.setdefault(key,len(newf)); hit=('N',newf[key])
        vec[hit]+=mm
    res[w]=dict(vec)
    if n%2000==0: print(n,len(newf),f"{time.time()-t0:.0f}s",flush=True)
pickle.dump((args,dict(W),res,newf),open('spur.pkl','wb'))
print("done; new irreducible factors:",len(newf))
