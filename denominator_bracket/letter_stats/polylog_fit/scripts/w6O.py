exec(open('core.py').read()); exec(open('sk2.py').read())
t0=time.time()
ALG={a for a in by if a[:6] in ('C1F0Z0','C7F0Z0','C9F0Z0')}
def d22_on_left(D):
    byy=collections.defaultdict(dict)
    for (x,y),c in D.items(): byy[y][x]=c
    out=collections.Counter()
    for y,F in byy.items():
        for (a,b),c in delta(F,2,2).items(): out[(a,b,y)]+=c
    return {k:v%P for k,v in out.items() if v%P}
def O(F): return d22_on_left(delta(F,4,2))
S=json.load(open('Li2_2_2.json'))
OA=O({tuple(w):c%P for w,c in S}); print("abstract O(Li222) terms:",len(OA),flush=True)
# also abstract O of depth-2 functions should vanish:
OD=O(dict((tuple(w),c%P) for w,c in C['symbols']['Li4_2'])); print("O(Li_{4,2}) abstract terms (expect 0):",len(OD))
OD=O(dict((tuple(w),c%P) for w,c in C['symbols']['Li3_3'])); print("O(Li_{3,3}) abstract terms (expect 0):",len(OD))
OD=O(dict((tuple(w),c%P) for w,c in C['symbols']['Li6'])); print("O(Li_6) abstract terms (expect 0):",len(OD))
T=json.load(open('triples.json'))
def substO(OA,vecs):
    sp_={a:[(LET[i],e) for i,e in enumerate(v) if e] for a,v in enumerate(vecs)}
    out=collections.Counter()
    for (x,y,z),c in OA.items():
        w=x+y+z
        for combo in itertools.product(*[sp_[a] for a in w]):
            m=c
            for _,e in combo: m*=e
            k=tuple(l for l,_ in combo); out[(k[0:2],k[2:4],k[4:6])]+=m
    return {k:v%P for k,v in out.items() if v%P}
K=3000; sk=Sketch(K,seed=5)
dv=data_vectors('',6,even_only=True,maxpts=10)
Dn=np.array([sk(O(v)) for a,v in dv if a not in ALG]); print("data depth-3 obstruction rank (non-alg):",rank_:=None,flush=True)
B=Basis(K); rows=[]; used=[]
for i,t in enumerate(T):
    v=sk(substO(OA,t['vecs']))
    if B.add(v): used.append(i)
    if i%600==599: print(i+1,"triples, span",len(B),f"{time.time()-t0:.0f}s",flush=True)
U=Basis(K); U.R=B.R.copy(); U.piv=list(B.piv)
print("Li222 depth-3 obstruction span:",len(B),"; data (non-alg) not covered:",U.add(Dn),flush=True)
json.dump(used,open('li222_used.json','w'))
