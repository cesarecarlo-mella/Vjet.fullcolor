exec(open('fast.py').read()); exec(open('select.py').read())
from fractions import Fraction as Fr
t0=time.time(); K=1500; sk=R1(K,5,seed=1)
Bcl=Basis(K); W,Cf=abstract_rho(C['symbols']['Li5'])
Bcl.add(np.array([sk.abstract(W,Cf,{0:[m(x) for x in z],2:[m(x) for x in o]}) for z,o in E['args']]))
W,Cf=abstract_rho(C['symbols']['Li3_2'])
c2=[sk.abstract(W,Cf,{0:[m(x) for x in x1],1:[m(x) for x in x2],2:[m(x) for x in o1],3:[m(x) for x in o2],4:[m(x) for x in o12]}) for x1,x2,o1,o2,o12 in E['pairs']]
T=json.load(open('triples.json')); D3=json.load(open('depth3_w5.json')); W,Cf=abstract_rho(D3['Li2_2_1'])
ext=lambda v:[m(x) for x in v]+[0,0]
c3=[sk.abstract(W,Cf,{a:ext(v) for a,v in enumerate(t['vecs'])}) for t in T]
data=[sk.data(v) for a,v in data_vectors('',5,even_only=False)]
print("built",f"{time.time()-t0:.0f}s",flush=True)
cands=c2+c3
Cq,Dq,n=quotient_coords(Bcl,cands,data,K); print("quotient dim",n,"data dim mod Li5",Basis(n).add(Dq),flush=True)
def cost(i):
    if i<len(c2):
        pr=E['pairs'][i]; alg=any(Fr(v)!=0 for vec in pr[:2] for v in vec[12:14])
        return (0,alg,sum(abs(Fr(v)) for vec in pr[:2] for v in vec))
    return (1,0,i)
order=sorted(range(len(cands)),key=cost)
best=None
for trial in range(10):
    o=order if trial==0 else sorted(order,key=lambda i:(cost(i)[0],random.Random(trial*7919+i).random()))
    used,_=support_basis(Cq,Dq,o)
    if best is None or len(used)<len(best): best=used
S=Basis(n); S.add(Cq[best]); T_=Basis(n); T_.R=S.R.copy(); T_.piv=list(S.piv); print("check uncovered:",T_.add(Dq))
p2=[i for i in best if i<len(c2)]; p3=[i-len(c2) for i in best if i>=len(c2)]
print("support: Li32",len(p2),"Li221",len(p3))
json.dump({"Li3_2":[E['pairs'][i][:2] for i in p2],"Li3_2_idx":p2,"Li2_2_1":p3},open('sel_w5.json','w'))
