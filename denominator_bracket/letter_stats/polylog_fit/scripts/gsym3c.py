exec(open('gsym3.py').read().split("S=symG(")[0])
import json, itertools
out={}
for ns in [c for c in itertools.product(range(1,5),repeat=3) if sum(c)==6]:
    n1,n2,n3=ns
    S=symG([0]*(n3-1)+[1/X3]+[0]*(n2-1)+[1/(X2*X3)]+[0]*(n1-1)+[1/(X1*X2*X3)],sp.Integer(1))
    out["Li"+"_".join(map(str,ns))]=[[list(w),int(c)] for w,c in S.items()]; print(ns,len(S),flush=True)
json.dump(out,open('depth3_w6.json','w'))
