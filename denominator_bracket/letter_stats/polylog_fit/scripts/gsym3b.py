exec(open('gsym3.py').read().split("S=symG(")[0])
import json
out={}
for ns in [(2,2,1),(2,1,2),(1,2,2),(3,1,1),(1,3,1),(1,1,3)]:
    n1,n2,n3=ns
    S=symG([0]*(n3-1)+[1/X3]+[0]*(n2-1)+[1/(X2*X3)]+[0]*(n1-1)+[1/(X1*X2*X3)],sp.Integer(1))
    out["Li"+"_".join(map(str,ns))]=[[list(w),int(c)] for w,c in S.items()]; print(ns,len(S))
json.dump(out,open('depth3_w5.json','w'))
