import json, collections
C=json.load(open('pairs.json'))
A={(a[0],tuple(a[1])):(b[0],tuple(b[1])) for a,b in C['args']}
def mul(a,b): return (a[0]*b[0],tuple(x+y for x,y in zip(a[1],b[1])))
keys=list(A)
pairs=collections.defaultdict(list)
for p in C['pairs']: pairs[(p['x1'][0],tuple(p['x1'][1]))].append((p['x2'][0],tuple(p['x2'][1])))
T=[]
for x1 in keys:
    for x2 in pairs[x1]:
        x12=mul(x1,x2)
        for x3 in pairs[x2]:
            x123=mul(x12,x3)
            if x123 in A and mul(x2,x3) in A:
                T.append((x1,x2,x3))
print("admissible triples:",len(T))
out=[]
for x1,x2,x3 in T:
    x12,x23,x123=mul(x1,x2),mul(x2,x3),mul(mul(x1,x2),x3)
    out.append({"vecs":[list(v[1]) for v in [x1,x2,x3,A[x1],A[x2],A[x3],A[x12],A[x23],A[x123]]],"x":[[v[0],list(v[1])] for v in (x1,x2,x3)]})
json.dump(out,open('triples.json','w'))
