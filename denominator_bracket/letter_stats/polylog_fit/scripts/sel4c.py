exec(open('sel4.py').read().split("ch,unc=greedy")[0])
ch,unc=greedy(Cq,Dq,order)
best=None
import random as rr
for trial in range(30):
    rest=[i for i in order if i not in ch]
    if trial: rr.Random(trial).shuffle(rest)
    used,_=support_basis(Cq,Dq,ch+rest)
    if best is None or len(used)<len(best): best=used
print("best support:",len(best))
json.dump([E['pairs'][i][:2] for i in best],open('sel_w4.json','w'))
