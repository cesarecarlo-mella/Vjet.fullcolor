"""evaluate the 60 generators H[1..60] at points of the Euclidean triangle with tridents_bs + generators_bs.bin.

    from gen_eval import evaluate
    H, info = evaluate([(0.3, 0.2), (0.1, 0.7)])     # H.shape == (2, 60), H[p, i] = H[i+1] at point p
    python3 gen_eval.py x y                           # prints the 60 values
"""
import os, pickle, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WT = pickle.load(open(os.path.join(HERE, 'out/system.pkl'), 'rb'))['wt']


def evaluate(points, tol=1e-10, binary=os.path.join(HERE, 'cpp/tridents_bs'), data=os.path.join(HERE, 'out/generators_bs.bin')):
    inp = '\n'.join('%.17g %.17g' % tuple(p) for p in points)
    r = subprocess.run([binary, data, str(tol), str(tol)], input=inp, capture_output=True, text=True, check=True)
    H, info, cur = [], [], None
    for l in r.stdout.split('\n'):
        if l.startswith('#'):
            h = l.split(); info.append(dict(x=float(h[1]), y=float(h[2]), steps=int(h[3]), evals=int(h[4]), time=float(h[5])))
            cur = []; H.append(cur)
        elif l:
            a, b = l.split(); cur.append(complex(float(a), float(b)))
    NT = len(WT)
    return np.array([[np.array(v).reshape(7, NT)[w, i] for i, w in enumerate(WT)] for v in H]), info


if __name__ == '__main__':
    H, info = evaluate([(float(sys.argv[1]), float(sys.argv[2]))])
    for i, v in enumerate(H[0]):
        print('H[%d] (weight %d) = %.15e %+.15e i' % (i + 1, WT[i], v.real, v.imag))
    print(info[0])
