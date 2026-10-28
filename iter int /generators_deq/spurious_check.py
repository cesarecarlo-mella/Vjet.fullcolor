"""end points ON the curves of the spurious letters 11-16 (where the direct integration fails and the rational functions
of letters 15, 16 have poles): the remainders from the stencil average across the curve (gen_eval.plan).
(i) stability: step H vs H/2;  (ii) correctness: 1e-4 off the curve, where the direct evaluation works, the direct
remainders vs the average forced at the same point."""
import numpy as np
import gen_eval
from remainder_eval import evaluate_remainders, remainders_from

PTS = [(0.21, 0.7), (0.25, 0.5), (0.16, 0.2), (0.3, 0.21), (0.5, 0.25), (0.6, 0.16), (0.24, 0.36), (0.36, 0.24),
       (0.16, 0.6)]


def rel(a, b):
    return (np.abs(a - b) / np.maximum(np.abs(b), 1)).max(axis=-1)


R, info = evaluate_remainders(PTS, 1e-12)
H = gen_eval.H
gen_eval.H = H / 2
R2, _ = evaluate_remainders(PTS, 1e-12)
gen_eval.H = H
Q = []
for P in PTS:
    n = gen_eval.spurious_near(P)[1]
    Q.append((P[0] + 1e-4 * n[0], P[1] + 1e-4 * n[1]))
T, _ = gen_eval.evaluate_direct(Q, 1e-12)
direct = remainders_from(Q, T)
near = gen_eval.NEAR
gen_eval.NEAR = {k: 1e-2 for k in near}               # force the average at Q (1e-4 off the curve)
avg, _ = evaluate_remainders(Q, 1e-12)
gen_eval.NEAR = near
for P, inf, a, b in zip(PTS, info, rel(R, R2), rel(direct, avg)):
    print('%-12s letter %d: H vs H/2 %.1e | 1e-4 off the curve, direct vs average %.1e' % (P, inf['averaged']['letter'], a, b))
