#!/usr/bin/env python3
"""accuracy study of the Zgqq evaluator -- an external job: it uses only the public interface of the package
(zgqq_eval.evaluate) and changes nothing in it.

N points, uniform in the decay triangle (fixed seed); this job evaluates chunk i of n of them:
  double precision with its error estimate (two evaluations: independent start values, different complex paths);
  with --quad also full quad precision with its own two-path error estimate (the reference).
Writes <out>/chunk_<i>.npz.

usage: accuracy_job.py --package DIR --npts 1000 --chunk 0 --nchunks 20 [--quad] [--seed 1] --out DIR"""
import argparse, os, sys, time
import numpy as np


def triangle_points(n, seed):
    rng = np.random.default_rng(seed)
    u = rng.random((n, 2))
    flip = u.sum(1) > 1
    u[flip] = 1 - u[flip]
    return u


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--package', required=True, help='directory of the built package (zgqq_eval.py)')
    ap.add_argument('--npts', type=int, default=1000)
    ap.add_argument('--chunk', type=int, default=0)
    ap.add_argument('--nchunks', type=int, default=1)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--quad', action='store_true')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    sys.path.insert(0, os.path.abspath(a.package))
    import zgqq_eval as Z

    P = triangle_points(a.npts, a.seed)
    idx = np.array_split(np.arange(a.npts), a.nchunks)[a.chunk]
    pts = [tuple(P[i]) for i in idx]
    out = dict(index=idx, points=P[idx], coefficients=np.array(Z.COEFFICIENTS), layers=np.array(Z.LAYERS),
               version=Z.__version__)
    t0 = time.time()
    rd = Z.evaluate(pts, error=True)
    out.update(vd=rd.values, ed=rd.errors, td=(time.time() - t0) / max(len(pts), 1))
    print('chunk %d: %d points, double %.2f s/pt' % (a.chunk, len(pts), out['td']), flush=True)
    if a.quad:
        t0 = time.time()
        rq = Z.evaluate(pts, error=True, precision='quad')
        out.update(vq=rq.values, vq_lo=rq.values_lo, eq=rq.errors, tq=(time.time() - t0) / max(len(pts), 1))
        print('chunk %d: quad %.1f s/pt' % (a.chunk, out['tq']), flush=True)
    os.makedirs(a.out, exist_ok=True)
    np.savez(os.path.join(a.out, 'chunk_%d.npz' % a.chunk), **out)


if __name__ == '__main__':
    main()
