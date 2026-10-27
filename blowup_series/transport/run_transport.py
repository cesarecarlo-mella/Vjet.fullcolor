#!/usr/bin/env python3
"""
run_transport.py -- boundary constants of the continued charts for all families and regions (cluster driver).

Per family F and region A > 1 (others b, c):
    edge     edge_transport.wl     bdy_R<A>_th  ->  bdy_R<A>_<b>v, bdy_R<A>_<c>v      (s = 0 charts)
    div_b    divisor_transport.wl  bdy_R<A>_<b>v ->  bdy_R<A>_v<b>                      (s = 1 chart)
    div_c    divisor_transport.wl  bdy_R<A>_<c>v ->  bdy_R<A>_v<c>
The threshold constants bdy_R<A>_th_w<w>.m must exist (threshold.py consts).  One Mathematica kernel per job,
-n jobs at a time; a divisor job starts when its edge job is done.  Resumable: a job whose last-weight output
exists is skipped.  Logs: transport/logs/<F>_<A>_<step>.log.

    python3 run_transport.py --dry-run
    nohup python3 run_transport.py -n 16 > logs/run_transport.log 2>&1 &
    python3 run_transport.py -n 4 --families H,I3 --regions z
    python3 run_transport.py -n 2 --steps edge          # only the (cheap) s = 0 charts
"""
import argparse
import concurrent.futures as cf
import os
import subprocess
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
FAMILIES = 'A B1 B2 C D E1 E2 F1 F2 G1 G2 G3 G4 H I1 I2 I3 L3 M3 N3 P3 R2'.split()

ap = argparse.ArgumentParser()
ap.add_argument('-n', type=int, default=4, help='parallel Mathematica kernels')
ap.add_argument('-w', type=int, default=6, help='top weight')
ap.add_argument('--families', default=','.join(FAMILIES))
ap.add_argument('--regions', default='z,x,y')
ap.add_argument('--fam-dir', default=os.path.join(HERE, '..', 'families'))
ap.add_argument('--wolframscript', default='wolframscript')
ap.add_argument('--steps', default='edge,div', help='edge, div or both')
ap.add_argument('--dry-run', action='store_true')
a = ap.parse_args()
os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
lock = threading.Lock()
status = {}


def inputs(F):
    return os.path.join(a.fam_dir, F, 'inputs')


def have(F, tag):
    return os.path.exists(os.path.join(inputs(F), 'bdy_%s_w%d.m' % (tag, a.w)))


def others(A):
    return [b for b in 'xyz' if b != A]


def run(F, A, step, cmd, outs):
    key = '%s %s %s' % (F, A, step)
    if all(have(F, t) for t in outs):
        with lock:
            status[key] = 'skipped (done)'
        return True
    log = os.path.join(HERE, 'logs', '%s_%s_%s.log' % (F, A, step))
    t0 = time.time()
    with lock:
        print('[%s] start %s' % (time.strftime('%H:%M:%S'), key), flush=True)
    with open(log, 'w') as fh:
        rc = subprocess.run([a.wolframscript, '-file'] + cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=HERE).returncode
    ok = rc == 0 and all(have(F, t) for t in outs)
    with lock:
        status[key] = '%s in %.0f s' % ('ok' if ok else 'FAILED (rc %d, see %s)' % (rc, log), time.time() - t0)
        print('[%s] %s: %s' % (time.strftime('%H:%M:%S'), key, status[key]), flush=True)
    return ok


fams = [F for F in a.families.split(',') if F]
regs = [A for A in a.regions.split(',') if A]
for F in fams:
    for A in regs:
        if not have(F, 'R%s_th' % A.upper()):
            raise SystemExit('missing threshold constants bdy_R%s_th_w%d.m for %s: run threshold.py consts' % (A.upper(), a.w, F))
if a.dry_run:
    for F in fams:
        for A in regs:
            R = A.upper()
            todo = [t for t in ['R%s_%sv' % (R, b) for b in others(A)] + ['R%s_v%s' % (R, b) for b in others(A)] if not have(F, t)]
            print('%-3s region %s > 1: %s' % (F, A, ('to do: ' + ' '.join(todo)) if todo else 'done'))
    raise SystemExit

pool = cf.ThreadPoolExecutor(max_workers=a.n)
futs = []


def chain(F, A):
    R = A.upper()
    b, c = others(A)
    edge_out = ['R%s_%sv' % (R, b), 'R%s_%sv' % (R, c)]
    if 'edge' in a.steps:
        if not run(F, A, 'edge', ['edge_transport.wl', inputs(F), A, str(a.w)], edge_out):
            return
    if 'div' not in a.steps or not all(have(F, t) for t in edge_out):
        return
    for d in (b, c):
        futs.append(pool.submit(run, F, A, 'div_' + d, ['divisor_transport.wl', inputs(F), A, d, str(a.w)], ['R%s_v%s' % (R, d)]))


t0 = time.time()
heads = [pool.submit(chain, F, A) for F in fams for A in regs]
for h in heads:
    h.result()
while True:                               # divisor jobs were added while the edge jobs finished
    pending = [f for f in futs if not f.done()]
    if not pending:
        break
    cf.wait(pending)
pool.shutdown()
print('\n==== summary (%.0f min) ====' % ((time.time() - t0) / 60))
for k in sorted(status):
    print('  %-14s %s' % (k, status[k]))
bad = [k for k, v in status.items() if v.startswith('FAILED')]
print('%d jobs, %d failed' % (len(status), len(bad)))
