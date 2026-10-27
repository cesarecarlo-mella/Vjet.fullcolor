#!/usr/bin/env python3
# =============================================================================
#  run_all_families.py -- blow-up series of EVERY family (weights 1..W, all
#  corners) with the array/modular recursion, one family after the other,
#  cheapest first (ordered by the size of families/<F>/inputs/atilde.m).
#
#  The output stays in blow-up variables (t, BIG): no phys conversion, since
#  the helicity coefficients are multiplied in these variables before going
#  back to (x, y).
#
#  Layout (self-contained -- copy the whole folder to the cluster):
#     code/                    the pipeline (run_array.py, 00/01/02, engine, build_letters.py)
#     letters/                 dlog series of the 20 letters per chart (build_letters.py)
#     families/<F>/inputs/     atilde.m, sol0.m, bdy_<tag>_w<w>.m
#     families/<F>/dist_<T>_<Y>/   results: exact_<tag>_w<w>.m, rec_<tag>_w<w>.pkl, logs/
#  Charts (code/charts.py): 1-3 decay corners, 4-6 their mirrors (tags c1..c6), and the
#  continued-region charts R<A>_<b>v, R<A>_v<b> (A > 1), e.g. -c RZ_yv,RZ_vy,RZ_xv,RZ_vx
#
#  Families are independent: a failure is reported and the next family runs.
#  Everything is resumable (finished weights / corners / families are skipped).
#
#  USAGE
#     ./run_all_families.py -n 64 --ordt 22 --ordy 66 -w 6
#     ./run_all_families.py -n 64 --ordt 22 --ordy 66 -w 6 --only M3,N3,H
#     ./run_all_families.py -n 64 --ordt 22 --ordy 66 -w 6 --resume-from E1
#     ./run_all_families.py --dry-run --ordt 22 --ordy 66
#  Any extra option after "--" is passed to run_array.py unchanged, e.g.
#     ./run_all_families.py -n 64 --ordt 22 --ordy 66 -- --nprimes 60
# =============================================================================
import argparse
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.join(HERE, "code")
LETTERS = os.path.join(HERE, "letters")
FAMILIES = os.path.join(HERE, "families")

argv = sys.argv[1:]
extra = []
if "--" in argv:
    i = argv.index("--"); extra = argv[i + 1:]; argv = argv[:i]
ap = argparse.ArgumentParser()
ap.add_argument("-n", "--ncores", type=int, default=os.cpu_count())
ap.add_argument("-w", "--weight", type=int, default=6)
ap.add_argument("-c", "--corners", default="1,2,3,4,5,6",
                help="chart ids of code/charts.py: 1-3 standard, 4-6 their mirrors, or continued charts such as RZ_yv")
ap.add_argument("--ordt", type=int, required=True)
ap.add_argument("--ordy", type=int, required=True)
ap.add_argument("--phys", type=int, default=None, help="also write phys_..._ord<ORD>.m (only for checks)")
ap.add_argument("--keep-residues", action="store_true",
                help="keep the prime residues (needed only to resume or to extend to weight 7)")
ap.add_argument("--only", default=None, help="comma list of families")
ap.add_argument("--resume-from", default=None, help="skip families before this one (in run order)")
ap.add_argument("--python", default=sys.executable)
ap.add_argument("--dry-run", action="store_true")
a = ap.parse_args(argv)

# fail fast (instead of once per family) if this interpreter lacks the packages
_chk = subprocess.run([a.python, "-c", "import numpy, gmpy2"], capture_output=True, text=True)
if _chk.returncode != 0:
    sys.exit("ERROR: %s cannot import numpy/gmpy2:\n  %s\n"
             "Run with the venv's interpreter, e.g.  venv/bin/python run_all_families.py ..."
             % (a.python, _chk.stderr.strip().splitlines()[-1] if _chk.stderr else "?"))

# the letters must reach the requested order (terms above it are dropped on load)
sys.path.insert(0, CODE)
import charts
for c in charts.parse_list(a.corners):
    f = os.path.join(LETTERS, "info_%s.m" % charts.tag(c))
    if not os.path.isfile(f):
        sys.exit("ERROR: no %s -- build the letters: python3 code/build_letters.py ORDT ORDY --corners %s" % (f, c))
    info = open(f).read()
    lt, ly = (int(re.search(r'"%s" -> (\d+)' % k, info).group(1)) for k in ("ORDt", "ORDy"))
    if lt < a.ordt or ly < a.ordy:
        sys.exit("ERROR: letters %s are built to %d/%d < requested %d/%d -- rebuild them with\n"
                 "  python3 code/build_letters.py %d %d --corners %s"
                 % (charts.tag(c), lt, ly, a.ordt, a.ordy, a.ordt, a.ordy, c))

fams = sorted((d for d in os.listdir(FAMILIES)
               if os.path.isfile(os.path.join(FAMILIES, d, "inputs", "atilde.m"))),
              key=lambda d: os.path.getsize(os.path.join(FAMILIES, d, "inputs", "atilde.m")))
if a.only:
    want = [f.strip() for f in a.only.split(",") if f.strip()]
    missing = [f for f in want if f not in fams]
    if missing:
        sys.exit("unknown families: %s  (have: %s)" % (missing, fams))
    fams = [f for f in fams if f in want]
if a.resume_from:
    if a.resume_from not in fams:
        sys.exit("--resume-from %s: not in %s" % (a.resume_from, fams))
    fams = fams[fams.index(a.resume_from):]

wd_name = "dist_%d_%d" % (a.ordt, a.ordy)
print("=== ALL FAMILIES  weights 1..%d | corners %s | order %d/%d | N = %d ===" %
      (a.weight, a.corners, a.ordt, a.ordy, a.ncores))
for f in fams:
    print("   %-4s  atilde %7.1f kB   -> families/%s/%s" %
          (f, os.path.getsize(os.path.join(FAMILIES, f, "inputs", "atilde.m")) / 1e3, f, wd_name))
if a.dry_run:
    sys.exit(0)

ok, bad = [], []
T0 = time.time()
for f in fams:
    t0 = time.time()
    cmd = [a.python, os.path.join(CODE, "run_array.py"), "-w", str(a.weight), "-n", str(a.ncores),
           "-c", a.corners, "--ordt", str(a.ordt), "--ordy", str(a.ordy),
           "--workdir", os.path.join(FAMILIES, f, wd_name),
           "--inputs", os.path.join(FAMILIES, f, "inputs"), "--letters", LETTERS]
    if a.phys is not None:
        cmd += ["--phys", str(a.phys)]
    if not a.keep_residues:
        cmd += ["--clean-residues"]
    cmd += extra
    print("\n################  FAMILY %s  (%s)  ################" % (f, time.strftime("%H:%M:%S")), flush=True)
    rc = subprocess.run(cmd).returncode
    (ok if rc == 0 else bad).append(f)
    print("################  %s %s  (%.1f min)" % (f, "OK" if rc == 0 else "FAILED rc=%d" % rc,
                                                 (time.time() - t0) / 60), flush=True)

print("\n================================================================")
print("  done in %.1f min   OK: %s" % ((time.time() - T0) / 60, " ".join(ok) or "-"))
print("  FAILED: %s" % (" ".join(bad) or "-"))
print("================================================================")
sys.exit(0 if not bad else 1)
