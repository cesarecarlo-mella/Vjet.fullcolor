"""remainder_eval vs the exact values of the leading-colour layer (C3F0Z0 = Nc^3, HA1-6) at the 406 grid points of
reproduceLC/HA_series_Nc3_vs_exact_digits (GPL form, >= 15 digits, cross-checked with GiNaC and the 2024 ancillary)."""
import json, os, time
import numpy as np
from remainder_eval import evaluate_remainders, LABELS

HERE = os.path.dirname(os.path.abspath(__file__))
ref = json.load(open(os.path.join(HERE, '..', '..', 'reproduceLC', 'HA_series_Nc3_vs_exact_digits', 'results',
                                  'digits_order18.json')))['rows']
pts = [(r['x'], r['y']) for r in ref]
E = np.array([[complex(float(a), float(b)) for a, b in r['exact']] for r in ref])
t0 = time.time()
R, info, err = evaluate_remainders(pts, 1e-12, error=True)
print('%d points in %.1f s, with the error estimate (median %.3f s per point and evaluation)'
      % (len(pts), time.time() - t0, np.median([i['time'] for i in info])))
cols = [LABELS.index('C3F0Z0_HA%d' % h) for h in range(1, 7)]
d = np.abs(R[:, cols] - E) / np.abs(E)
est = err[:, cols] / np.abs(E)
big = d > 1e-11                                       # where the error is above the noise of the reference
print('error estimate vs actual error, where the actual error is > 1e-11 (%d values): estimate >= actual / 3 for %.1f%%,'
      ' median ratio %.2f' % (big.sum(), 100 * np.mean(est[big] >= d[big] / 3), np.median(est[big] / d[big])))
for h in range(6):
    i = int(d[:, h].argmax())
    print('HA%d: relative difference median %.1e, 90%% %.1e, worst %.1e at (%.3f, %.3f)'
          % (h + 1, np.median(d[:, h]), np.quantile(d[:, h], 0.9), d[i, h], *pts[i]))
print('other path deformation needed at:', [(i['x'], i['y'], i['delta']) for i in info if i['delta'] != (0.1, 0.2)])
