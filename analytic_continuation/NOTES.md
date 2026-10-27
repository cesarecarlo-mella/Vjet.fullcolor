# Analytic continuation: how few series cover the production region?

Test bed: the LC finite remainders continued to the production region
(`reproduceLC/2024_Vjet_3L_LC_ancillary/FiniteRemainders/{lmL,lpL,lrP,lrM,rmR,rpR}.m`), which are exact
GPLs. This is planar only (letters 0, 1, 1-v, -v in u and 0, 1 in v), so it tests the line letters but
not the non-planar ones.

## Conventions

The ancillary region is y = s13/q2 > 1 (x, z < 0), with

    v = 1/y,   u = -z/y,   w = -x/y = 1 - u - v,      u, v, w > 0  (a triangle).

The CC work (`series.all/analytic.continuation`) used z > 1 instead: v = 1/z, u = -y/z, w = -x/z. It is
the same with y <-> z, and the corner labels are the same:

| corner | small | in (x, y, z) here | line through the corner (letter) |
|---|---|---|---|
| th (threshold, CC1) | u, w | y -> 1, x, z -> 0 (continued decay corner) | u + w = 0 (1 - v) |
| he2 (high energy, CC2) | u, v | y -> oo at fixed z = -u/v | u + v = 0 (1 - w), i.e. z = 1 |
| he3 (high energy, CC3) | w, v | y -> oo at fixed x = -w/v | w + v = 0 (1 - u), i.e. x = 1 |

## 1. The line through a high-energy corner is a genuine branch line (also at amplitude level)

`xyseries.py` expands iteratively at he2 (u << v first, then v). If the line u + v = 0 were spurious, the
negative powers of v from the towers in u/v would cancel and the lowest power would not depend on the
order. All 18 coefficients (6 channels x 3) give

| order N | lowest power of v |
|---|---|
| 6 | -5 or -6 |
| 10 | -9 or -10 |
| 16 | -15 |

i.e. about -(N-1): the towers survive. The amplitude depends on the ratio z = -u/v along the corner line
(y -> oo at fixed s23), with a branch point at z = 1 (s23 = q2), outside the region but at distance 1
from z = 0. So the ratio series covers only -1 < z < 0, and the old CC strategy needed three charts per
high-energy corner (both orderings and a band chart).

## 2. The square chart at a high-energy corner

    rho = u + v,   s = u/(u+v) in [0, 1]      (u = s rho, v = (1-s) rho)

maps the triangle onto the square [0,1]^2 and sends the branch line to s = oo. The series in (rho, s) is
an exact re-expansion of the iterated (u, v) data (the u << v expansion contains it), so no new
computation is needed.

Result (lmL alpha1, Nc^2 Nf layer, order 16, 406 points, digits vs the exact GPLs, `he_corner.py`):

| representation | median | >= 4 digits | >= 8 digits |
|---|---|---|---|
| ratio (old chart, series in u/v) | -1.9 | 9% | 1% |
| square (rho, s) | 0.7 | 16% | 2% |
| Bernoulli in rho | -1.0 | 21% | 7% |
| Bernoulli in s | -1.1 | 11% | 2% |
| Bernoulli in both | -7.7 | 11% | 3% |

- The square chart fixes the band u ~ v (s = 1/2) where the ratio series fails: e.g. at u = v = 0.1 the
  ratio series gives nothing, the square chart 6 digits.
- But one high-energy series cannot cover the triangle. It converges only on the side u < w. The reason
  is the letter 1 - v = u + w, the line through the threshold corner: at fixed rho its zero sits at
  s = -(1 - rho)/rho, which reaches s = 0 as rho -> 1. A Bernoulli map cannot move it. Bernoulli in s
  makes it worse (the point goes to sigma = log rho, close to 0), and Bernoulli in both turns 1 - u into
  the curve e^beta + e^sigma = 1, which passes through beta = sigma = -log 2.

### A toy with the same behaviour (`toy_example.wl`, Mathematica)

    F(u, v) = Log[1 + u/(v (1 - v))] = Log[u + v - v^2] - Log[v] - Log[1 - v]

Its branch curve u + v - v^2 = 0 starts at the high-energy corner, tangent to u + v = 0 (ratio
dependence, branch point u/v = -1), and ends at the threshold corner (0, 1), outside the triangle in
between. In the square chart F = Log[1 - (1-s)^2 rho] - Log[1 - s] - Log[1 - (1-s) rho]: the branch point
is at s = oo, but the curve comes back at s = 1 - 1/Sqrt[rho], so the series converges only for
(2u + v)^2 < u + v. Order 12, digits vs exact:

| u, v | ratio | square | Bernoulli rho | Bernoulli s | barrier |
|---|---|---|---|---|---|
| 0.05, 0.1 | 4.3 | 6.8 | 6.8 | 5.7 | inside |
| 0.1, 0.1 (u = v) | 0.7 | 4.6 | 4.6 | 2.1 | inside |
| 0.2, 0.1 (u > v) | -2.8 | 2.8 | 2.6 | -1.7 | inside |
| 0.1, 0.5 | 4.4 | 4.3 | 5.3 | 4.1 | inside |
| 0.3, 0.5 | 0.3 | 1.3 | -1.6 | -1.7 | outside |
| 0.4, 0.2 | -3.4 | -0.3 | -1.8 | -4.7 | outside |
| 0.5, 0.1 | -7.5 | -1.3 | -2.8 | -7.2 | outside |

Ratio series fails at u > v, the square chart fixes the band and u > v near the corner, it stops at the
barrier, and Bernoulli in s makes it worse: the same pattern as the amplitude. A plain product such as
Log[1 + u/v] Log[1 - v] does NOT show the barrier: a function analytic in (u, v) has polynomial
coefficients in s at each power of rho, so the truncation in s costs nothing. The barrier needs the
branch structure to connect the two corners.

## 3. Three corners (in progress)

Each corner series is limited by the lines through the other two corners. Expected:
- he2 covers u < w, and he3 w < u (by symmetry);
- th is a clean double series in (w, u): its own line u + w = 0 should be spurious there, as at the decay
  corner. With Bernoulli in w and u the lines w = 1 and u = 1 go away, so it should cover the whole
  triangle and be slow only near the high-energy edge v = 0;
- the middle of the edge v = 0 (u ~ w, v -> 0), far from both high-energy corners, is the weakest place.

th and he3 need the coefficient in their own fibration: `refibrate.wls` (PolyLogTools `ToFibrationBasis`,
v = 1 - u - w with {w, u}, or u = 1 - v - w with {w, v}). In both cases the letters are again
{0, 1, 1-y, -y} and {0, 1}, so `xyseries.py` expands them unchanged.

## 3b. Result: five series cover the whole triangle

Charts (here in the LC region; the CC labels are the same):
- **th**: threshold corner, double series in (w, u), Bernoulli in both (w -> 1 - e^-bw, u -> 1 - e^-bu);
- **he2 s=0 / he2 s=1**: corner u, v -> 0, square chart p = u + v with s = u/(u+v) around s = 0 (u << v),
  and around s = 1 (v << u), both with Bernoulli in p only;
- **he3 s=0 / he3 s=1**: the same at the corner w, v -> 0 (u <-> w).

Test: LC remainders, colour layer Nc Nf^2 (GPLs up to weight 4, refibration ~1 min per corner), order 16,
406 points of the triangle, digits vs the exact GPLs:

| | lmL alpha1: median, >= 4, >= 8, worst | lrP alpha1: median, >= 4, >= 8, worst |
|---|---|---|
| th alone | 7.6, 71%, 47%, -5.8 | 8.1, 76%, 51%, -6.4 |
| three corners (th + both s = 0 charts) | 9.2, 86%, 60%, 1.4 | 9.2, 88%, 59%, 0.1 |
| **five series** | **10.2, 100%, 80%, 5.6** | **9.7, 100%, 65%, 4.9** |

- The threshold series covers the whole upper part of the triangle (v >~ 0.25) at >= 9 digits: Bernoulli
  works there as in the decay region (plain: median 5.0, Bernoulli 7.6).
- The two s = 1 charts are essential: they cover the high-energy edge v -> 0, which no other series
  reaches (worst point 0-1 digits without them, ~5 with them).
- At the high-energy corners Bernoulli in p (only) is the best representation: with the plain square
  charts the five series give median 8.6-8.9 and worst 2.9-3.8.
- Weakest region at order 16: around v ~ 0.15-0.25, u ~ 0.15-0.3 (5-6 digits), where th and both charts
  of the (u, v) corner meet. Higher order improves it.
- The Nc^2 Nf layer (weight 5) needed ~30 min per refibration with ToFibrationBasis and was stopped;
  the test is about the convergence domains, which depend on the letters, not on the layer.

## 3c. A fixed (democratic) split of the triangle (`democratic.py`)

A fixed geometric split replaces the point-by-point best of five. It is symmetric under u <-> w:
- **threshold** series (Bernoulli in w, u) where v >= lambda max(u, w);
- elsewhere, for u <= w, the (u, v) corner: the s = 0 chart if u <= v and the s = 1 chart if u > v;
  for w < u, the same at the (w, v) corner.

Every chart stays inside its convergence domain (for u <= w, v < lambda w < w). In the (s, p) square the
threshold cut is p = lambda/(1 - s + lambda) for u <= w, and the vertical line s = 1/(1 + lambda) for
u > w. The high-energy splits are s = 1/2 and p = 1/(1 + s).

Scan of the cut, worst digits over the four coefficients below (triangle / (s, p) grid):

| cut | band v >= 0.15 | band 0.25 | lambda 0.40 | **lambda 0.45** | lambda 0.50 | lambda 0.60 |
|---|---|---|---|---|---|---|
| worst | -1.1 / -1.2 | 1.4 / 1.2 | 2.1 / 2.0 | **2.3 / 2.1** | 2.0 / 1.8 | 1.4 / 1.2 |

The best horizontal band v >= v0 depends on the coefficient: v0 = 0.15 for lpL gamma1, 0.25-0.30 for the
others. One lambda works for all four. **Choice: lambda = 0.45.**

Order 16, Nc Nf^2 layer, digits vs the exact GPLs (median, >= 4, >= 8, worst):

| coefficient | triangle (406 points) | (s, p) grid (841 points) |
|---|---|---|
| lmL alpha1 | 10.1, 98%, 72%, 2.3 | 10.2, 99%, 72%, 2.8 |
| lrP alpha1 | 9.7, 99%, 64%, 3.4 | 9.8, 99%, 68%, 2.1 |
| lpL gamma1 | 9.6, 98%, 76%, 2.8 | 9.6, 98%, 73%, 2.5 |
| lrM gamma1 | 9.7, 100%, 64%, 4.1 | 8.8, 99%, 59%, 3.2 |

Medians per region: threshold 10.3-11.7, s = 0 charts 7.6-9.9, s = 1 charts 6.1-10.6. Compared with the
best of five (3b), the median is the same and the worst point loses 1.5-3 digits. Weak spots:
- the junction where th and both charts of a corner meet, u ~ v ~ 0.2-0.33 (w ~ 0.5): 4-5 digits, and
  2.8 for lpL gamma1 (s = 1 chart at (u, v, w) = (0.33, 0.20, 0.47)). Here s ~ 1/2, the limit of both
  charts in s, so only higher order helps;
- the threshold region along its cut at large u, small w: lmL alpha1 gives 2.3 at (0.67, 0.30, 0.03);
- the (s, p) worst for lrP alpha1, 2.1 at (0.48, 0.02, 0.50), is relative: the Nc Nf^2 layer is
  suppressed there (|exact| ~ 3.5e-4), so the absolute error is ~3e-6.

The order-16 figures have been removed; the order-18 ones below supersede them.

### Order 18, digits measured against the size of the function

The coefficients are not O(1) everywhere. Over the triangle (Nc Nf^2 layer) the median |f| is 0.7-3.2,
the minimum 4e-4 to 4e-3 and the maximum 20-1200 (log growth towards the edges). Digits are now
-log10 |S - E| / max(|f|, median |f|): relative where f is large, absolute in units of its typical size
where it dips (`democratic.py ... scaled`, from the `err` field of `he_corner.py`). When choosing the rule,
lattice points on a boundary are scored on the worse side.

- Picking one series per point, the same for all four coefficients (best worst case), already gives
  >= 6 digits everywhere at order 18: worst 6.2 (triangle), 5.9 ((s, p) grid).
  See `figures/common_best_NcNf2_N18.png`.
- Rule refined to match (`cut:0.5,-0.2,0.55`): threshold where v >= 0.5 max(u, w) - 0.2 min(u, w)
  (v = 0.33 at the edges, 0.13 at u = w; the threshold series beats the corner charts at the junction).
  The s split moves to s0 = 0.55, because at s = 1/2 the s = 0 charts beat the s = 1 charts.
  This supersedes lambda = 0.45, which has worst 4.7 at order 18.

| coefficient | order 16, triangle worst | order 18, triangle: median, >= 6, worst | order 18, (s, p): median, >= 6, worst |
|---|---|---|---|
| lmL alpha1 | 4.2 | 12.2, 99%, 5.7 | 12.4, 99%, 5.4 |
| lrP alpha1 | 4.3 | 11.8, 99%, 5.8 | 12.4, 99%, 5.5 |
| lpL gamma1 | 5.0 | 11.6, 99%, 6.0 | 11.7, 99%, 5.9 |
| lrM gamma1 | 4.3 | 11.8, 99%, 5.8 | 11.4, 99%, 5.7 |

At the weakest points the gain is 0.45-0.75 digits per order (orders 16 to 18). Extrapolating linearly,
6 digits everywhere needs order ~19 and 7 digits order ~20-21. The weakest points are the threshold series
next to the edge w -> 0 (u -> 0 for lrM) at v ~ 0.33, and the junction near (0.27, 0.20, 0.53).
Figures: `figures/democratic_<coefficient>_NcNf2_N18.png`.

**Near the threshold corner (v -> 1) the measured digits stop at 12-13 because of the reference, not the
series.** At (u, v, w) = (1/30, 14/15, 1/30), evaluated at 30 digits:
- the threshold series at orders 16 and 18 agree to 1e-27;
- its double-precision value is good to ~1e-15 relative;
- the "exact" value is off by ~1e-13 (lmL) to 2e-12 (lrP) relative.

Two reasons. On Apple Silicon numpy `longdouble` is plain double, while `gpl_eval.py` was validated to 1e-17
with 80-bit long double. And the sum of rational x GPL terms cancels by a factor 30-600 there, because of the
spurious 1/(1-v)^k poles. This doesn't matter for the 6-7 digit target. For a cleaner reference, use GiNaC
at higher precision or an x86 machine.

## 4. Boundary constants at s = 1: transport along the divisor (`divisor_transport.wls`)

The s = 1 charts (v << u at the corner u, v -> 0, the old CC2e ordering) need their own constants. They
follow from the s = 0 ones (CC2, known at all weights from `transport.boundaries`) by transport along
the exceptional divisor of the corner. In the CC chart (x = -w/v, y = -u/v, z = 1/v) with
u = s p, v = (1-s) p, w = 1 - p, as p -> 0 every letter becomes p^k g_k(s):

| letters | on the divisor |
|---|---|
| l1-l6, l11-l16, l18, l20 | (1 - s) |
| l2, l9, l19 | s, 1 - s |
| l10, l17 | 1 - 2s (the band y = -1) |
| l8 (odd, r = -Sqrt|xyz|) | Log[(q + I)/(q - I)], q = Sqrt[s/(1-s)] |
| l7 (odd) | -> 0 |

The log-free part F0(s) obeys dF0/ds = A(s) F0 with A = Sum_k atilde_k d/ds Log[g_k] (the powers of p only
give the Log[p] tower). With s = q^2/(1+q^2) (s in [0, 1) <-> q in [0, oo)) the alphabet is {0, +-1, +-I}
and F0 is an iterated integral from the tangential base point q = 0 (the CC2 constants), weight by
weight. The s = 1 constants are the regularised limit q -> oo (Log[q] -> 0, i.e. Log[1 - s] -> 0).

The limit q -> oo is exact: with q = 1/r, dq/(q - a) = -dr/r + dr/(r - 1/a) (a != 0) and dq/q = -dr/r, so
T(w; r) = G(w; 1/r) obeys T(a, w'; r) = G(a, w'; 1) - S(1) + S(r) with S(r) = Int_0^r k_a T(w'; r') dr', and
the regularised limit (Log[r] -> 0, every other G(..; r) -> 0) is R(a, w') = G(a, w'; 1) - S(1). All
constants are G(...; 1) with letters {0, +-I} (closed under a -> 1/a), i.e. HPLs at I.

**H family, weights 1 and 2: the transported constants equal the fitted CC2e constants
(`cc_series/inputs/bdy_CC2e_w{1,2}.m`, recognised from the exact continued solution) exactly, in all
371 components.** At weight 2, 119 components differ between the CC2 and CC2e orderings (e.g.
-11 Pi^2/18 -> -13 Pi^2/18), and all are reproduced. The band letters +-1 cancel: only {0, +-I} survive,
so y = -1 is spurious (at least to weight 2). Output: `results/bdy_CC2e_H_w<w>_transported.m`.
The corner u -> 1 is the same with u <-> w.

Weights 3 and 4 (no reference): Catalan and Log[2] appear in intermediate steps and cancel in every
component; the constants are pure i Pi / Zeta[3] combinations, as a genuine boundary value should be.

**The band y = -1 at weights >= 5.** The letter l10 = 1 - x(1+y) enters the H atilde in a single column,
master 182, feeding 15 rows. So the band letter (q = +-1, s = 1/2) appears in the transported functions as
soon as F_182 is nonzero: its CC2 constants vanish at weights 0-3 and not from weight 4 on, hence band
words from weight 5. The 15 rows are regular at y = -1 iff F_182 vanishes there, weight by weight:
F_182^(4)(q) is a sum of 17 GPLs with Pi^4, i Pi^3, ... coefficients, and **it vanishes at q = 1 to 15
digits**. So the band stays spurious at weight 5 (a non-trivial check of the CC2 constants that feed
F_182). The words with letter 1 then make G(1, ...; 1) individually divergent at the split point q = 1 of
the q -> oo limit; the total is regular, so PolyLogTools `ShuffleRegulate` with HPL[{1}, 1] -> 0 is exact.
**The coefficients of the divergence Log[1 - 1] cancel exactly in every component at weights 1-6: the
band is regular (spurious) through weight 6.**

**Result: the CC2e constants of H at weights 1-6** (`results/bdy_CC2e_H_w<w>.m`, by
`divisor_transport.wls 6` + `finalise_cc2e.wls`, about 17 min on one kernel):
- weights 1-2: equal to the fitted constants, all 371 components;
- every component at every weight is a rational combination of i Pi, Zeta[3], Zeta[5] (MZVs): the
  level-4 constants (Catalan, Log[2], Hurwitz Zeta[5, 1/4], Zeta[5, 3/4]) that appear in intermediate steps
  cancel. The 15 (w5) + 23 (w6) components that PolyLogTools left as G(...; 1) were fitted on the MZV
  basis (GiNaC values at 75 digits, 3-term relations, residual < 1e-60), e.g. master 268 at weight 6:
  -22360469/8709120 Pi^6 + 12947/48 Zeta[3]^2 + I (4465/144 Pi^3 Zeta[3] - 13569/20 Pi Zeta[5]).

This replaces the open problem of `cc_series/CONSTANTS.md` (CC2e/CC3e constants at weights >= 3) for H.
The same script works for every family and for the corner u -> 1 (u <-> w), given that family's CC2/CC3
constants at s = 0 (from `transport.boundaries`, to be generalised beyond H).

## Files

- `to_xy.wls`: LC remainders -> `xyseries.py` input, (u, v) -> (x, y); q2 -> 1; I -> ii.
- `refibrate.wls`: one coefficient and colour layer in the fibration of the th or he3 corner.
- `xyseries.py`, `gpl_eval.py`: copied unchanged from `reproduceLC` (iterated GPL expansion; fast
  numerical GPLs, ~1e-17).
- `he_corner.py`: one corner's series in its representations vs the exact GPLs on a lattice of the triangle.
- `coverage.py`: best of several corner series, point by point, with the map on the triangle.
- `plot_convergence.py`: digit maps of the five series and their best, on the (s, p) square and the triangle.
- `democratic.py`: the fixed split of 3c (`bary:<lambda>` or `band:<v0>`, or a scan), stats and maps.
