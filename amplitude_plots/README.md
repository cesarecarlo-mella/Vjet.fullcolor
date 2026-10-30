# Amplitude plots

Plots of the full-colour three-loop finite remainders (not of the iterated integrals).

| file | content |
|---|---|
| `region_z_alpha1_uv_re.png`, `region_z_alpha1_uv_im.png` | α1, real / imaginary part, production region z > 1, all 14 colour layers, in (u, v) |
| `plot_region_layers.py` | the script: `plot_region_layers.py --region z --coef alpha1 [--n 30]` (any region and coefficient) |
| `region_z_alpha1_uv.npz` | the evaluated values (cache of the script) |

- Variables: v = 1/z (inverse energy), u = −y/z, w = −x/z = 1 − u − v. v → 1 is threshold, v → 0 high energy.
- Shown: 0.02 ≤ v ≤ 0.3 (z ≥ 3.3), 0.02 from the edges u = 0 and w = 0 (`--vmax` changes the cut).
- Height: the value itself (q² = 1), no rescaling.
- Layers in the order of the colour table: Nc³ first, Nc⁻³ last.
- Data: the double evaluator of the repository, with the region's start values from the continued corner series
  (commit 1043025); 756 points.
- Why the cut: towards threshold (v → 1) the coefficients grow like (1−v)⁻ᵏ (e.g. (x+y)⁻⁵ = −v⁵/(1−v)⁵), and towards
  high energy (v → 0) they go to zero. On a linear scale threshold dominated everything (50–100× the bulk), so it
  is cut. The remaining large values are the collinear edges u → 0 (y → 0) and w → 0 (x → 0).
