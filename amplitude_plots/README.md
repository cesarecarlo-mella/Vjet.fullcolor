# Amplitude plots

Plots of the full-colour three-loop finite remainders (not of the iterated integrals).

| file | content |
|---|---|
| `region_z_alpha1_re.png`, `region_z_alpha1_im.png` | α1, real / imaginary part, production region z > 1, all 14 colour layers, over the region's triangle (X, Y) |
| `region_z_alpha1_phys_re.png`, `region_z_alpha1_phys_im.png` | the same in physical variables: energy log10 z (z = s23/q², 1.05 to 1000) against the collinear variable ξ = x/(x+y) = s12/(s12+s13) |

- Height: asinh(value/s), with s the median |value| of the layer (linear near 0, logarithmic in the tails).
- Data: 900 points from the double region evaluator, made on 2026-10-29 before the quad start values of the regions.
- The script that made these was a one-off and is not in the repository; to be redone with a saved script.
