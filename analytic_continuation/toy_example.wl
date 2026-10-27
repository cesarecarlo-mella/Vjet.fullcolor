(* ::Package:: *)
(* toy_example.wl -- a one-line function with the singularity structure that matters for the
   expansions in the production triangle.  Run it cell by cell in a notebook.

   Triangle: v = 1/y, u = -z/y, w = 1 - u - v  (u, v, w > 0).
   High-energy corner he2: u, v -> 0 (w -> 1).  Threshold corner th: v -> 1 (u, w -> 0).

     F(u, v) = Log[1 + u/(v (1 - v))] = Log[u + v - v^2] - Log[v] - Log[1 - v]

   Its branch curve u + v - v^2 = 0 STARTS at the high-energy corner, where it is tangent to the line
   u + v = 0 (so F depends on the ratio u/v there, as the amplitude does: branch point u/v = -1), and
   ENDS at the threshold corner (u, v) = (0, 1).  In between it stays outside the triangle (u < 0).
   F is real inside the triangle. *)

F[u_, v_] := Log[1 + u/(v (1 - v))];
n = 12;                                         (* expansion order *)

(* ---- 1. the old chart at he2: iterated expansion, u << v first (a series in the ratio u/v) ---- *)
serRatio = Normal[Series[Normal[Series[F[u, v], {u, 0, n}]], {v, 0, n}]];
(* contains (u/v)^k: converges only for |u/v| < 1, i.e. u < v *)

(* ---- 2. the square chart: rho = u + v, s = u/(u + v)  (u = s rho, v = (1 - s) rho) ---- *)
Fsq[rho_, s_] := F[s rho, (1 - s) rho];
Fsq2[rho_, s_] := Log[1 - (1 - s)^2 rho] - Log[1 - s] - Log[1 - (1 - s) rho];
Chop[Fsq[0.3, 0.4] - Fsq2[0.3, 0.4]]                     (* the same function: 0 *)
(* the branch point u/v = -1 is now at s = oo (only Log[1 - s] is left of it), but the branch curve
   comes back at s = 1 - 1/Sqrt[rho] (negative, -> 0 as rho -> 1): the double series converges for
   rho (1 + s)^2 < 1, i.e. (2u + v)^2 < u + v.  One high-energy series cannot reach the far side. *)
serSquare = Normal[Series[Fsq[rho, s], {rho, 0, n}, {s, 0, n}]];

(* ---- 3. Bernoulli variables in the square chart: rho = 1 - E^-beta, s = 1 - E^-sigma ---- *)
serBernRho = Normal[Series[Fsq[1 - E^-beta, s], {beta, 0, n}, {s, 0, n}]];
serBernS = Normal[Series[Fsq[rho, 1 - E^-sigma], {rho, 0, n}, {sigma, 0, n}]];

(* ---- 4. the threshold corner: series in (w, u) with v = 1 - u - w ---- *)
Fth[w_, uu_] := F[uu, 1 - uu - w];
serTh = Normal[Series[Normal[Series[Fth[w, uu], {w, 0, n}]], {uu, 0, n}]] /. uu -> u;
(* here the curve passes THROUGH the corner (it is the toy's analogue of a line through th), so the
   toy also needs the ratio there; in the amplitude that line should be spurious at th *)

(* ---- digits against the exact function ---- *)
dig[ser_, pt_] := Min[16, -Log10[Abs[(ser /. pt) - F[u, v] /. pt] / Abs[F[u, v] /. pt] + 10^-30]];
toPts[u0_, v0_] := {u -> u0, v -> v0, w -> 1 - u0 - v0, rho -> u0 + v0, s -> u0/(u0 + v0),
   beta -> -Log[1 - u0 - v0], sigma -> -Log[1 - u0/(u0 + v0)]};
pts = {{1/20, 1/10}, {1/10, 1/10}, {1/5, 1/10}, {1/10, 1/2}, {1/5, 1/2}, {3/10, 1/2}, {2/5, 1/5}, {1/2, 1/10}};
TableForm[Table[Join[p, N[dig[#, toPts @@ p] & /@ {serRatio, serSquare, serBernRho, serBernS}, 3],
     {If[(2 p[[1]] + p[[2]])^2 < Total[p], "inside", "outside"]}], {p, pts}],
  TableHeadings -> {None, {"u", "v", "ratio", "square", "bern rho", "bern s", "(2u+v)^2 < u+v"}}]

(* ---- maps on the triangle ---- *)
mapOf[ser_, title_] := DensityPlot[dig[ser, toPts[uu, vv]], {uu, 0.01, 0.98}, {vv, 0.01, 0.98},
   RegionFunction -> Function[{uu, vv}, uu + vv < 0.99], PlotPoints -> 60, PlotRange -> {0, 16},
   ColorFunction -> "TemperatureMap", PlotLegends -> Automatic, FrameLabel -> {"u", "v"}, PlotLabel -> title,
   Epilog -> {Dashed, Line[{{0, 0}, {1/2, 1/2}}],
     Line[Table[{uu, vv /. Last[Solve[(2 uu + vv)^2 == uu + vv && vv > 0, vv]]}, {uu, 0.001, 0.25, 0.005}]]}];
(* dashed: u = v (where the ratio series stops) and (2u + v)^2 = u + v (where the square chart stops) *)
GraphicsRow[{mapOf[serRatio, "ratio series (old chart)"], mapOf[serSquare, "square chart (rho, s)"],
  mapOf[serBernRho, "Bernoulli in rho"], mapOf[serBernS, "Bernoulli in s"]}, ImageSize -> 1200]
