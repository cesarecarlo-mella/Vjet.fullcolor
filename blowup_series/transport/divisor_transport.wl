(* divisor_transport.wl -- boundary constants of the s = 1 chart R<A>_v<b> of a continued region (one family):
   transport along the exceptional divisor of the high-energy corner bh, v -> 0, from the s = 0 chart R<A>_<b>v.
   Generalises analytic_continuation/divisor_transport.wls + finalise_cc2e.wls (H, A = z, b = y: R<Z>_vy = CC2e).

   Region A > 1 (charts.py): A = 1/v, b = -bh/v, c = -ch/v, v + bh + ch = 1.  Blow-up of the corner:
       bh = s p ,  v = (1 - s) p ,  ch = 1 - p ;   s -> 0: bh << v (R<A>_<b>v),  s -> 1: v << bh (R<A>_v<b>).
   On the divisor p -> 0 every even letter becomes p^k g_k(s); an odd letter (a - I r)/(a + I r) with
   rho = r/a, r = -Sqrt|xyz| (continued region), tends to +-1 (form 0) or to (1 - I rho0)/(1 + I rho0) with rho0(s)
   algebraic.  With s = q^2/(1 + q^2) all forms are rational in q (checked).  The log-free part obeys
       F^(w)(q) = Sum_a M_a . Int_0^q dq'/(q' - a) F^(w-1)(q') + C^(w)_{s=0}      (M_a = Sum_k atilde_k res_a dlog_k)
   from the tangential base point q = 0; the s = 1 constants are its regularised limit q -> oo (Log[q] -> 0),
   taken exactly by the inversion q = 1/r (all constants become G(...; 1)).  Divergent G(1, ...; 1) are
   shuffle-regularised (their coefficients must cancel: checked).  Constants that PolyLogTools leaves as G(...; 1),
   Catalan, Log[2], ... are fitted on Q[i Pi, MZV] (GiNaC values, 75 digits, accepted below 1e-60).

   usage: wolframscript -file divisor_transport.wl <family inputs dir> <region x|y|z> <b> [WMAX = 6]
   reads <inputs>/bdy_R<A>_<b>v_w<w>.m (edge_transport.wl), writes <inputs>/bdy_R<A>_v<b>_w<w>.m *)
SetEnvironment["PATH" -> GetEnvironment["PATH"][[2]] <> ":/usr/local/bin:/opt/homebrew/bin"];
Get["PolyLogTools`"];                     (* always first *)
argv = Rest[$ScriptCommandLine];
inp = argv[[1]]; region = argv[[2]]; bb = argv[[3]]; wmax = If[Length[argv] > 3, ToExpression[argv[[4]]], 6];
here = DirectoryName[$InputFileName];
t0 = AbsoluteTime[];
alph = Get[FileNameJoin[{here, "alphabet.wl"}]];          (* copy of data/trident_definition/alphabet.wl *)
W = Association[Table[k -> (l[k] /. alph /. Log[a_] :> a), {k, 20}]];
atilde = Get[FileNameJoin[{inp, "atilde.m"}]];
nf = Length[atilde];
R = ToUpperCase[region];
seedTag = "R" <> R <> "_" <> bb <> "v"; outTag = "R" <> R <> "_v" <> bb;
seed[w_] := Get[FileNameJoin[{inp, "bdy_" <> seedTag <> "_w" <> ToString[w] <> ".m"}]];
Print["family dir ", inp, ": ", seedTag, " -> ", outTag, ", ", nf, " masters, weights <= ", wmax];

(* ---- letters on the divisor, as dlogs in q ---- *)
vars = {"x", "y", "z"};
others = DeleteCases[vars, region];
If[! MemberQ[others, bb], Print["!! b must be one of ", others]; Exit[1]];
cc = First[DeleteCases[others, bb]];
hat = <|bb -> s p, cc -> 1 - p|>;
vv = (1 - s) p;
val = Association[region -> 1/vv, Table[b -> -hat[b]/vv, {b, others}]];
X = val["x"]; Y = val["y"];
sgn = Association[region -> 1, Table[b -> -1, {b, others}]];      (* signs of x, y, z in the region *)
sOfq = q^2/(1 + q^2);
lowest[e_, v_] := Module[{n = Numerator[e], d = Denominator[e]},
   {Exponent[n, v, Min] - Exponent[d, v, Min], Together[Coefficient[n, v, Exponent[n, v, Min]]/Coefficient[d, v, Exponent[d, v, Min]]]}];
dlogq = Association[Table[k -> If[MemberQ[{7, 8}, k],
     Module[{a = If[k == 7, x (1 - x - y), x y], sa = If[k == 7, sgn["x"] sgn["z"], sgn["x"] sgn["y"]], r2, rho0},
       r2 = lowest[Together[(x y (1 - x - y)/a^2 /. {x -> X, y -> Y})], p];
       If[r2[[1]] != 0, 0,
         rho0 = -sa Simplify[Sqrt[Factor[Together[r2[[2]] /. s -> sOfq]]], q > 0];
         If[! FreeQ[rho0, Power[_, _Rational] | Abs | Sqrt], Print["!! odd letter l", k, " not rational in q: ", rho0]; Exit[1]];
         Print["  l", k, " on the divisor: (1 - I rho)/(1 + I rho), rho = ", rho0];
         Apart[Together[D[Log[(1 - I rho0)/(1 + I rho0)], q]], q]]],
     Module[{g = lowest[Together[W[k] /. {x -> X, y -> Y}], p][[2]]},
       Apart[Together[D[Log[g /. s -> sOfq], q]], q]]],
   {k, 20}]];
lets = Union[Flatten[Table[If[dlogq[k] === 0, {}, q /. Solve[Denominator[Together[dlogq[k]]] == 0, q]], {k, 20}]]];
coefOf[f_, a_] := Together[Limit[(q - a) f, q -> a]];
dl = Association[Table[k -> Association[Table[a -> If[dlogq[k] === 0, 0, coefOf[dlogq[k], a]], {a, lets}]], {k, 20}]];
Do[If[Simplify[dlogq[k] - Total[Table[dl[k][a]/(q - a), {a, lets}]]] =!= 0,
     Print["!! dlog of l", k, " is not a sum of simple poles"]; Exit[1]], {k, 20}];
Print["alphabet in q: ", lets];
Do[If[dlogq[k] =!= 0, Print["  l", k, ": ", Total[Table[dl[k][a]/(q - a), {a, lets}]]]], {k, 20}];

(* ---- M_a = Sum_k atilde_k res_a(dlog_k) ---- *)
Mc = Association[Table[a -> SparseArray[atilde /. l[k_] :> dl[k][a]], {a, lets}]];

(* ---- transport: F^(w) as vectors of G[..., q] ---- *)
integrate[vec_, a_] := Module[{v = Expand[vec]},
   Map[Function[e, If[e === 0, 0,
       Total[Map[If[FreeQ[#, G], # G[a, q], # /. G[ws__, q] :> G[a, ws, q]] &, If[Head[e] === Plus, List @@ e, {e}]]]]], v]];
F[0] = seed[0];
Put[F[0], FileNameJoin[{inp, "bdy_" <> outTag <> "_w0.m"}]];
Do[
  F[w] = Expand[Sum[Mc[a] . integrate[F[w - 1], a], {a, lets}] + seed[w]];
  Print["weight ", w, ": ", Length[Union[Cases[F[w], _G, Infinity]]], " distinct G  (", Round[AbsoluteTime[] - t0, 0.1], " s)"],
  {w, 1, wmax}];

(* ---- regularised limit q -> oo by inversion q = 1/r ---- *)
Do[F[w] = Collect[F[w], _G, Simplify], {w, 1, wmax}];
words = Select[Union[Cases[Table[F[w], {w, 1, wmax}], _G, Infinity]],
   Function[g, ! AllTrue[Flatten[Table[Coefficient[F[w], g], {w, 1, wmax}]], PossibleZeroQ]]];
Print["letters in the transported functions: ", Union[Flatten[List @@@ words /. q -> Nothing]]];
intR[a_, e_] := Module[{tt = Expand[e]},
   Total[Map[Function[term,
       If[FreeQ[term, G[___, r]],
         term (-G[0, r] + If[a === 0, 0, G[1/a, r]]),
         term /. G[u___, r] :> (-G[0, u, r] + If[a === 0, 0, G[1/a, u, r]])]],
     If[Head[tt] === Plus, List @@ tt, {tt}]]]];
Tinv[{}] = 1;
Tinv[{a_, rest___}] := Tinv[{a, rest}] = Module[{S = intR[a, Tinv[{rest}]]},
    Expand[G[a, rest, 1] - (S /. G[u__, r] :> G[u, 1]) + S]];
Rexact[G[ws__, q]] := Expand[Tinv[{ws}] /. G[__, r] -> 0] /. G[0 .., 1] -> 0;
rules = Table[g -> Rexact[g], {g, words}];
Print["regularised limits of ", Length[words], " G's  (", Round[AbsoluteTime[] - t0, 0.1], " s)"];
rules = rules /. g_G :> ShuffleRegulate[g] /; MatchQ[g, G[1, ___, 1]];
Do[
  dv = Expand[F[w] /. rules];
  dmax = Max[Abs[N[Flatten[Table[Coefficient[dv, HPL[{1}, 1], k], {k, 1, 3}]] /. g_G :> Ginsh[g, {}], 20]] /. Indeterminate -> 99];
  If[dmax > 10^-15, Print["!! weight ", w, ": the divergence Log[1 - 1] does not cancel (", dmax, ")"]; Exit[1]],
  {w, 1, wmax}];
rules = rules /. HPL[{1}, 1] -> 0;

(* ---- constants: simplify, fit what is left on Q[i Pi, MZV] ---- *)
(* real and imaginary parts each on the full weight-w basis of Pi^k x MZV (the odd masters carry an extra I);
   a relation is accepted only with small integers (true ones here are < 10^10; spurious ones at 100 digits
   are ~ 10^20) and if it reproduces the value to 1e-80 *)
basis = <|1 -> {Pi}, 2 -> {Pi^2}, 3 -> {Zeta[3], Pi^3}, 4 -> {Pi^4, Pi Zeta[3]},
          5 -> {Pi^5, Pi^2 Zeta[3], Zeta[5]}, 6 -> {Pi^6, Zeta[3]^2, Pi^3 Zeta[3], Pi Zeta[5]}|>;
fit[xv_, b_] := Which[Abs[xv] < 10^-80, 0, Precision[xv] < 90, $Failed,
   True, Module[{nv = FindIntegerNullVector[Join[{xv}, N[b, 120]], WorkingPrecision -> 95]},
     If[Head[nv] === FindIntegerNullVector || nv[[1]] == 0 || Max[Abs[nv]] > 10^13, $Failed, -Rest[nv] . b/nv[[1]]]]];
other = G | Catalan | _Log | _PolyLog | _HPL | Zeta[_, _] | _PolyGamma;
Do[
  cst = Expand[F[w] /. rules] /. Log[n_Integer /; n > 1] :> Total[#2 Log[#1] & @@@ FactorInteger[n]];
  cst = cst /. {Zeta[5, 1/4] -> 496 Zeta[5] + 5 Pi^5/3, Zeta[5, 3/4] -> 496 Zeta[5] - 5 Pi^5/3,
                PolyGamma[3, 3/4] -> 16 Pi^4 - PolyGamma[3, 1/4]};
  cst = Simplify /@ Expand[cst];
  left = Select[Range[nf], ! FreeQ[cst[[#]], other] &];
  If[left =!= {},
    gs = Union[Cases[cst[[left]], G[__, 1], Infinity]];
    gv = Table[g -> Ginsh[g, {}, PrecisionGoal -> 110], {g, gs}];
    Do[
      xv = N[cst[[k]] /. gv, 100];
      f = fit[Re[xv], basis[w]] + I fit[Im[xv], basis[w]];
      If[! FreeQ[f, $Failed] || Abs[N[f - xv, 100]] > 10^-80, Print["!! weight ", w, " component ", k, " not fitted on Q[i Pi, MZV]"]; Exit[1]];
      cst[[k]] = Expand[f],
      {k, left}]];
  If[! FreeQ[cst, other | q | s | p], Print["!! weight ", w, ": constants left: ", Union[Cases[cst, other, Infinity]]]; Exit[1]];
  Put[cst, FileNameJoin[{inp, "bdy_" <> outTag <> "_w" <> ToString[w] <> ".m"}]];
  Print["weight ", w, ": ", Length[left], " components fitted; constants ", Union[Cases[cst, Pi | _Zeta, Infinity]],
        "; differ from ", seedTag, " in ", Count[Expand[cst - seed[w]], d_ /; d =!= 0], " components  (", Round[AbsoluteTime[] - t0, 0.1], " s)"],
  {w, 1, wmax}];
Print["done (", Round[AbsoluteTime[] - t0, 0.1], " s)"];
