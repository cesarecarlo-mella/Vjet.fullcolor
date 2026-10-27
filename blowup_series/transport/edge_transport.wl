(* edge_transport.wl -- boundary constants of the s = 0 charts R<A>_<b>v of a continued region (one family):
   transport along the edge bh = 0 of the region's triangle, from the threshold corner to the high-energy corner.

   Region A > 1 (charts.py):  A = 1/v,  b = -bh/v for the other two variables,  v + bh + ch = 1.
   The edge bh = 0 runs from the threshold corner (ch = 0, v = 1: the decay corner continued, constants
   bdy_R<A>_th_w<w>.m from threshold.py) to the high-energy corner (ch = 1, bh = v = 0), reached with bh -> 0
   first: the ordering bh << v of chart R<A>_<b>v (bh = t v^2).  (For A = z this is transport.boundaries:
   R<Z>_yv = CC2 along u = yh = 0, R<Z>_xv = CC3 along w = xh = 0.)

   On the edge (ch = s, bh = eps -> 0) every even letter becomes eps^k g_k(s) with g_k in {s, 1 - s}; the
   odd letters tend to +-1 (xyz -> 0), so their forms vanish.  The eps^0 part obeys
       dF/ds = A(s) F ,   A = Sum_k atilde_k d/ds Log[g_k]       (weight by weight, eps absorbed)
   integrated with GIntegrate from the tangential base point s = 0 (seed: the threshold constants), then
   s -> 1 - t, change of fibration (fibrtto1.m), drop G[0,..,0,t], t -> 0.

   usage: wolframscript -file edge_transport.wl <family inputs dir> <region x|y|z> [WMAX = 6]
   writes <inputs>/bdy_R<A>_<b>v_w<w>.m for both edges (w = 0..WMAX) *)
SetEnvironment["PATH" -> GetEnvironment["PATH"][[2]] <> ":/usr/local/bin:/opt/homebrew/bin"];
Get["PolyLogTools`"];                     (* always first *)
argv = Rest[$ScriptCommandLine];
inp = argv[[1]]; region = argv[[2]]; wmax = If[Length[argv] > 2, ToExpression[argv[[3]]], 6];
here = DirectoryName[$InputFileName];
t0 = AbsoluteTime[];
alph = Get[FileNameJoin[{here, "alphabet.wl"}]];          (* copy of data/trident_definition/alphabet.wl *)
W = Association[Table[k -> (l[k] /. alph /. Log[a_] :> a), {k, 20}]];
atilde = Get[FileNameJoin[{inp, "atilde.m"}]];
nf = Length[atilde];
fibr = Get[FileNameJoin[{here, "fibrtto1.m"}]];
zeroRules = Table[G @@ Join[ConstantArray[0, n], {t}] -> 0, {n, 1, 8}];
seed[w_] := Get[FileNameJoin[{inp, "bdy_R" <> ToUpperCase[region] <> "_th_w" <> ToString[w] <> ".m"}]];
Print["family dir ", inp, ", region ", region, " > 1, ", nf, " masters, weights <= ", wmax];

(* region map: physical (x, y) in terms of v and the two hatted variables *)
vars = {"x", "y", "z"};
others = DeleteCases[vars, region];
hat = <|"x" -> xh, "y" -> yh, "z" -> zh|>;
val = Association[region -> 1/vv, Table[b -> -hat[b]/vv, {b, others}]];
vsub = vv -> 1 - Total[hat /@ others];
X = val["x"] /. vsub; Y = val["y"] /. vsub;
If[Simplify[(1 - X - Y) - (val["z"] /. vsub)] =!= 0, Print["!! region map inconsistent"]; Exit[1]];

lowest[e_, p_] := Module[{n = Numerator[e], d = Denominator[e]},
   {Exponent[n, p, Min] - Exponent[d, p, Min], Together[Coefficient[n, p, Exponent[n, p, Min]]/Coefficient[d, p, Exponent[d, p, Min]]]}];

Do[
  c = First[DeleteCases[others, b]];
  tag = "R" <> ToUpperCase[region] <> "_" <> b <> "v";
  sub = {hat[b] -> eps, hat[c] -> s};
  forms = Association[Table[
     k -> If[MemberQ[{7, 8}, k],
       Module[{a = If[k == 7, x (1 - x - y), x y], r2},
         r2 = lowest[Together[(x y (1 - x - y)/a^2 /. {x -> X, y -> Y}) /. sub], eps];
         If[r2[[1]] == 0, Print["!! odd letter l", k, " not constant on the edge of ", tag]; Exit[1]];
         0],
       Module[{g = lowest[Together[(W[k] /. {x -> X, y -> Y}) /. sub], eps][[2]], f},
         f = Apart[Together[D[Log[g], s]], s];
         If[! PolynomialQ[Together[f s (1 - s)], s] || Exponent[Together[f s (1 - s)], s] > 1,
           Print["!! letter l", k, " has a pole off {0, 1} on the edge of ", tag, ": ", f]; Exit[1]];
         f]],
     {k, 20}]];
  Print[tag, ": edge ", hat[b], " = 0, s = ", hat[c], "; forms ",
        Table[k -> forms[k], {k, Select[Range[20], forms[#] =!= 0 &]}]];
  conn = atilde /. l[k_] :> forms[k];
  int = seed[0];
  Put[int, FileNameJoin[{inp, "bdy_" <> tag <> "_w0.m"}]];
  Do[
    toint = Expand[conn . int];
    int = Map[GIntegrate[#, s] &, toint] + seed[w];
    after = Expand[(int /. s -> 1 - t) /. fibr];
    bad = Select[Union[Cases[after, _G, Infinity]], Last[List @@ #] =!= t &];
    If[bad =!= {}, Print["!! ", tag, " weight ", w, ": G not refibrated: ", Take[bad, UpTo[5]]]; Exit[1]];
    cst = Expand[after /. zeroRules /. t -> 0];
    If[! FreeQ[cst, G | t | s], Print["!! ", tag, " weight ", w, ": not a constant"]; Exit[1]];
    Put[cst, FileNameJoin[{inp, "bdy_" <> tag <> "_w" <> ToString[w] <> ".m"}]];
    Print["  ", tag, " weight ", w, ": nonzero ", Count[cst, q_ /; q =!= 0], "/", nf, ", constants ",
          Union[Cases[cst, Pi | _Zeta | Catalan | _Log | _PolyLog, Infinity]], "  (", Round[AbsoluteTime[] - t0, 0.1], " s)"],
    {w, 1, wmax}],
  {b, others}];
Print["done (", Round[AbsoluteTime[] - t0, 0.1], " s)"];
