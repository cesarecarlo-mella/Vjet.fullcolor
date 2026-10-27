(* Checks for polesII_polesonly.wl / polesII_eps0.wl output, as checkPoles.wl but on one part of the eps expansion:
   usage: checkPolesOnly.wl fc_poles_<i>.m [poles|eps0]   (default poles)
   (1) pole part (eps^k, k<0) of the original vs <name>_poles_collected, every II -> memoised random rational
       (eps kept symbolic in the original, exact Laurent truncation)
   (2) _poles_collected vs _poles_shuffled, II[a1..an] -> r[a1]..r[an]/n! (shuffle homomorphism)
   Exact rational arithmetic: deviations are exact. *)
base = FileBaseName[$ScriptCommandLine[[2]]];
part = If[Length[$ScriptCommandLine] > 2, $ScriptCommandLine[[3]], "poles"];
keep = If[part === "eps0", # == 0 &, # < 0 &];
dir = DirectoryName[$InputFileName];
f[s_] := FileNameJoin[{dir, If[s === "", Nothing, "FC_POLES"], base <> s <> ".m"}];  (* original here, results in FC_POLES/ *)
pts = {{1371/10000, 2913/10000, 3, 5, 1/7, 1}, {2/7, 3/11, 4, 2, -3/13, 2}, {41/100, 13/100, 5, 3, 2/9, 3}};
rnd[] := RandomInteger[{-1000, 1000}]/RandomInteger[{1, 1000}];
rel[a_, b_] := If[Expand[a - b] === 0, "0 (exact)", ToString@NumberForm[N[Abs[a - b]/Abs[a]], 2]];
Do[
  {x, y, Nc, Nf} = p[[1 ;; 4]]; e0 = p[[5]];
  SeedRandom[p[[6]]]; ClearAll[II]; II[args__] := II[args] = rnd[];
  Clear[eps]; orig = Expand@Get[f[""]];
  poles = Select[If[Head[orig] === Plus, List @@ orig, {orig}], keep[Exponent[#, eps]] &];
  If[! FreeQ[Exponent[#, eps] & /@ poles, -Infinity], Print["WARNING: eps not monomial"]];
  v0 = Total[poles] /. eps -> e0;
  eps = e0; v1 = Get[f["_" <> part <> "_collected"]];
  SeedRandom[p[[6]]]; ClearAll[II, r]; r[l_] := r[l] = rnd[];
  II[w__] := Times @@ (r /@ {w})/Length[{w}]!;
  w1 = Get[f["_" <> part <> "_collected"]]; w2 = Get[f["_" <> part <> "_shuffled"]];
  Clear[eps];
  Print[base, "  point ", p[[1 ;; 5]], ":  orig ", part, " part vs collected ", rel[v0, v1],
    ",  collected vs shuffled ", rel[w1, w2]],
  {p, pts}];
