(* finalise_cc2e.wls -- writes results/bdy_CC2e_H_w<w>.m (w = 1..6): the transported CC2e constants of H
   with every remaining G(...; 1) (weights 5, 6) replaced by its MZV form.  The MZV form is found by an
   integer-relation fit of the real and imaginary part (GiNaC values at ~75 digits) on
     w5: Re {Pi^2 Zeta[3], Zeta[5]}, Im {Pi^5};   w6: Re {Pi^6, Zeta[3]^2}, Im {Pi^3 Zeta[3], Pi Zeta[5]}
   and accepted only if it reproduces the value to 1e-60. *)
SetEnvironment["PATH" -> GetEnvironment["PATH"][[2]] <> ":/usr/local/bin:/opt/homebrew/bin"];
Get["PolyLogTools`"];
dir = FileNameJoin[{DirectoryName[$InputFileName], "results"}] <> "/";
basis = <|5 -> {{Pi^2 Zeta[3], Zeta[5]}, {Pi^5}}, 6 -> {{Pi^6, Zeta[3]^2}, {Pi^3 Zeta[3], Pi Zeta[5]}}|>;
fit[x_, b_] := Module[{nv = FindIntegerNullVector[Join[{x}, N[b, 100]], WorkingPrecision -> Floor[Precision[x]] - 3]},
   If[Head[nv] === FindIntegerNullVector || nv[[1]] == 0, $Failed, -Rest[nv] . b/nv[[1]]]];
If[Abs[N[Zeta[5, 1/4] - 496 Zeta[5] - 5 Pi^5/3, 40]] > 10^-35 || Abs[N[Zeta[5, 3/4] - 496 Zeta[5] + 5 Pi^5/3, 40]] > 10^-35,
  Print["!! Hurwitz identity"]; Exit[1]];
Do[
  c = Get[dir <> "bdy_CC2e_H_w" <> ToString[w] <> "_transported.m"];
  withG = Select[Range[Length[c]], ! FreeQ[c[[#]], G] &];
  If[withG =!= {},
    gs = Union[Cases[c[[withG]], G[__, 1], Infinity]];
    gv = Table[g -> Ginsh[g, {}, PrecisionGoal -> 70], {g, gs}];
    Do[
      x = N[c[[k]] /. gv, 80];
      f = fit[Re[x], basis[w][[1]]] + I fit[Im[x], basis[w][[2]]];
      If[! FreeQ[f, $Failed] || Abs[N[f - x, 90]] > 10^-60, Print["!! component ", k, " not fitted"]; Exit[1]];
      c[[k]] = Expand[f],
      {k, withG}]];
  If[! FreeQ[c, G], Print["!! G left at weight ", w]; Exit[1]];
  (* Hurwitz zetas and logs of numbers from PolyLogTools' own reductions: back to Zeta[5], Pi, Log[2] *)
  c = c /. {Zeta[5, 1/4] -> 496 Zeta[5] + 5 Pi^5/3, Zeta[5, 3/4] -> 496 Zeta[5] - 5 Pi^5/3} /.
        Log[n_?NumericQ] :> PowerExpand[Log[Factor[n]]];
  c = Simplify /@ Expand[c];
  Put[c, dir <> "bdy_CC2e_H_w" <> ToString[w] <> ".m"];
  Print["weight ", w, ": ", Length[withG], " components fitted; constants: ", Union[Cases[c, Pi | _Zeta | Catalan | _Log, Infinity]],
        "  -> bdy_CC2e_H_w", w, ".m"],
  {w, 1, 6}];
