(* compare_exact.wls -- exact comparison of two families of constant files after reducing the identities that old files
   may leave unreduced: PolyGamma[3,1/4] + PolyGamma[3,3/4] = 16 Pi^4, Hurwitz Zeta[5, 1/4], Zeta[5, 3/4].
   usage: wolframscript -file compare_exact.wls MINE_PATTERN REF_PATTERN W1 W2   (patterns with {w}) *)
{mine, ref, w1, w2} = Rest[$ScriptCommandLine]; {w1, w2} = ToExpression /@ {w1, w2};
red = {PolyGamma[3, 3/4] -> 16 Pi^4 - PolyGamma[3, 1/4], Zeta[5, 1/4] -> 496 Zeta[5] + 5 Pi^5/3, Zeta[5, 3/4] -> 496 Zeta[5] - 5 Pi^5/3};
ok = True;
Do[
  a = Get[StringReplace[mine, "{w}" -> ToString[w]]]; b = Get[StringReplace[ref, "{w}" -> ToString[w]]];
  d = Select[Range[Length[a]], Simplify[Expand[(a[[#]] - b[[#]]) /. red]] =!= 0 &];
  If[d =!= {}, ok = False];
  Print["weight ", w, ": ", Length[d], " of ", Length[a], " components differ", If[d =!= {}, ", e.g. " <> ToString[Take[d, UpTo[5]]], ""]],
  {w, w1, w2}];
Print[If[ok, "IDENTICAL", "DIFFERENT"]];
