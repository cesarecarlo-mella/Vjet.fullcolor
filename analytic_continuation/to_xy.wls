(* to_xy.wls -- LC production-region finite remainders (2024_Vjet_3L_LC_ancillary/FiniteRemainders,
   region y = s13/q2 > 1, variables v = 1/y, u = -z/y; G[..., u] with letters 0, 1, 1-v, -v and
   G[..., v] with letters 0, 1) -> input format of xyseries.py, with (u, v) renamed (x, y).

   Every line is  + (rational coefficient in x, y)*constants*G[..., x]*G[..., y].
   Constants: Pi, Zeta[3], Zeta[5], Nc, Nf, and ii = I; q2 -> 1 (overall factor).

   usage: wolframscript -file to_xy.wls <outdir> <channel> [<channel> ...]     e.g. lrP lmL *)
Get["PolyLogTools`"];                     (* always first, before anything that contains G *)
dir = FileNameJoin[{DirectoryName[$InputFileName], "..", "reproduceLC", "2024_Vjet_3L_LC_ancillary", "FiniteRemainders"}];
{out, channels} = {$ScriptCommandLine[[2]], Drop[$ScriptCommandLine, 2]};
str[e_] := ToString[e, InputForm, PageWidth -> Infinity];
Do[
  rem = Get[FileNameJoin[{dir, ch <> ".m"}]];
  Do[
    key = rule[[1]]; name = ch <> "_" <> ToString[Head[key]] <> ToString[key[[1]]];
    (* numeric constants become plain symbols, so that CoefficientRules accepts them *)
    e = rule[[2]] /. Complex[a_, b_] :> a + b ii /. {u -> x, v -> y} /.
        {Pi -> pp, Zeta[3] -> z3, Zeta[5] -> z5, q2 -> 1};
    gs = Union[Cases[e, _G, Infinity]];
    vars = Join[Array[gg, Length[gs]], {pp, z3, z5, Nc, Nf, ii}];
    cr = CoefficientRules[e /. Thread[gs -> Array[gg, Length[gs]]], vars];
    back = Join[Thread[Array[gg, Length[gs]] -> gs], {pp -> Pi, z3 -> Zeta[3], z5 -> Zeta[5]}];
    lines = Table[
      Module[{facs = Flatten[MapThread[Table[#1, #2] &, {vars, r[[1]]}]] /. back},
        If[Count[facs, G[__, x]] > 1 || Count[facs, G[__, y]] > 1,
          Print["!! two G in the same variable in ", name]; Exit[1]];
        " + (" <> str[Together[r[[2]]]] <> ")" <> StringJoin[("*" <> str[#]) & /@ facs]],
      {r, cr}];
    Export[FileNameJoin[{out, name <> "_G.m"}], name <> " = (0\n" <> StringRiffle[lines, "\n"] <> "\n);\n", "Text"];
    Print[name, ": ", Length[cr], " terms, ", Length[gs], " distinct G"],
    {rule, rem}],
  {ch, channels}];
