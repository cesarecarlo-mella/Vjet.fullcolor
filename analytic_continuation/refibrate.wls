(* refibrate.wls -- one LC remainder coefficient (one colour layer), re-expressed for the expansion
   at another corner of the production triangle (v = 1/y, u = -z/y, w = 1 - u - v), and written in
   the xyseries.py input format.

     corner th   (threshold, v -> 1):              v = 1 - u - w, fibration {w, u}, (w, u) -> (x, y)
     corner u1   (u -> 1, chart around s~ = 0):    u = 1 - v - w, fibration {w, v}, (w, v) -> (x, y)
     corner he2b (w -> 1, chart around s = 1):     fibration {v, u}, (v, u) -> (x, y)
     corner he3b (u -> 1, chart around s~ = 1):    u = 1 - v - w, fibration {v, w}, (v, w) -> (x, y)

   In both cases the G's in the first variable have letters {0, 1, 1 - y, -y} and those in the second
   {0, 1}, exactly the structure xyseries.py expands (first x at fixed y, then y).

   usage: wolframscript -file refibrate.wls <channel> <coef, e.g. alpha1> <layer, e.g. Nc2Nf1> <corner> <out.m> *)
Get["PolyLogTools`"];                     (* always first, before anything that contains G *)
{ch, coef, layer, corner, out} = $ScriptCommandLine[[2 ;; 6]];
dir = FileNameJoin[{DirectoryName[$InputFileName], "..", "reproduceLC", "2024_Vjet_3L_LC_ancillary", "FiniteRemainders"}];
rem = Get[FileNameJoin[{dir, ch <> ".m"}]];
{hd, idx} = StringCases[coef, (h : LetterCharacter ..) ~~ (n : DigitCharacter ..) :> {h, ToExpression[n]}][[1]];
e = (ToExpression[hd][idx, 3] /. rem) /. q2 -> 1;
{nc, nf} = ToExpression /@ StringCases[layer, "Nc" ~~ a : DigitCharacter .. ~~ "Nf" ~~ b : DigitCharacter .. :> {a, b}][[1]];
e = Coefficient[Coefficient[Expand[e, Nc | Nf], Nc, nc], Nf, nf];
t0 = AbsoluteTime[];
{sub, fib, ren} = Switch[corner,
   "th", {v -> 1 - u - w, {w, u}, {w -> x, u -> y}},          (* threshold corner, u, w -> 0 *)
   "u1", {u -> 1 - v - w, {w, v}, {w -> x, v -> y}},          (* corner u -> 1, chart around s~ = 0 *)
   "he2b", {{}, {v, u}, {v -> x, u -> y}},                    (* corner w -> 1, chart around s = 1 *)
   "he3b", {u -> 1 - v - w, {v, w}, {v -> x, w -> y}}];       (* corner u -> 1, chart around s~ = 1 *)
mx = StringReplace[out, ".m" ~~ EndOfString -> ".refib.m"];
If[FileExistsQ[mx],
  f = Get[mx]; Print["   refibrated expression read from ", mx],
  f = ToFibrationBasis[e /. sub, fib];
  (* G's whose argument is not a fibration variable (e.g. G[.., 1 - v - w]) occasionally survive one
     pass: refibrate those again *)
  Do[
    left = Union[Cases[f, G[___, a_] /; ! MemberQ[fib, a] && ! NumericQ[a], Infinity]];
    If[left === {}, Break[]];
    Print["   pass ", pass, ": ", Length[left], " G's with a compound argument, refibrating them again: ", InputForm[Take[left, UpTo[3]]], " -> ", InputForm[ToFibrationBasis[First[left], fib]]];
    f = f /. Thread[left -> (ToFibrationBasis[#, fib] & /@ left)],
    {pass, 1, 3}];
  If[Cases[f, G[___, a_] /; ! MemberQ[fib, a] && ! NumericQ[a], Infinity] =!= {}, Print["!! compound arguments left"]; Exit[1]];
  f = ShuffleG[f];                      (* products of G with the same argument -> single G's *)
  f = f /. G[z : 0 .., c_?NumericQ] :> Log[c]^Length[{z}]/Length[{z}]!;     (* e.g. G[0, 2] = Log[2] *)
  If[Cases[f, G[___, _?NumericQ], Infinity] =!= {}, Print["!! constant G's left: ", Cases[f, G[___, _?NumericQ], Infinity]]; Exit[1]];
  f = Collect[f, {_G}, Together];
  Put[f, mx];                           (* keep it: the refibration is the expensive step *)
  Print[ch, " ", coef, " ", layer, " corner ", corner, ": refibrated in ", Round[AbsoluteTime[] - t0, 0.1], " s, ",
        Length[Union[Cases[f, _G, Infinity]]], " distinct G  -> ", mx]];
(* write in the xyseries input format *)
str[x_] := ToString[x, InputForm, PageWidth -> Infinity];
g = f /. ren /. Complex[a_, b_] :> a + b ii /. {Pi -> pp, Zeta[3] -> z3, Zeta[5] -> z5, Log[2] -> lg2};
gs = Union[Cases[g, _G, Infinity]];
vars = Join[Array[gg, Length[gs]], {pp, z3, z5, ii, lg2}];
cr = CoefficientRules[g /. Thread[gs -> Array[gg, Length[gs]]], vars];
back = Join[Thread[Array[gg, Length[gs]] -> gs], {pp -> Pi, z3 -> Zeta[3], z5 -> Zeta[5], lg2 -> Log[2]}];
bad = Select[gs, ! SubsetQ[If[Last[#] === x, {0, 1, 1 - y, -y}, {0, 1}], Most[List @@ #]] &];
If[bad =!= {}, Print["!! letters outside {0, 1, 1-y, -y} in x / {0, 1} in y: ", Take[bad, UpTo[5]]]; Exit[1]];
lines = Table[
   Module[{facs = Flatten[MapThread[Table[#1, #2] &, {vars, r[[1]]}]] /. back},
     If[Count[facs, G[__, x]] > 1 || Count[facs, G[__, y]] > 1, Print["!! two G in one variable"]; Exit[1]];
     " + (" <> str[Together[r[[2]]]] <> ")" <> StringJoin[("*" <> str[#]) & /@ facs]],
   {r, cr}];
name = ch <> "_" <> coef <> "_" <> layer <> "_" <> corner;
Export[out, name <> " = (0\n" <> StringRiffle[lines, "\n"] <> "\n);\n", "Text"];
Print["   wrote ", out, ": ", Length[cr], " terms, ", Length[gs], " distinct G  (", Round[AbsoluteTime[] - t0, 0.1], " s)"];
