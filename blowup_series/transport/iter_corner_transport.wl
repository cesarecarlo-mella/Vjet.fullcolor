(* iter_corner_transport.wl -- exact boundary constants of the iterated-integral system (family ITER) at the decay
   corners c3 and c2, transported from c1 (where every constant of weight >= 1 vanishes: the shuffle-regularised base
   point) along the two edges:
     edge x = 0, s = y from c1 to c3 (x, z -> 0): the edge is the t = 0 line of chart 3 (x = t z^2)   -> bdy_c3_w<w>.m
     edge y = 0, s = x from c1 to c2 (y, z -> 0): the edge is the t = 0 line of chart 5 (y = t z^2)   -> bdy_c2_w<w>.m
     (the pipeline reads the constants of the point c2 for both of its charts, 2 and the mirror 5)
   On an edge the system is  dF/ds = (A0/s + A1/(s - 1)) . F  (edge_x0.m, edge_y0.m from family_inputs.py: every letter
   in use is c s^a (1 - s)^b there; x, resp. y, and the odd letters drop out).
   GIntegrate weight by weight from the regularised base point s = 0 (weight 0: sol0, all higher constants 0), then
   s -> 1 - t, change of fibration (fibrtto1.m), drop G[0,..,0,t] (the powers of log t = log z, regularised to 0 in the
   chart, like log t of the chart), t -> 0.  The constants are MZVs.

   usage: wolframscript -file iter_corner_transport.wl <ITER inputs dir> [WMAX = 6] *)
SetEnvironment["PATH" -> GetEnvironment["PATH"][[2]] <> ":/usr/local/bin:/opt/homebrew/bin"];
Get["PolyLogTools`"];                     (* always first *)
argv = Rest[$ScriptCommandLine];
inp = argv[[1]]; wmax = If[Length[argv] > 1, ToExpression[argv[[2]]], 6];
here = DirectoryName[$InputFileName];
t0 = AbsoluteTime[];
fibr = Get[FileNameJoin[{here, "fibrtto1.m"}]];
zeroRules = Table[G @@ Join[ConstantArray[0, n], {t}] -> 0, {n, 1, 8}];
sol0 = Get[FileNameJoin[{inp, "sol0.m"}]];
nf = Length[sol0];
Print["ITER inputs ", inp, ": ", nf, " functions, weights <= ", wmax];

Do[
  {edge, tag} = et;
  {A0, A1} = Get[FileNameJoin[{inp, edge <> ".m"}]];
  If[Dimensions[A0] =!= {nf, nf} || Dimensions[A1] =!= {nf, nf}, Print["!! ", edge, ": wrong dimensions"]; Exit[1]];
  int = sol0;
  Do[
    toint = Expand[(A0 . int)/s + (A1 . int)/(s - 1)];
    int = Map[GIntegrate[#, s] &, toint];
    after = Expand[(int /. s -> 1 - t) /. fibr];
    bad = Select[Union[Cases[after, _G, Infinity]], Last[List @@ #] =!= t &];
    If[bad =!= {}, Print["!! ", tag, " weight ", w, ": G not refibrated: ", Take[bad, UpTo[5]]]; Exit[1]];
    cst = Expand[after /. zeroRules /. t -> 0];
    If[! FreeQ[cst, G | t | s], Print["!! ", tag, " weight ", w, ": not a constant"]; Exit[1]];
    Put[cst, FileNameJoin[{inp, "bdy_" <> tag <> "_w" <> ToString[w] <> ".m"}]];
    Print["  ", tag, " (", edge, ") weight ", w, ": nonzero ", Count[cst, q_ /; q =!= 0], "/", nf, ", constants ",
          Union[Cases[cst, Pi | _Zeta | Catalan | _Log | _PolyLog, Infinity]], "  (", Round[AbsoluteTime[] - t0, 0.1], " s)"],
    {w, 1, wmax}],
  {et, {{"edge_x0", "c3"}, {"edge_y0", "c2"}}}];
Print["done (", Round[AbsoluteTime[] - t0, 0.1], " s)"];
