(* EPS^0 ONLY: as polesII.wl but keeps only the eps^0 pieces (all colour layers);
   writes FC_POLES/<name>_eps0_collected.m and FC_POLES/<name>_eps0_shuffled.m.
   For one file (single expression):
   split by color Nc^a Nf^b -> by power eps^k -> collect II monomials (coefficients
   Collect[.., {Pi,_Zeta}, Together]) -> linearise II products with the shuffle product.
   Writes <name>_collected.m and <name>_shuffled.m. Single kernel. *)
file = $ScriptCommandLine[[2]];
dir = DirectoryName[$InputFileName];
base = FileBaseName[file] <> "_eps0";
t0 = AbsoluteTime[];
expr = Get[FileNameJoin[{dir, file}]];
Print[base, ": loaded in ", Round[AbsoluteTime[] - t0], "s"];

topTerms[e_] := If[Head[e] === Plus, List @@ e, {e}];
factorSplit[term_, patt_] := Module[{f = If[Head[term] === Times, List @@ term, {term}], sel},
  sel = FreeQ[#, patt] & /@ f;
  {Times @@ Pick[f, sel, False], Times @@ Pick[f, sel, True]}];

(* Laurent coefficients in eps without expanding anything else: Association k -> coef *)
conv[p_, q_] := Merge[Flatten@KeyValueMap[Function[{i, a},
      KeyValueMap[Function[{j, b}, (i + j) -> a b], q]], p], Total];
epsL[e_] /; FreeQ[e, eps] := <|0 -> e|>;
epsL[eps] := <|1 -> 1|>;
epsL[Power[eps, k_Integer]] := <|k -> 1|>;
epsL[e_Plus] := Merge[epsL /@ List @@ e, Total];
epsL[e_Times] := Fold[conv, <|0 -> 1|>, epsL /@ List @@ e];
epsL[Power[b_, n_Integer?Positive]] := Nest[conv[#, epsL[b]] &, <|0 -> 1|>, n];
epsL[e_] := (Print["WARNING: cannot split in eps: ", Short[e]]; Abort[]);

(* collect in II *)
allII = DeleteDuplicates@Cases[expr, _II, Infinity];
(* II <-> iiSym via hashed Association lookup: e /. Dispatch[...] with ~40k II rules is not hashed
   and cost ~10 min per piece on the full-colour files (polesII_orig.wl); same index assignment *)
iiIdx = Association@MapIndexed[#1 -> iiSym[First@#2] &, allII];
simpCoef[c_] := Collect[c, {Pi, _Zeta}, Together];
collectII[e_] := Module[{x1, grp},
  x1 = Expand[e /. ii_II :> iiIdx[ii], iiSym[_]];
  grp = GroupBy[factorSplit[#, iiSym] & /@ topTerms[x1], First -> Last, Total];
  Plus @@ KeyValueMap[(#1 /. iiSym[n_] :> allII[[n]]) simpCoef[#2] &, grp]];

(* shuffle *)
sh[{}, w_List] := {w};
sh[w_List, {}] := {w};
sh[u : {a_, as___}, v : {b_, bs___}] := sh[u, v] =
   Join[Prepend[#, a] & /@ sh[{as}, v], Prepend[#, b] & /@ sh[u, {bs}]];
shLin[p_Association, q_Association] := Merge[Flatten@KeyValueMap[
     Function[{u, cu}, KeyValueMap[Function[{v, cv},
        (# -> cu cv) & /@ sh[u, v]], q]], p], Total];
monoWords[m_] := monoWords[m] = Module[{f = If[Head[m] === Times, List @@ m, {m}]},
  f = Cases[Flatten[f /. Power[ii_II, k_Integer?Positive] :> ConstantArray[ii, k]], _II];
  Fold[shLin, <|{} -> 1|>, <|List @@ # -> 1|> & /@ f]];
If[! (sh[{a}, {b}] === {{a, b}, {b, a}} && Length[sh[{1, 2, 3}, {4, 5}]] == 10 &&
     monoWords[II[a]^2] === <|{a, a} -> 2|> && monoWords[1] === <|{} -> 1|>),
  Print["shuffle self-test FAILED"]; Exit[1]];
shuffle[e_] := Module[{acc},
  acc = GroupBy[Flatten[Function[{m, c}, KeyValueMap[(#1 -> #2 c) &, monoWords[m]]] @@@
       (factorSplit[#, II] & /@ topTerms[e])], First -> Last];
  acc = If[Length[#] == 1, First[#], simpCoef[Total[#]]] & /@ acc;
  Plus @@ KeyValueMap[If[#1 === {}, 1, II @@ #1] #2 &, acc]];

(* pieces: {color, k, expr} *)
colors = Normal@GroupBy[factorSplit[#, Nc | Nf] & /@ topTerms[expr], First -> Last, Total];
If[! FreeQ[colors[[All, 2]], Nc | Nf], Print["WARNING: Nc/Nf inside a color coefficient"]];
Print["   colors: ", colors[[All, 1]]];
pieces = Flatten[Function[{col, ce}, KeyValueMap[{col, #1, #2} &, KeySort@epsL[ce]]] @@@ colors, 1];
pieces = Select[pieces, #[[2]] == 0 &];
Print["   (color, eps power) pieces: ", pieces[[All, 1 ;; 2]]];

coll = Table[Module[{t = AbsoluteTime[], r},
    r = collectII[p[[3]]];
    Print["   collected ", p[[1]], " eps^", p[[2]], "  (", Round[AbsoluteTime[] - t], "s)"];
    {p[[1]], p[[2]], r}], {p, pieces}];
Put[Plus @@ (#[[1]] eps^#[[2]] #[[3]] & /@ coll), FileNameJoin[{dir, "FC_POLES", base <> "_collected.m"}]];
Print["   wrote ", base, "_collected.m  (", Round[AbsoluteTime[] - t0], "s)"];

shuf = Table[Module[{t = AbsoluteTime[], r},
    r = shuffle[p[[3]]];
    Print["   shuffled ", p[[1]], " eps^", p[[2]], "  (", Round[AbsoluteTime[] - t], "s)"];
    {p[[1]], p[[2]], r}], {p, coll}];
Print["   nonlinear II left: ", ! FreeQ[shuf, II[__]^_ | HoldPattern[Times[___, _II, ___, _II, ___]]]];
Put[Plus @@ (#[[1]] eps^#[[2]] #[[3]] & /@ shuf), FileNameJoin[{dir, "FC_POLES", base <> "_shuffled.m"}]];
Print["   wrote ", base, "_shuffled.m  (total ", Round[AbsoluteTime[] - t0], "s)"];
