(* v2: singlet. As ../XY_POLES/convertFC.wl but on ../poles/PolePrediction/AmpPS<i>.m (Nc, Nf layers; the
   overall NfZ is NOT in the files and is multiplied in at the combine step), fibration fibrXY + ext + ext2.
   Writes fc_ps_<i>.m. *)
SetEnvironment["PATH" -> Environment["PATH"] <> ":/usr/local/bin"];
Quiet@Get["PolyLogTools`"];
dir = DirectoryName[$InputFileName];
lc = FileNameJoin[{ParentDirectory[ParentDirectory[dir]], "LC_Z_convert_to_II"}];
t0 = AbsoluteTime[];
(* fibrXY_ext.m: rules for the G[..., 1-x-y] of the subleading colour layers missing in fibrXY.m (extendFibr.wl) *)
dispatched = Dispatch@Join[Get[FileNameJoin[{lc, "fibrXY.m"}]], Get[FileNameJoin[{lc, "fibrXY_ext.m"}]], Get[FileNameJoin[{dir, "fibrXY_ext2.m"}]]];
Print["fibrXY loaded (", Round[AbsoluteTime[] - t0], "s)"];

(* Laurent coefficients in (Nc, Nf) without expanding anything else: Association {a,b} -> coef *)
conv[p_, q_] := Merge[Flatten@KeyValueMap[Function[{i, a},
      KeyValueMap[Function[{j, b}, (i + j) -> a b], q]], p], Total];
colL[e_] /; FreeQ[e, Nc | Nf] := <|{0, 0} -> e|>;
colL[Nc] := <|{1, 0} -> 1|>;
colL[Nf] := <|{0, 1} -> 1|>;
colL[Power[Nc, k_Integer]] := <|{k, 0} -> 1|>;
colL[Power[Nf, k_Integer?Positive]] := <|{0, k} -> 1|>;
colL[e_Plus] := Merge[colL /@ List @@ e, Total];
colL[e_Times] := Fold[conv, <|{0, 0} -> 1|>, colL /@ List @@ e];
colL[Power[b_, n_Integer?Positive]] := Nest[conv[#, colL[b]] &, <|{0, 0} -> 1|>, n];
colL[e_] := (Print["WARNING: cannot split in Nc, Nf: ", Short[e]]; Exit[1]);

(* colour split *)
poles = Table[Module[{t = AbsoluteTime[], tmp, cl},
    tmp = Get[FileNameJoin[{ParentDirectory[dir], "PolePrediction", "AmpPS" <> ToString[i] <> ".m"}]];
    cl = KeySort@DeleteCases[colL[tmp], 0];
    Print["AmpPS", i, ": colour layers {Nc^a, Nf^b} = ", Keys[cl], "  (", Round[AbsoluteTime[] - t], "s)"];
    cl], {i, 6}];

(* z -> 1-x-y and fibration, over all G's of all files and colour layers *)
gs = Union @@ (Cases[Values[#], _G, Infinity] & /@ poles);
gz = Select[gs, ! FreeQ[#, z] &];
mapgs = Dispatch@Thread[(gz /. z -> 1 - x - y) -> (gz /. z -> 1 - x - y /. dispatched)];
notFibr = Select[gz /. z -> 1 - x - y, (# /. dispatched) === # &];
Print[Length[gs], " distinct G (", Length[gz], " with z); not in fibrXY: ", Length[notFibr],
  If[notFibr === {}, "", {"  e.g. ", Take[notFibr, UpTo[5]]}]];
poles = Map[(# /. z -> 1 - x - y /. mapgs) &, poles, {2}];
If[! FreeQ[poles, z], Print["WARNING: z left after substitution"]];

(* G -> II *)
gs = Union @@ (Cases[Values[#], _G, Infinity] & /@ poles);
Print[Length[gs], " distinct G after fibration; computing symbols"];
t = AbsoluteTime[];
symbGSRules = Dispatch@Thread[gs -> (ComputeSymbol[gs] /. CiTi -> II)];
Print["   symbols done (", Round[AbsoluteTime[] - t], "s)"];

Do[Module[{res = Map[# /. symbGSRules &, poles[[i]]]},
   If[! FreeQ[res, G], Print["WARNING: G left in file ", i]];
   Put[Plus @@ KeyValueMap[Nc^#1[[1]] Nf^#1[[2]] #2 &, res],
    FileNameJoin[{dir, "fc_ps_" <> ToString[i] <> ".m"}]];
   Print["wrote fc_ps_", i, ".m"]], {i, 6}];
Print["total ", Round[AbsoluteTime[] - t0], "s"];
