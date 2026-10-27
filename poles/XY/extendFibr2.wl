(* Extend fibrXY.m to the G[..., 1-x-y] of the full-colour poles (poles_delta/Amp*.m, z -> 1-x-y)
   (v2: source = ../poles/PolePrediction/AmpPS<i>.m, singlet; known = fibrXY.m + fibrXY_ext.m) that are not yet in it, closing the set under d/dx, d/dy and working up in weight.
   Same procedure as the commented cells of conversion.nb:
     tofix = Int dx (DG[g,x] /. fibr);  toadd = Int dy (DG[g,y] - DG[ShuffleG[tofix],y] /. fibr);
     constant c = g - tofix - toadd fixed by Ginsh at 50 digits + PSLQ on the weight-w MZVs,
     then checked at a second point.
   Writes fibrXY_ext2.m (new rules only). *)
SetEnvironment["PATH" -> Environment["PATH"] <> ":/usr/local/bin"];
Quiet@Get["PolyLogTools`"];
dir = DirectoryName[$InputFileName];  (* poles/XY *)
lc = FileNameJoin[{ParentDirectory[ParentDirectory[dir]], "LC_Z_convert_to_II"}];
t0 = AbsoluteTime[];
fibr = Join[Get[FileNameJoin[{lc, "fibrXY.m"}]], Get[FileNameJoin[{lc, "fibrXY_ext.m"}]]];
known = Association[# -> True & /@ fibr[[All, 1]]];
disp = Dispatch[fibr];
Print["fibrXY: ", Length[fibr], " rules (", Round[AbsoluteTime[] - t0], "s)"];

gz = Union @@ Table[Cases[Get[FileNameJoin[{ParentDirectory[dir], "PolePrediction", "AmpPS" <> ToString[i] <> ".m"}]],
     G[__, z], Infinity], {i, 6}];
needQ[g_] := ! FreeQ[Last[g], 1 - x - y] && ! KeyExistsQ[known, g];
todo = Select[gz /. z -> 1 - x - y, needQ];
(* close under derivatives *)
new = todo;
While[new =!= {},
  new = Complement[Select[Union@Cases[Expand[{DG[#, x], DG[#, y]} & /@ new], _G, Infinity], needQ], todo];
  todo = Union[todo, new]];
todo = SortBy[todo, {Length, Identity}];
Print[Length[todo], " G's to fibrate; weights ", Counts[Length[#] - 1 & /@ todo]];

consts[w_] := Switch[w, 1, {}, 2, {Pi^2}, 3, {Zeta[3]}, 4, {Pi^4}, 5, {Zeta[5], Pi^2 Zeta[3]},
   6, {Pi^6, Zeta[3]^2}, 7, {Zeta[7], Pi^2 Zeta[5], Pi^4 Zeta[3]}];
pt1 = {x -> 3/17, y -> Zeta[7]/9};
pt2 = {x -> 3/21, y -> Zeta[5]/9};
gn[e_, pt_] := Ginsh[e, pt, PrecisionGoal -> 50];

fibrate[g_] := Module[{w = Length[g] - 1, tofix, dy, cand, d1, cb, nv, c, chk},
  tofix = GIntegrate[#, x] &@Expand[DG[g, x] /. disp];
  dy = Expand@Together[DG[g, y] - DG[ShuffleG[tofix], y] /. disp];
  If[! FreeQ[dy /. G[__, y] -> 1, x], Print["  WARNING: x left in y-derivative of ", g]];
  cand = tofix + GIntegrate[dy, y];
  d1 = N[gn[g - cand, pt1], 50];
  If[Abs[Im[d1]] > 10^-40, Print["  WARNING: complex constant ", N[d1, 10], " for ", g]];
  cb = consts[w];
  c = If[Abs[d1] < 10^-40, 0,
    nv = FindIntegerNullVector[N[Prepend[cb, Re[d1]], 50], 10^6];
    If[nv[[1]] == 0, Print["  WARNING: PSLQ failed for ", g]; Return[$Failed]];
    -Rest[nv].cb/nv[[1]]];
  cand = cand + c;
  chk = Abs@N[gn[g - cand, pt2], 30];
  If[chk > 10^-25, Print["  WARNING: check at 2nd point failed (", N[chk, 3], ") for ", g]; Return[$Failed]];
  cand];

ext = {};
Do[Module[{t = AbsoluteTime[], r = fibrate[todo[[k]]]},
   If[r =!= $Failed,
    AppendTo[ext, todo[[k]] -> r]; known[todo[[k]]] = True;
    disp = Dispatch[Join[fibr, ext]]];
   Print[k, "/", Length[todo], "  w", Length[todo[[k]]] - 1, "  ", todo[[k]], "  ",
    If[r === $Failed, "FAILED", "ok"], "  (", Round[AbsoluteTime[] - t], "s)"];
   If[Mod[k, 10] == 0, Put[ext, FileNameJoin[{dir, "fibrXY_ext2.m"}]]]],
  {k, Length[todo]}];
Put[ext, FileNameJoin[{dir, "fibrXY_ext2.m"}]];
Print["wrote fibrXY_ext2.m: ", Length[ext], " of ", Length[todo], " rules; total ", Round[AbsoluteTime[] - t0], "s"];
