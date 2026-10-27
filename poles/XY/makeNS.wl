(* Non-singlet, new pole prediction (../poles/PolePrediction/Amp<i>.m) vs old (poles_delta/Amp<i>.m):
   the difference is G-free (Delta^(3) quartic term, eps^-1, Nc layer only), so the new II-file is
   fc_poles_new_<i>.m = LC_Z_convert_to_II/XY_POLES/fc_poles_<i>.m + (new - old /. z -> 1-x-y). Asserts G-free and eps^-1 only. *)
dir = DirectoryName[$InputFileName];
lc = FileNameJoin[{ParentDirectory[ParentDirectory[dir]], "LC_Z_convert_to_II"}];
Do[Module[{new, old, d, fc},
   new = Get[FileNameJoin[{ParentDirectory[dir], "PolePrediction", "Amp" <> ToString[i] <> ".m"}]];
   old = Get[FileNameJoin[{lc, "poles_delta", "Amp" <> ToString[i] <> ".m"}]];
   d = Factor[Expand[new - old]];
   If[! FreeQ[d, G], Print["ERROR: difference contains G in Amp", i]; Exit[1]];
   If[d =!= 0 && Union[Exponent[#, eps] & /@ If[Head[Expand[d]] === Plus, List @@ Expand[d], {Expand[d]}]] =!= {-1},
     Print["ERROR: difference not pure eps^-1 in Amp", i]; Exit[1]];
   Print["Amp", i, ": new - old = ", d];
   fc = Get[FileNameJoin[{lc, "XY_POLES", "fc_poles_" <> ToString[i] <> ".m"}]];
   Put[fc + Together[d /. z -> 1 - x - y], FileNameJoin[{dir, "fc_poles_new_" <> ToString[i] <> ".m"}]];
   Print["wrote fc_poles_new_", i, ".m"]], {i, 6}];
