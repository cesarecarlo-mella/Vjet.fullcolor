(* For each i:  fc_ns_<i>_shuffled.m  = non-singlet poles (new, with Delta^(3)) + eps^0
                fc_ps_<i>_shuffled.m  = NfZ * (singlet poles + eps^0)   [NfZ not in the PolePrediction files]
                fc_poles_<i>_shuffled.m = sum of the two.   Sources: ../FC_POLES/ *)
dir = DirectoryName[$InputFileName]; src = FileNameJoin[{dir, "..", "FC_POLES"}];
g[s_] := Get[FileNameJoin[{src, s <> ".m"}]];
put[e_, s_] := (Put[e, FileNameJoin[{dir, s <> ".m"}]]; Get[FileNameJoin[{dir, s <> ".m"}]] === e);
Do[Module[{n = ToString[i], ns, ps, ok},
   ns = g["fc_poles_new_" <> n <> "_poles_shuffled"] + g["fc_poles_new_" <> n <> "_eps0_shuffled"];
   ps = NfZ (g["fc_ps_" <> n <> "_poles_shuffled"] + g["fc_ps_" <> n <> "_eps0_shuffled"]);
   ok = {put[ns, "fc_ns_" <> n <> "_shuffled"], put[ps, "fc_ps_" <> n <> "_shuffled"], put[ns + ps, "fc_poles_" <> n <> "_shuffled"]};
   Print["fc_poles_", i, "_shuffled.m  reread identical (ns, ps, total): ", ok]], {i, 6}];
