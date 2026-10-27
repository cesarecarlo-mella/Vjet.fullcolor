(* divisor_transport.wls -- boundary constants of a high-energy corner in the OTHER ordering, by
   transport along the exceptional divisor of the corner (prototype for the H family).

   CC chart (region z > 1, x, y < 0):  x = -w/v, y = -u/v, z = 1/v, u + v + w = 1.
   Corner CC2: u, v -> 0 (w -> 1).  Blow-up:  u = s p, v = (1 - s) p, w = 1 - p.
     s -> 0 : u << v  (the CC2 ordering; constants known at all weights, transport.boundaries)
     s -> 1 : v << u  (the CC2e ordering; constants wanted)
   On the divisor p -> 0 every letter becomes p^k g_k(s) (the even ones: g in {s, 1-s, 1-2s}; the odd
   l8 -> Log[(q + I)/(q - I)], q = Sqrt[s/(1-s)], with r = -Sqrt|xyz|; l7 -> 0).  The log-free part
   F0(s) obeys  dF0/ds = A(s) F0,  A = Sum_k atilde_k d/ds Log[g_k]   (eps absorbed in the weight).
   With s = q^2/(1+q^2) the alphabet is {0, 1, -1, I, -I} in q, so
       F0^(w)(q) = Sum_a M_a . Int_0^q dq'/(q' - a) F0^(w-1)(q') + C^(w)_CC2
   (the integral prepends the letter a to every G word; tangential base point at q = 0 = s = 0).
   CC2e constants: regularised q -> oo limit (drop Log[q], i.e. Log[1 - s] -> 0), evaluated
   numerically (30 digits) with  R(a, w) = G(a, w; 1) + Int_1^oo [Delta_w(t)/t + a G(w; t)/(t (t - a))] dt,
   Delta_w = G(w; t) minus its Log[t] asymptotics (weights <= 2 here; letters off the positive real
   axis).  Compared with the fitted cc_series/inputs/bdy_CC2e_w<w>.m (w <= 2).

   usage: wolframscript -file divisor_transport.wls [WMAX = 2] *)
SetEnvironment["PATH" -> GetEnvironment["PATH"][[2]] <> ":/usr/local/bin:/opt/homebrew/bin"];
Get["PolyLogTools`"];                     (* always first *)
wmax = If[Length[$ScriptCommandLine] > 1, ToExpression[$ScriptCommandLine[[2]]], 2];
cc = "/Users/cesarecarlomella/Desktop/Physics/PhD/XjetProject/INTEGRALS/4P1M/3L1M/series.all/analytic.continuation/";
atilde = Get[cc <> "cc_series/inputs/atilde.m"];
nf = Length[atilde];
bdy[w_] := Get[If[w == 0, cc <> "transport.boundaries/H/dist/bdy_CC2_w0.m", cc <> "cc_series/inputs/bdy_CC2_w" <> ToString[w] <> ".m"]];
t0 = AbsoluteTime[];

(* ---- the letters on the divisor, as dlogs in q ---- *)
L = <|1 -> x, 2 -> y, 3 -> 1 - x - y, 4 -> 1 - x, 5 -> 1 - y, 6 -> x + y, 9 -> 1 - (1 + x) y, 10 -> 1 - x (1 + y),
      11 -> x - y + y^2, 12 -> -x + x^2 + y, 13 -> (x - 1)^2 - y, 16 -> -x + (1 - y)^2|>;
chart = {x -> -w/v, y -> -u/v} /. {u -> s p, v -> (1 - s) p, w -> 1 - p};
sOfq = q^2/(1 + q^2);
dlogq = Association[KeyValueMap[Function[{k, l},
     Module[{e = Together[l /. chart], num, den, g},
       num = Numerator[e]; den = Denominator[e];
       g = Coefficient[num, p, Exponent[num, p, Min]] / Coefficient[den, p, Exponent[den, p, Min]];
       k -> Apart[Together[D[Log[g /. s -> sOfq], q]], q]]], L]];
dlogq[7] = 0;
dlogq[8] = Apart[D[Log[(q + I)/(q - I)], q], q];
letters = {0, 1, -1, I, -I};
coefOf[expr_, a_] := Together[Limit[(q - a) expr, q -> a]];
dl = Table[Association[Table[a -> coefOf[dlogq[k], a], {a, letters}]], {k, 1, 16}];
Do[If[KeyExistsQ[dlogq, k] && Simplify[dlogq[k] - Total[Table[dl[[k]][a]/(q - a), {a, letters}]]] =!= 0,
     Print["!! dlog of l", k, " not in the alphabet"]; Exit[1]], {k, Keys[dlogq]}];
Print["divisor dlogs in q:"];
Do[Print["  l", k, ": ", Total[Table[dl[[k]][a]/(q - a), {a, letters}]]], {k, Sort[Keys[dlogq]]}];

(* ---- M_a = Sum_k atilde_k dl_k(a) ---- *)
Mof[a_] := SparseArray[atilde /. l[k_] :> dl[[k]][a]];
Mc = Association[Table[a -> Mof[a], {a, letters}]];
Print["connection matrices built (", Round[AbsoluteTime[] - t0, 0.1], " s)"];

(* ---- transport: F^(w) as vectors of G[..., q] ---- *)
integrate[vec_, a_] := Module[{v = Expand[vec]},
   Map[Function[e, If[e === 0, 0,
       Total[Map[If[FreeQ[#, G], # G[a, q], # /. G[ws__, q] :> G[a, ws, q]] &, If[Head[e] === Plus, List @@ e, {e}]]]]], v]];
F[0] = bdy[0];
Do[
  F[w] = Expand[Sum[Mc[a] . integrate[F[w - 1], a], {a, letters}] + bdy[w]];
  Print["weight ", w, ": ", Length[Union[Cases[F[w], _G, Infinity]]], " distinct G  (", Round[AbsoluteTime[] - t0, 0.1], " s)"],
  {w, 1, wmax}];

(* ---- regularised limit q -> oo ---- *)
Do[F[w] = Collect[F[w], _G, Simplify], {w, 1, wmax}];
words = Select[Union[Cases[Table[F[w], {w, 1, wmax}], _G, Infinity]],
   Function[g, ! AllTrue[Flatten[Table[Coefficient[F[w], g], {w, 1, wmax}]], PossibleZeroQ]]];
Print["letters in the transported functions: ", Union[Flatten[List @@@ words /. q -> Nothing]]];
band = Select[words, MemberQ[Most[List @@ #], 1] &];
If[band =!= {}, Print["!! words with the band letter 1 (q = 1, y = -1): ", band,
   "; their coefficients vanish in F: ", AllTrue[band, Function[g, AllTrue[Flatten[Table[Coefficient[F[w], g], {w, 1, wmax}]], PossibleZeroQ]]]]];
(* exact regularised limit at q = oo: T(w; r) = G(w; 1/r) by inversion,
     dq/(q - a) = -dr/r + dr/(r - 1/a)  (a != 0),   dq/q = -dr/r,
   T(a, w'; r) = G(a, w'; 1) - S(1) + S(r),  S(r) = Int_0^r k_a T(w'; r'),  and the regularised
   r -> 0 limit (Log[r] -> 0, all other G(..; r) -> 0) is  R(a, w') = G(a, w'; 1) - S(1).
   Everything ends up as G(...; 1) with letters {0, +-I} (closed under a -> 1/a). *)
intR[a_, e_] := Module[{t = Expand[e]},
   Total[Map[Function[term,
       If[FreeQ[term, G[___, r]],
         term (-G[0, r] + If[a === 0, 0, G[1/a, r]]),
         term /. G[u___, r] :> (-G[0, u, r] + If[a === 0, 0, G[1/a, u, r]])]],
     If[Head[t] === Plus, List @@ t, {t}]]]];
Tinv[{}] = 1;
Tinv[{a_, rest___}] := Tinv[{a, rest}] = Module[{S = intR[a, Tinv[{rest}]]},
    Expand[G[a, rest, 1] - (S /. G[u__, r] :> G[u, 1]) + S]];
Rexact[G[ws__, q]] := Expand[Tinv[{ws}] /. G[__, r] -> 0] /. G[0 .., 1] -> 0;
regInf[g_] := Rexact[g];
rules = Table[g -> regInf[g], {g, words}];
Print["regularised limits of ", Length[words], " G's  (", Round[AbsoluteTime[] - t0, 0.1], " s)"];
(* words with the band letter 1 give divergent G(1, ...; 1) at the split point q = 1 = r; the total is
   regular there (F_182 vanishes at y = -1), so a consistent shuffle regularisation is exact:
   ShuffleRegulate, then the divergent HPL[{1}, 1] = -Log[1 - 1] -> 0.  Its coefficients must cancel. *)
rules = rules /. g_G :> ShuffleRegulate[g] /; MatchQ[g, G[1, ___, 1]];
Do[
  dv = Expand[F[w] /. rules];
  dcoef = Table[Coefficient[dv, HPL[{1}, 1], k], {k, 1, 3}];
  dmax = Max[Abs[N[Flatten[dcoef] /. g_G :> Ginsh[g, {}], 20]] /. Indeterminate -> 99];
  Print["weight ", w, ": coefficient of the divergence Log[1 - 1] (must cancel): max ", dmax],
  {w, 1, wmax}];
rules = rules /. HPL[{1}, 1] -> 0;
consts = Union[Cases[Table[F[w] /. rules, {w, 1, wmax}], G[__, 1], Infinity]];
cval = Table[g -> Ginsh[g, {}], {g, consts}];
Print[Length[consts], " constants G(...; 1), e.g. ", Take[consts, UpTo[3]], " -> HPL: ",
      Take[GToHPL /@ consts, UpTo[3]], "  (", Round[AbsoluteTime[] - t0, 0.1], " s)"];
Do[
  c = Expand[F[w] /. rules] /. Log[n_Integer /; n > 1] :> Total[#2 Log[#1] & @@@ FactorInteger[n]];
  c = Simplify /@ Expand[c];                (* level-4 constants (Catalan, Log[2]) must cancel *)
  Print["weight ", w, ": constants in the result: ",
        Union[Cases[c, Pi | Catalan | _Zeta | _Log | _PolyLog | _G | _HPL, Infinity]]];
  Put[c, FileNameJoin[{DirectoryName[$InputFileName], "results", "bdy_CC2e_H_w" <> ToString[w] <> "_transported.m"}]];
  cn = N[c /. cval, 20];
  refFile = cc <> "cc_series/inputs/bdy_CC2e_w" <> ToString[w] <> ".m";
  If[FileExistsQ[refFile],
    d = Chop[cn - N[Get[refFile], 20], 10^-12];
    nd = Count[d, x_ /; x =!= 0];
    Print["weight ", w, ": transported vs fitted CC2e constants: ", nd, " of ", nf, " components differ",
          If[nd > 0, "  e.g. " <> ToString[Take[Select[Transpose[{Range[nf], cn, N[Get[refFile], 20]}], Chop[#[[2]] - #[[3]], 10^-12] =!= 0 &], UpTo[3]], InputForm], ""]],
    Print["weight ", w, ": no fitted reference (written to results/); e.g. components 1-3: ", N[Take[cn, 3], 12]]],
  {w, 1, wmax}];
Print["done (", Round[AbsoluteTime[] - t0, 0.1], " s)"];
