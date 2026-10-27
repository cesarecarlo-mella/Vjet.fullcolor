Symbol-level checks used in sections 9 of ../eps0_letter_stats.txt
  python3 parse.py                 -> rows.pkl (eps^0 terms; run first, ~10 s)
  python3 integ.py AMP [AMP ...]   integrability of the eps^0 symbols with letters.dict.m (even letters)
  python3 classfit.py CONST W [lin] is rho(weight-W part) in span{rho Li_W(x)}, x over args.json (36 orbits)?
                                     e.g.  classfit.py 1 3   |  classfit.py Pi^2 2
  classical.py                     helpers (rho, delta, pure functions at random points mod 2^61-1)
  args.json                        the 36 admissible orbits as exponent vectors over letters [1..6, 11..16]
