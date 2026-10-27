Pipeline used for polylog_fit_results.txt (ran in a cloud workspace; paths inside the scripts point to
/home/claude/fit and to a staged copy of the .wl files in /mnt/user-data/uploads/denominator_bracket;
change DATA in core.py and the json paths to run locally).
  core.py      parse eps^0 terms, rho, delta, pure-function sampling, candidate substitution
  fast.py      rank-one sketches for abstract symbols and data; sk2.py/sk3.py: incremental rref mod p
  gsym*.py     symbols of Li_{n1,n2}, Li_{n1,n2,n3} via Goncharov's I-recursion (abstract letters)
  algargs.py / extargs.py  rational + algebraic admissible arguments, pairs;  triples.py  triples
  full5.py     weight-5 sufficiency test;  w6O.py  depth-3 obstruction at weight 6
  dlevel6d.py / final6c.py  weight-6 tests (delta level / full, incl. spurious-letter combos via spur.py)
  sel4c.py / sel5.py  selection of the sufficient Li_{2,2} / Li_{3,2}, Li_{2,2,1} lists
  pairs.json   (rational args/pairs + abstract symbols), ext_args.json (all args with names, pairs),
  triples.json, Li2_2_2.json, depth3_w5.json, depth3_w6.json, sel_w4.json, sel_w5.json
