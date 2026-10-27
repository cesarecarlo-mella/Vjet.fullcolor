#!/usr/bin/env python3
"""Letter statistics of the eps^0 (finite) parts of the denominator-bracket export.

Usage: python3 letter_stats.py [folder]   (default: parent folder of this script)
Writes summary.txt next to this script and checks every constraint below on
every single term (exit code 1 if any counterexample is found).

Conventions: II[a1,...,an] has weight n, a1 = first integration (leftmost).
Letters 1-6 linear, 7-8 algebraic (odd), >=9 quadratic.
"Effective depth" of an entry = position in the word + transcendental weight
of the constant prefactor (Pi^2 counts 2), i.e. constants are read as the
leftmost (first) integrations.
"""
import re, glob, os, sys, collections
D = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'summary.txt')
CW = {'': 0, 'I': 0, 'Pi^2': 2, 'Zeta[3]': 3, 'Pi^4': 4, 'Zeta[5]': 5, 'Pi^2*Zeta[3]': 5, 'Zeta[3]^2': 6, 'Pi^6': 6}
LIN, ALG = set(range(1, 7)), {7, 8}
PARTNER = {11: 5, 16: 5, 12: 4, 13: 4, 14: 6, 15: 6}

rows = []  # (amp, const, word, denominator)
for f in sorted(glob.glob(os.path.join(D, '*.wl'))):
    amp = os.path.basename(f)[:-3]
    L = open(f).read().splitlines()
    i = next(k for k, l in enumerate(L) if l.startswith('+eps^(0)*'))
    buf = []
    for l in L[i + 1:]:
        m = re.match(r'^\+\(([-0-9/]+)\)\*\((.*)\)\*\((.*)\)$', l)
        if m:
            c = m.group(3); w = re.search(r'II\[([0-9,]*)\]', c).group(1)
            buf.append((c[:c.index('II[')].rstrip('*'), tuple(int(x) for x in w.split(',')) if w else ()))
        elif l.startswith(')/('):
            rows += [(amp, c, w, l[3:-1]) for c, w in buf]; buf = []
    assert not buf, amp

out = []
def p(*a): out.append(' '.join(str(x) for x in a))
W = collections.defaultdict(collections.Counter)
for a, c, w, d in rows: W[w][c] += 1
ALL = sorted({l for w in W for l in w})
QUAD = {l for l in ALL if l >= 9}
p(f"eps^0 terms: {len(rows)} in {len({r[0] for r in rows})} amplitudes; distinct II words: {len(W)}")
p("distinct words by weight:", dict(sorted(collections.Counter(len(w) for w in W).items())))
p("letters present:", ALL, " absent among 1..16:", sorted(set(range(1, 17)) - set(ALL)))
p("constants x II-weight:", dict(collections.Counter((c, len(w)) for a, c, w, d in rows).most_common()))
p("\n# distinct words containing letter, weights 1..6")
for l in ALL: p(f"  {l:>2}", [sum(1 for w in W if len(w) == n and l in w) for n in range(1, 7)])
p("\n# position (from left) histogram per letter, distinct words: {weight: [pos1..posn]}")
for l in ALL:
    h = collections.Counter((len(w), i + 1) for w in W for i, x in enumerate(w) if x == l)
    p(f"  {l:>2}", {n: [h[(n, i)] for i in range(1, n + 1)] for n in range(1, 7)})
p("\n# term occurrences per letter (a term counts once per letter it contains)")
p(" ", dict(sorted(collections.Counter(l for a, c, w, d in rows for l in set(w)).items())))
p("\n# adjacency counts (row = entry i, column = entry i+1), distinct words")
A = collections.Counter((w[i], w[i + 1]) for w in W for i in range(len(w) - 1))
p("     " + "".join(f"{c:>6}" for c in ALL))
for r in ALL: p(f"  {r:>2} " + "".join(f"{A[(r, c)]:>6}" for c in ALL))
P = {(w, c) for w in W for c in W[w]}
p("\n# first two entries, words with no constant prefactor:", dict(sorted(collections.Counter(w[:2] for w, c in P if c == '' and len(w) > 1).items())))
p("# first entry by constant prefactor:", {c: dict(collections.Counter(w[0] for w, cc in P if cc == c and w)) for c in CW})
p("\n# per structure / helicity")
H = collections.defaultdict(set); N = collections.Counter()
for a, c, w, d in rows: H[a] |= set(w); N[a] += 1
for s in sorted({a[:-4] for a in H}):
    sets = [tuple(sorted(H[f"{s}_HA{h}"] - LIN)) for h in range(1, 7)]
    p(f"  {s}: terms {sum(N[f'{s}_HA{h}'] for h in range(1, 7)):>7}  non-linear letters HA1..6:", sets)
p("\n# algebraic letters 7, 8")
alg = [w for w in W if ALG & set(w)]
p("  words:", len(alg), " (#l7,#l8):", dict(collections.Counter((w.count(7), w.count(8)) for w in alg)))
p("  entries 1-2 -> entry 3:", dict(collections.Counter((w[:2], w[2]) for w in alg)))
p("  entries 3-4:", dict(collections.Counter(w[2:4] for w in alg if len(w) > 3)))
p("  entries 5-6:", dict(collections.Counter(w[4:] for w in alg if len(w) == 6)))
p("  prefactors:", dict(sum((W[w] for w in alg), collections.Counter())))
p("\n# quadratic letters")
for q in sorted(QUAD):
    ws = [(w, c) for w, c in P if q in w]
    p(f"  {q}: words {len({w for w, c in ws})}, eff. depth {dict(collections.Counter(w.index(q) + 1 + CW[c] for w, c in ws))},"
      f" before {sorted({w[w.index(q)-1] for w, c in ws if w.index(q) > 0})}, after {sorted({w[w.index(q)+1] for w, c in ws if w.index(q) < len(w)-1})},"
      f" prefactors {dict(collections.Counter(c for w, c in ws))}, #distinct prefixes {len({(c, w[:w.index(q)]) for w, c in ws})}")
def facs(d): return {re.sub(r'\^\d+$', '', x) for x in re.findall(r'\(([^()]*)\)', d)}
p("\n# denominator factors by word type (term counts)")
for name, sel in [("linear-only", lambda w: not (QUAD | ALG) & set(w)), ("with l7/l8", lambda w: ALG & set(w)), ("with quadratic", lambda w: QUAD & set(w))]:
    p(f"  {name}:", dict(collections.Counter(f for a, c, w, d in rows if sel(w) for f in facs(d))))

p("\n# bare II length (constant prefactor ignored)")
T = collections.Counter(len(w) for a, c, w, d in rows)
p("  terms by length:", dict(sorted(T.items())))
p("  distinct words by length:", dict(sorted(collections.Counter(len(w) for w in W).items())))
p("  terms by (length, prefactor):", dict(sorted(collections.Counter((len(w), c) for a, c, w, d in rows).items())))
WS = collections.defaultdict(set)
for a, c, w, d in rows: WS[a[:-4]].add(w)
p("  per structure, distinct words by length 0..6 and max length:")
for s in sorted(WS):
    h = collections.Counter(len(w) for w in WS[s]); p(f"    {s}", [h[n] for n in range(7)], "max", max(h))
for name, sel in [("linear only", lambda w: not (QUAD | ALG) & set(w)), ("with l7/l8", lambda w: ALG & set(w)), ("with quadratic", lambda w: QUAD & set(w))]:
    p(f"  {name}: distinct words by length", dict(sorted(collections.Counter(len(w) for w in W if sel(w)).items())))
for q in sorted(QUAD): p(f"  letter {q}: distinct words by length", dict(sorted(collections.Counter(len(w) for w in W if q in w).items())))

# ---------------- constraints, checked term by term ----------------
def eff(w, c, i): return i + 1 + CW[c]
C = [
 ("letters 9 and 10 never appear", lambda a, c, w, d: not {9, 10} & set(w)),
 ("no constant => first entry in {1,2,3}", lambda a, c, w, d: c != '' or not w or w[0] in (1, 2, 3)),
 ("first entry never 6, 7, 8 or quadratic", lambda a, c, w, d: not w or w[0] <= 5),
 ("no constant: second entry 4 only after 1, 5 only after 2, 6 only after 3",
     lambda a, c, w, d: c != '' or len(w) < 2 or w[1] not in (4, 5, 6) or w[0] == w[1] - 3),
 ("at most one quadratic letter per word", lambda a, c, w, d: sum(l in QUAD for l in w) <= 1),
 ("quadratic and algebraic letters never in the same word", lambda a, c, w, d: not (QUAD & set(w) and ALG & set(w))),
 ("quadratic words carry no constant or Pi^2 only", lambda a, c, w, d: not QUAD & set(w) or c in ('', 'Pi^2')),
 ("quadratic letter sits at effective depth 4", lambda a, c, w, d: all(eff(w, c, i) == 4 for i, l in enumerate(w) if l in QUAD)),
 ("neighbours of quadratic letter q are in {1,2,3,partner(q)}; partner 11,16->5, 12,13->4, 14,15->6",
     lambda a, c, w, d: all(x in (1, 2, 3, PARTNER[l]) for i, l in enumerate(w) if l in QUAD for x in w[max(i-1, 0):i] + w[i+1:i+2])),
 ("l7/l8: weight>=4 words have exactly two odd letters, at entries 3 and 4",
     lambda a, c, w, d: not ALG & set(w) or len(w) < 4 or [i for i, l in enumerate(w) if l in ALG] == [2, 3]),
 ("l7/l8: weight-3 words have one odd letter, at entry 3, prefactor I; I only there",
     lambda a, c, w, d: (c == 'I') == (len(w) == 3 and bool(ALG & set(w))) and (len(w) != 3 or not ALG & set(w) or [i for i, l in enumerate(w) if l in ALG] == [2])),
 ("l7/l8 words: entries 1,2 are two different letters of {1,2,3}",
     lambda a, c, w, d: not ALG & set(w) or (w[0] in (1, 2, 3) and w[1] in (1, 2, 3) and w[0] != w[1])),
 ("entry 3 = l7 => 3 in entries 1-2; entry 3 = l8 => 2 in entries 1-2",
     lambda a, c, w, d: not ALG & set(w) or ((w[2] == 7 and 3 in w[:2]) or (w[2] == 8 and 2 in w[:2]))),
 ("l7/l8 words: entry 5 in {1,2,3}", lambda a, c, w, d: not ALG & set(w) or len(w) < 5 or w[4] in (1, 2, 3)),
 ("l7/l8 only in C1F0Z0, C7F0Z0, C9F0Z0", lambda a, c, w, d: not ALG & set(w) or a[:-4] in ('C1F0Z0', 'C7F0Z0', 'C9F0Z0')),
 ("quadratic denominators t^2-t+u, t+u^2-u multiply linear-only words",
     lambda a, c, w, d: not ('t^2-t+u' in d or 't+u^2-u' in d) or set(w) <= LIN),
]
p("\n# CONSTRAINTS (checked on every term)")
fail = False
for name, f in C:
    bad = [r for r in rows if not f(*r)]
    p(f"  [{'OK' if not bad else 'FAIL'}] {name}" + (f"  -- {len(bad)} counterexamples, e.g. {bad[:3]}" if bad else ""))
    fail |= bool(bad)
open(OUT, 'w').write('\n'.join(out) + '\n')
print('\n'.join(out[-len(C)-1:])); print("written", OUT)
sys.exit(1 if fail else 0)
