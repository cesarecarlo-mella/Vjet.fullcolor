"""
iisys.py -- canonical system closed under d for a set of pure functions sum c II[a1..an] over the 20-letter
alphabet (2610.anc/alphabet.wl order; d II[a1..an] = dlog W_an II[a1..a(n-1)]), and its numerical evaluation with
the Czakon-Tancredi method of trident_deq (tridents_bs.cpp, unchanged).

    closure : V_w = span(functions of weight w, strips of V_{w+1}),  w = 5..1,   V_0 = {1}
    numerics: y_w = Q_w^T II_w with Q_w an orthonormal basis of V_w in word space  ->
              dy_w = sum_k dlog W_k A[k,w] y_{w-1},   A[k,w] = (Q_{w-1}^T S_k Q_w)^T,   |A| <= 1
              (S_k strips a trailing k); the functions are read off: dF = sum_k dlog W_k R[k,w] y_{w-1}
    base    : corner c1, shuffle regularised (II[] = 1, every non-empty II -> 0), nested: x -> 0 first
              ('x'), or y -> 0 first ('y'); start values by 1-d log-power series along the edge and then
              transversally (in sqrt of the small variable: the odd letters go like sqrt there).
"""
import collections, re, struct
from fractions import Fraction
import numpy as np
import mpmath as mp

P = 2**61 - 1
NL = 20
# plain letters as {(i, j): c} for x^i y^j  (as in tridents_bs.cpp)
PLAIN = {1: {(1, 0): 1}, 2: {(0, 1): 1}, 3: {(0, 0): 1, (1, 0): -1, (0, 1): -1}, 4: {(0, 0): 1, (1, 0): -1},
         5: {(0, 0): 1, (0, 1): -1}, 6: {(1, 0): 1, (0, 1): 1}, 9: {(0, 0): 1, (0, 1): -1, (1, 1): -1},
         10: {(0, 0): 1, (1, 0): -1, (1, 1): -1}, 11: {(1, 0): 1, (0, 1): -1, (0, 2): 1},
         12: {(1, 0): -1, (2, 0): 1, (0, 1): 1}, 13: {(0, 0): 1, (1, 0): -2, (2, 0): 1, (0, 1): -1},
         14: {(0, 1): -1, (2, 0): 1, (1, 1): 2, (0, 2): 1}, 15: {(1, 0): -1, (2, 0): 1, (1, 1): 2, (0, 2): 1},
         16: {(1, 0): -1, (0, 0): 1, (0, 1): -2, (0, 2): 1}, 17: {(1, 0): 1, (1, 1): 1, (0, 2): 1},
         18: {(2, 0): 1, (0, 1): 1, (1, 1): 1}, 19: {(0, 0): 1, (0, 1): -2, (0, 2): 1, (1, 1): 1},
         20: {(0, 0): 1, (1, 0): -2, (2, 0): 1, (1, 1): 1}}
TERM = re.compile(r'\+\((-?[\d/]+)\)\*\(II\[([\d, ]*)\]\)')


def modp(c):
    return c.numerator % P * pow(c.denominator % P, P - 2, P) % P


def parse_generators(path):
    """Generators.wl: '(* NAME *)' followed by '+(c)*(II[...])' terms -> [(name, {word: Fraction})]"""
    src = open(path).read()
    blocks = re.split(r'\(\*\s*(\w+\[\d+\])\s*\*\)', src)[1:]
    out = []
    for name, body in zip(blocks[::2], blocks[1::2]):
        f = collections.defaultdict(Fraction)
        for c, w in TERM.findall(body):
            f[tuple(int(a) for a in w.split(',')) if w.strip() else ()] += Fraction(c)
        out.append((name, {w: c for w, c in f.items() if c}))
    return out


def strip(f, k):
    out = collections.defaultdict(Fraction)
    for w, c in f.items():
        if w and w[-1] == k:
            out[w[:-1]] += c
    return {w: c for w, c in out.items() if c}


def build(funcs, verbose=True):
    """funcs: [(name, {word: Fraction})], each of one weight <= 6.  Returns a dict with the numerical system."""
    wt = [len(next(iter(f))) for _, f in funcs]
    top = max(wt)
    assert top <= 6
    basis = {}                                       # w -> list of exact vectors (dict word -> Fraction)
    origin = {}                                      # w -> index of the input each basis function is (None: a strip)
    above = []
    for w in range(top, 0, -1):
        gidx = [i for i, n in enumerate(wt) if n == w]
        cand = [funcs[i][1] for i in gidx]
        ngen = len(cand)
        cand += [g for f in above for k in range(1, NL + 1) for g in [strip(f, k)] if g]
        piv, keep, src = [], [], []
        for ic, f in enumerate(cand):
            r = {word: modp(c) for word, c in f.items()}
            for col, pr in piv:
                if r.get(col):
                    a = r[col]
                    for c2, v2 in pr.items():
                        r[c2] = (r.get(c2, 0) - a * v2) % P
            r = {c: v for c, v in r.items() if v}
            if r:
                col = min(r)
                inv = pow(r[col], P - 2, P)
                piv.append((col, {c: v * inv % P for c, v in r.items()}))
                keep.append(f)
                src.append(gidx[ic] if ic < ngen else None)
        basis[w] = keep
        origin[w] = src
        above = keep
        if verbose:
            print('  weight %d: %3d functions (%d inputs)' % (w, len(keep), ngen), flush=True)
    for w in range(top + 1, 6):
        basis[w] = []; origin[w] = []
    basis[0] = [{(): Fraction(1)}]
    origin[0] = [None]
    # word supports and orthonormal bases (weights 0..5)
    sup, Q = {}, {}
    for w in range(0, 6):
        words = sorted({word for f in basis[w] for word in f})
        sup[w] = {word: i for i, word in enumerate(words)}
        if not basis[w]:
            Q[w] = np.zeros((len(words), 0))
            continue
        M = np.zeros((len(words), len(basis[w])))
        for j, f in enumerate(basis[w]):
            for word, c in f.items():
                M[sup[w][word], j] = float(c)
        q, r = np.linalg.qr(M)
        d = np.abs(np.diag(r))
        if verbose and w:
            print('  weight %d: QR |r_ii| in [%.2e, %.2e]' % (w, d.min(), d.max()))
        Q[w] = q

    def strip_mat(w, k, X):
        """S_k X for a matrix X on the support of weight w (rows) -> rows on the support of weight w-1"""
        out = np.zeros((len(sup[w - 1]), X.shape[1]))
        for word, i in sup[w].items():
            if word[-1] == k:
                j = sup[w - 1].get(word[:-1])
                if j is None:
                    raise ValueError('support not prefix closed: %s' % (word,))
                out[j] += X[i]
        return out

    dims = [Q[w].shape[1] for w in range(6)]
    A, worst = {}, 0.0
    for w in range(1, 6):
        for k in range(1, NL + 1):
            if dims[w] == 0:
                A[k, w] = np.zeros((0, dims[w - 1]))
                continue
            S = strip_mat(w, k, Q[w])
            C = Q[w - 1].T @ S
            worst = max(worst, np.abs(S - Q[w - 1] @ C).max() if S.size else 0)
            A[k, w] = C.T
    # read-off of the inputs
    NT = len(funcs)
    R = {(k, w): np.zeros((NT, dims[w - 1])) for k in range(1, NL + 1) for w in range(1, 7)}
    const = np.zeros(NT)
    for i, ((name, f), n) in enumerate(zip(funcs, wt)):
        if n == 0:
            const[i] = float(f[()])
            continue
        for k in range(1, NL + 1):
            g = strip(f, k)
            if not g:
                continue
            v = np.zeros(len(sup[n - 1]))
            for word, c in g.items():
                v[sup[n - 1][word]] = float(c)
            cf = Q[n - 1].T @ v
            worst = max(worst, np.abs(v - Q[n - 1] @ cf).max())
            R[k, n][i] = cf
    if verbose:
        print('  dims', dims, ' closure residual %.1e' % worst)
        print('  max |A| entry %.3f, max |R| entry %.3g' % (max(np.abs(a).max() if a.size else 0 for a in A.values()),
                                                           max(np.abs(r).max() for r in R.values())))
    return dict(dims=dims, A=A, R=R, wt=np.array(wt), const=const, names=[n for n, _ in funcs], basis=basis, origin=origin)


# ---------------------------------------------------------------- 1-d series of the forms along a line
def _ser(a, N):
    a = list(a)[:N]
    return a + [mp.mpf(0)] * (N - len(a))


def _mul(a, b, N):
    return [mp.fsum(a[i] * b[n - i] for i in range(n + 1)) for n in range(N)]


def _inv(a, N):
    assert a[0] != 0
    out = [1 / a[0]]
    for n in range(1, N):
        out.append(-mp.fsum(a[i] * out[n - i] for i in range(1, n + 1)) / a[0])
    return out


def _sqrt(a, N):
    """a[0] > 0"""
    out = [mp.sqrt(a[0])]
    for n in range(1, N):
        s = mp.fsum(out[i] * out[n - i] for i in range(1, n))
        out.append((a[n] - s) / (2 * out[0]))
    return out


def _der(a, N):
    return [(n + 1) * a[n + 1] for n in range(N - 1)] + [mp.mpf(0)]


def _val(a):
    for i, c in enumerate(a):
        if abs(c) > mp.mpf(10) ** (-(mp.mp.dps - 5)):
            return i
    return None


def _poly(coef, X, Y, N):
    """sum c x^i y^j for series X, Y"""
    out = [mp.mpf(0)] * N
    Xp, Yp = [[mp.mpf(1)] + [mp.mpf(0)] * (N - 1)], [[mp.mpf(1)] + [mp.mpf(0)] * (N - 1)]
    for _ in range(2):
        Xp.append(_mul(Xp[-1], X, N)); Yp.append(_mul(Yp[-1], Y, N))
    for (i, j), c in coef.items():
        t = _mul(Xp[i], Yp[j], N)
        out = [o + c * u for o, u in zip(out, t)]
    return out


def _dlog(W, N):
    """Laurent coefficients of W'/W: list of length N+1, index 0 <-> s^-1; None if W == 0"""
    v = _val(W)
    if v is None:
        return None
    U = W[v:] + [mp.mpf(0)] * v
    L = _mul(_der(U, N), _inv(U, N), N)
    return [mp.mpf(v)] + L[:N]


def forms_on_line(X, Y, N):
    """Laurent series (in s) of dlog W_k along x = X(s), y = Y(s) (power series): {k: complex array (N+1)},
    index 0 <-> s^-1.  Odd letters with r = +sqrt(xyz) on the real positive side; the form of a letter that is
    identically constant (or 0/0 -> constant) on the line is 0."""
    M = N + 8
    X, Y = _ser(X, M), _ser(Y, M)
    Z = [(1 if i == 0 else 0) - X[i] - Y[i] for i in range(M)]
    out = {}
    for k, coef in PLAIN.items():
        L = _dlog(_poly(coef, X, Y, M), M)
        out[k] = np.zeros(N + 1, complex) if L is None else np.array([complex(c) for c in L[:N + 1]])
    b = _mul(_mul(X, Y, M), Z, M)
    vb = _val(b)
    for k in (7, 8):
        a = _mul(X, Z, M) if k == 7 else _mul(X, Y, M)
        if vb is None:
            out[k] = np.zeros(N + 1, complex)
            continue
        assert vb % 2 == 0, 'odd valuation of xyz on the line'
        beta = b[vb:] + [mp.mpf(0)] * vb
        r = [mp.mpf(0)] * (vb // 2) + _sqrt(beta, M)[:M - vb // 2]           # r = s^(vb/2) sqrt(beta)
        num = [2 * u - v for u, v in zip(_mul(b, _der(a, M), M), _mul(a, _der(b, M), M))]
        den = _mul(r, [u + v for u, v in zip(_mul(a, a, M), b)], M)
        vn, vd = _val(num), _val(den)
        if vn is None:
            out[k] = np.zeros(N + 1, complex)
            continue
        q = _mul(num[vn:] + [0] * vn, _inv(den[vd:] + [0] * vd, M), M)       # s^(vn-vd) q(s)
        sh = vn - vd                                                          # >= -1 for a dlog form
        assert sh >= -1, sh
        L = [mp.mpf(0)] * (N + 1)
        for i in range(N + 1):
            j = i - 1 - sh
            if 0 <= j < M:
                L[i] = q[j]
        out[k] = 1j * np.array([complex(c) for c in L])
    return out


# ---------------------------------------------------------------- log-power series of the system on a line
def propagate(S, forms, s_end, y_init, T_init, N):
    """regularised solution along a line from its (tangential) start s = 0, where y_w = y_init[w], T = T_init
    (T_init: (7, NT)); returns y_w(s_end), T(s_end).  Series in s with powers of L = log s up to 6."""
    dims, A, R = S['dims'], S['A'], S['R']
    NT = R[1, 1].shape[0]
    Lmax = 7
    ls = np.log(s_end)
    sp = s_end ** np.arange(N)
    lp = ls ** np.arange(Lmax)
    # B_j = sum_k c_k[j] A[k,w]: Laurent index j = 0 <-> s^-1
    def integrate(I):
        """I[p+1, m] for powers p = -1..N-2 (arrays of shape (N, Lmax, ...)) -> antiderivative coefficients
        out[p, m] for p = 0..N-1 (regularised: no constant)"""
        out = np.zeros((N, Lmax) + I.shape[2:], complex)
        for m in range(Lmax - 1):
            out[0, m + 1] += I[0, m] / (m + 1)
        for p in range(0, N - 1):
            for m in range(Lmax):
                if not np.any(I[p + 1, m]):
                    continue
                fac = 1.0
                for i in range(m + 1):                    # int s^p L^m = s^(p+1) sum_i (-1)^i m!/(m-i)! L^(m-i)/(p+1)^(i+1)
                    out[p + 1, m - i] += (-1) ** i * fac * I[p + 1, m] / (p + 1) ** (i + 1)
                    fac *= (m - i)
        return out

    def rhs(blocks, Yprev):
        """sum_k forms_k(s) blocks[k] Yprev(s): Yprev (N, Lmax, d) -> (N, Lmax, rows) Laurent-shifted by 1"""
        rows = next(iter(blocks.values())).shape[0]
        I = np.zeros((N, Lmax, rows), complex)
        for k, Bk in blocks.items():
            if not np.any(Bk) or not np.any(forms[k]):
                continue
            Z = Yprev @ Bk.T                               # (N, Lmax, rows)
            c = forms[k]
            for j in range(N):                             # c[j] s^(j-1) * Z[p] s^p -> power p + j - 1 -> index p + j
                if c[j] == 0:
                    continue
                I[j:] += c[j] * Z[:N - j]
        return I

    Y = {0: np.zeros((N, Lmax, 1), complex)}
    Y[0][0, 0, 0] = 1.0
    vals = {0: np.array([1.0 + 0j])}
    T = np.array(T_init, complex).copy()
    for w in range(1, 7):
        if w <= 5 and dims[w]:
            Iy = rhs({k: A[k, w] for k in range(1, NL + 1)}, Y[w - 1])
            Y[w] = integrate(Iy)
            Y[w][0, 0] += y_init[w]
            vals[w] = np.einsum('pmd,p,m->d', Y[w], sp, lp)
        elif w <= 5:
            Y[w] = np.zeros((N, Lmax, 0), complex); vals[w] = np.zeros(0, complex)
        if dims[w - 1]:
            It = rhs({k: R[k, w] for k in range(1, NL + 1)}, Y[w - 1])
            T[w] += np.einsum('pmd,p,m->d', integrate(It), sp, lp)
    return vals, T


def start_values(S, base='x', edge=0.5, eps=0.02, N=60, dps=40, y_init=None, T_init=None):
    """values at the point near c1 reached by the nested regularised path:
       base 'x': x -> 0 first: along x = 0 from y = 0 to y = edge, then x = s^2 from 0 to eps at y = edge
       base 'y': the mirror.   Returns (point, y_w values, T (7, NT)).
       y_init, T_init: the regularised values at the corner (default: those of the base point, 0 above weight 0);
       non-zero for a production region, whose corner is the decay corner continued."""
    mp.mp.dps = dps
    NT = S['R'][1, 1].shape[0]
    zero = [mp.mpf(0)]
    one_s = [mp.mpf(0), mp.mpf(1)]
    if T_init is None:
        T0 = np.zeros((7, NT), complex)
        T0[0] = S['const']
    else:
        T0 = np.array(T_init, complex)
    if y_init is None:
        y0 = {w: np.zeros(S['dims'][w], complex) for w in range(6)}
        y0[0] = np.array([1.0 + 0j])
    else:
        y0 = {w: np.array(y_init[w], complex) for w in range(6)}
    e = mp.mpf(edge)
    s2 = [mp.mpf(0), mp.mpf(0), mp.mpf(1)]
    if base == 'x':
        leg1 = forms_on_line(zero, one_s, N)               # x = 0, y = s
        leg2 = forms_on_line(s2, [e], N)                   # x = s^2, y = edge
        pt = (eps, edge)
    else:
        leg1 = forms_on_line(one_s, zero, N)
        leg2 = forms_on_line([e], s2, N)
        pt = (edge, eps)
    v1, T1 = propagate(S, leg1, edge, y0, T0, N)
    v2, T2 = propagate(S, leg2, np.sqrt(eps), v1, T1, N)
    return pt, v2, T2


# ---------------------------------------------------------------- export for tridents_bs.cpp (TRIDBS01)
def export(S, fn, pt, yv, T, readout_y=False):
    """data file for tridents_bs.cpp; with readout_y the y_w (w = 1..5) are appended to the read-off functions
    (to carry the start values to another point)."""
    dims = S['dims']
    R = {key: v.copy() for key, v in S['R'].items()}
    T = np.array(T, complex)
    if readout_y:
        for w in range(1, 6):
            for k in range(1, NL + 1):
                ext = np.zeros((sum(dims[1:6]), dims[w - 1]))
                o = sum(dims[1:w])
                ext[o:o + dims[w]] = S['A'][k, w]
                R[k, w] = np.vstack([R[k, w], ext])
            for k in range(1, NL + 1):
                pass
        for w2 in (6,):
            for k in range(1, NL + 1):
                R[k, w2] = np.vstack([R[k, w2], np.zeros((sum(dims[1:6]), dims[5]))])
        Ty = np.zeros((7, sum(dims[1:6])), complex)
        for w in range(1, 6):
            o = sum(dims[1:w])
            Ty[w, o:o + dims[w]] = yv[w]
        T = np.hstack([T, Ty])
    NT = T.shape[1]
    x0, y0 = pt
    out = open(fn, 'wb')
    I = lambda *v: out.write(struct.pack('<%di' % len(v), *v))
    C = lambda a: out.write(np.ascontiguousarray(a, dtype=np.complex128).tobytes())

    def block(M):
        if M.size == 0 or not np.any(M):
            I(0, M.shape[0], M.shape[1], 0)
        else:
            I(1, M.shape[0], M.shape[1], 0); C(M)
    out.write(b'TRIDBS01'); I(NL, 6, NT); I(*dims)
    out.write(struct.pack('<2d', x0, y0))
    C(np.array([np.sqrt(x0 * y0 * (1 - x0 - y0)) + 0j]))
    for w in range(6):
        C(yv[w])
    C(T.reshape(-1))
    for k in range(1, NL + 1):
        for w in range(1, 6):
            block(S['A'][k, w])
    for k in range(1, NL + 1):
        for w in range(1, 7):
            block(R[k, w])
    out.close()
    return NT
