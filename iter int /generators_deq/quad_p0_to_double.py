"""put the quad start values at P0 (out/iter_quad_P0.bin, 27 digits) into the double evaluator: the head of
out/generators_bs.bin (and out/start_y.head, so that the second evaluation differs only by its path).

The double system (iisys.build) integrates y_w = Q_w^T II_w (w = 0..5; Q_w orthonormal in word space, 479 functions)
and reads off the 579 generators.  The 586 functions of the exact system span the same spaces, so Q_w = M_w K_w with
M_w their word vectors, and y_w(P0) = K_w^T F_w(P0).  Q_w is recomputed as in iisys.build and checked against the
matrices A[k, w] of the data file.  The old heads are kept as *_edgeP0.
usage: quad_p0_to_double.py [--write]"""
import os, re, struct, sys, pickle, shutil
from fractions import Fraction
import numpy as np
import mpmath as mp
import flint
import iisys
import quad_data as Q

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'out')
mp.mp.dps = 40
S = pickle.load(open(os.path.join(OUT, 'system.pkl'), 'rb'))
dims = S['dims']

# Q_w as in iisys.build
sup, Qw = {}, {}
for w in range(6):
    words = sorted({word for f in S['basis'][w] for word in f})
    sup[w] = {word: i for i, word in enumerate(words)}
    M = np.zeros((len(words), len(S['basis'][w])))
    for j, f in enumerate(S['basis'][w]):
        for word, c in f.items():
            M[sup[w][word], j] = float(c)
    Qw[w] = np.linalg.qr(M)[0]
worst = 0
for w in range(1, 6):
    for k in range(1, iisys.NL + 1):
        Sk = np.zeros((len(sup[w - 1]), dims[w]))
        for word, i in sup[w].items():
            if word[-1] == k:
                Sk[sup[w - 1][word[:-1]]] += Qw[w][i]
        worst = max(worst, np.abs((Qw[w - 1].T @ Sk).T - S['A'][k, w]).max())
print('recomputed Q_w reproduce the matrices of the data file: max |A - A_file| = %.1e' % worst)
assert worst < 1e-13

def _pivots(A):
    """indices of a maximal set of independent columns of the fmpq_mat A (pivot columns of its rref)"""
    R, rank = A.rref()
    piv, r = [], 0
    for c in range(A.ncols()):
        if r < rank and R[r, c] != 0:
            piv.append(c); r += 1
    return piv


# word vectors of the 586 functions
src = re.sub(r'\(\*.*?\*\)', '', open(os.path.join(OUT, 'basis.wl')).read(), flags=re.S)
vec = {}
for lab, e in re.findall(r'([A-Z]\[[\d,]+\]) -> (.*?)(?=,\s*\n[A-Z]\[|\s*\}\s*$)', src, re.S):
    vec[lab] = {tuple(int(a) for a in word.split(',') if a.strip()): Fraction(c)
                for c, word in re.findall(r'\((-?[\d/]+)\)\*II\[([\d, ]*)\]', e)}
assert len(vec) == Q.NF, len(vec)

def head_values(F, verbose=True):
    """the 503-system start values from the 586 functions F (mp): (states y_0..y_5 concatenated, T (7 x 579))"""
    yv = []
    for w in range(6):
        idx = [i for i, l in enumerate(Q._labels) if vec[l] and len(next(iter(vec[l]))) == w] if w else \
              [Q._labels.index('B[0,1]')]
        if w == 0:
            yv.append(np.array([complex(F[idx[0]] * 96)]))                 # y_0 = II[] = 1 = 96 B[0,1]
            continue
        C, piv, Rm = _map(w, idx)
        b = [mp.fsum(mp.mpf(int(C[j, a].p)) / int(C[j, a].q) * F[idx[c]] for j, c in enumerate(piv) if C[j, a] != 0)
             for a in range(dims[w])]
        yw = mp.lu_solve(Rm.T, mp.matrix(b))
        yv.append(np.array([complex(yw[a]) for a in range(dims[w])]))
        if verbose:
            print('weight %d: %d states from %d functions (%d independent), exact C' % (w, dims[w], len(idx), len(piv)), flush=True)
    T = Q.generators(F)
    Tm = np.zeros((7, len(T)), complex)
    for i, (wgt, v) in enumerate(zip(S['wt'], T)):
        Tm[wgt, i] = complex(v)
    return np.concatenate(yv), Tm


_MAPS = {}


def _map(w, idx):
    """exact: the iisys basis functions b = C^T F (C rational, flint); y_w = R_w^{-T} b with Q_w R_w = M_iisys"""
    if w in _MAPS:
        return _MAPS[w]
    nw = len(sup[w])
    B586 = flint.fmpq_mat(nw, len(idx))
    for j, i in enumerate(idx):
        for word, c in vec[Q._labels[i]].items():
            B586[sup[w][word], j] = flint.fmpq(c.numerator, c.denominator)
    Mi = flint.fmpq_mat(nw, dims[w])
    Mf = np.zeros((nw, dims[w]))
    for j, f in enumerate(S['basis'][w]):
        for word, c in f.items():
            Mi[sup[w][word], j] = flint.fmpq(c.numerator, c.denominator); Mf[sup[w][word], j] = float(c)
    piv = _pivots(B586)
    Bs = flint.fmpq_mat([[B586[r, c] for c in piv] for r in range(nw)])
    Bt = Bs.transpose()
    C = (Bt * Bs).solve(Bt * Mi)                                         # exact: Bs C = M_iisys
    assert Bs * C == Mi, 'not in the span'
    Rw = np.linalg.qr(Mf)[1]
    qr_res = np.abs(Qw[w] @ Rw - Mf).max() / np.abs(Mf).max()
    assert qr_res < 1e-13, qr_res
    _MAPS[w] = (C, piv, mp.matrix([[mp.mpf(float(v)) for v in row] for row in Rw]))
    return _MAPS[w]


def patch_head(fn, new_y, Tm, hub=(0.3, 0.3), write=False, also=()):
    """compare with the head of the data file fn (TRIDBS01: magic, (NL, 6, NT), dims, x0 y0, r0, y_0..y_5, T) and,
    with write, replace it (the old file kept as *_edgeP0.bin); also: further files to receive the new head alone"""
    raw = open(fn, 'rb').read()
    nl, six, nt = struct.unpack('<3i', raw[8:20])
    assert raw[:8] == b'TRIDBS01' and nt == Tm.shape[1] and list(struct.unpack('<6i', raw[20:44])) == dims
    assert struct.unpack('<2d', raw[44:60]) == tuple(hub)
    o = 60 + 16
    old_y = np.frombuffer(raw[o:o + 16 * sum(dims)], complex)
    old_T = np.frombuffer(raw[o + 16 * sum(dims):o + 16 * (sum(dims) + 7 * nt)], complex).reshape(7, nt)
    print('%s: old vs quad start values: states max |diff| %.1e (max |y| %.2g); generators max |diff| %.1e (max |T| %.2g)'
          % (os.path.basename(fn), np.abs(old_y - new_y).max(), np.abs(new_y).max(), np.abs(old_T - Tm).max(), np.abs(Tm).max()))
    off = np.abs(old_T[np.arange(7)[:, None] != S['wt'][None, :]]).max()
    assert off == 0, 'old head: values outside the weight slot of a generator'
    head_len = o + 16 * (sum(dims) + 7 * nt)
    new_head = raw[:o] + new_y.astype(np.complex128).tobytes() + Tm.astype(np.complex128).tobytes()
    assert len(new_head) == head_len
    if write:
        for f in (fn,) + tuple(also):
            bak = f.replace('.bin', '_edgeP0.bin').replace('.head', '_edgeP0.head')
            if not os.path.exists(bak):
                shutil.copy2(f, bak)
        open(fn, 'wb').write(new_head + raw[head_len:])
        for f in also:
            open(f, 'wb').write(new_head)
        if os.path.exists(fn[:-4] + '_y.bin'):                         # cached second data file: rebuilt on use
            os.remove(fn[:-4] + '_y.bin')
        print('written: %s (head)%s; old ones kept as *_edgeP0' % (fn, ''.join(', ' + f for f in also)))


if __name__ == '__main__':
    (inf, F), = Q.run([(0.3, 0.3)], os.path.join(OUT, 'iter_quad_P0.bin'), 'quad')   # zero-length path: the stored values
    new_y, Tm = head_values(F)
    patch_head(os.path.join(OUT, 'generators_bs.bin'), new_y, Tm, write='--write' in sys.argv,
               also=(os.path.join(OUT, 'start_y.head'),))
