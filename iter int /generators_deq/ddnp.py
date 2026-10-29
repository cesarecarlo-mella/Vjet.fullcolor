"""double-double arithmetic on numpy arrays (about 32 digits): DD (real) and CD (complex), for the quad assembly of the
remainders.  Products by Dekker splitting (numpy has no fma).  Sums over an axis by pairwise (tree) reduction."""
import numpy as np
import mpmath as mp

_SPLIT = 134217729.0                                    # 2^27 + 1


def _two_sum(a, b):
    s = a + b
    bb = s - a
    return s, (a - (s - bb)) + (b - bb)


def _split(a):
    t = _SPLIT * a
    hi = t - (t - a)
    return hi, a - hi


def _two_prod(a, b):
    p = a * b
    ah, al = _split(a); bh, bl = _split(b)
    return p, ((ah * bh - p) + ah * bl + al * bh) + al * bl


def _norm(s, e):
    h = s + e
    return h, e - (h - s)


class DD:
    __slots__ = ('hi', 'lo')
    __array_priority__ = 100

    def __init__(self, hi, lo=None):
        self.hi = np.asarray(hi, float)
        self.lo = np.zeros_like(self.hi) if lo is None else np.asarray(lo, float)

    @staticmethod
    def of(v):
        return v if isinstance(v, DD) else DD(v)

    @staticmethod
    def from_mp(v):
        hi = float(v)
        return DD(hi, float(mp.mpf(v) - hi))

    def __add__(self, o):
        o = DD.of(o)
        s, e = _two_sum(self.hi, o.hi)
        t, f = _two_sum(self.lo, o.lo)
        s, e = _norm(s, e + t)
        return DD(*_norm(s, e + f))
    __radd__ = __add__

    def __neg__(self):
        return DD(-self.hi, -self.lo)

    def __sub__(self, o):
        return self + (-DD.of(o))

    def __rsub__(self, o):
        return DD.of(o) + (-self)

    def __mul__(self, o):
        o = DD.of(o)
        p, e = _two_prod(self.hi, o.hi)
        return DD(*_norm(p, e + (self.hi * o.lo + self.lo * o.hi)))
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = DD.of(o)
        q1 = self.hi / o.hi
        r = self - o * q1
        q2 = r.hi / o.hi
        r = r - o * q2
        q3 = r.hi / o.hi
        return DD(*_norm(q1, q2)) + q3

    def __rtruediv__(self, o):
        return DD.of(o) / self

    def __pow__(self, n):
        assert isinstance(n, int)
        if n < 0:
            return 1.0 / (self ** (-n))
        out, b = DD(np.ones_like(self.hi)), self
        while n:
            if n & 1:
                out = out * b
            b, n = b * b, n >> 1
        return out

    def sqrt(self):
        x = 1.0 / np.sqrt(self.hi)
        ax = self.hi * x
        p, e = _two_prod(ax, ax)
        return DD(*_two_sum(ax, (self - DD(p, e)).hi * (x * 0.5)))

    def __getitem__(self, k):
        return DD(self.hi[k], self.lo[k])

    def to_float(self):
        return self.hi + self.lo


class CD:
    """complex double-double"""
    __slots__ = ('re', 'im')

    def __init__(self, re, im=None):
        self.re = DD.of(re)
        self.im = DD(np.zeros_like(self.re.hi)) if im is None else DD.of(im)

    @staticmethod
    def from_mp(vals):
        vals = [mp.mpc(v) for v in vals]
        R = [DD.from_mp(v.real) for v in vals]; I = [DD.from_mp(v.imag) for v in vals]
        return CD(DD([r.hi for r in R], [r.lo for r in R]), DD([i.hi for i in I], [i.lo for i in I]))

    def __add__(self, o):
        return CD(self.re + o.re, self.im + o.im)

    def __sub__(self, o):
        return CD(self.re - o.re, self.im - o.im)

    def __mul__(self, o):
        if isinstance(o, CD):
            return CD(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
        return CD(self.re * o, self.im * o)                # o real (DD or float)
    __rmul__ = __mul__

    def __getitem__(self, k):
        return CD(self.re[k], self.im[k])

    def to_complex(self):
        return self.re.to_float() + 1j * self.im.to_float()

    def to_mp(self):
        return np.vectorize(lambda a, b, c, d: mp.mpc(mp.mpf(a) + b, mp.mpf(c) + d), otypes=[object])(
            self.re.hi, self.re.lo, self.im.hi, self.im.lo)


def dd_sum(x, axis=-1):
    """sum of a DD or CD array along an axis (tree reduction)"""
    if isinstance(x, CD):
        return CD(dd_sum(x.re, axis), dd_sum(x.im, axis))
    hi, lo = np.moveaxis(x.hi, axis, 0), np.moveaxis(x.lo, axis, 0)
    v = DD(hi, lo)
    while v.hi.shape[0] > 1:
        n = v.hi.shape[0]
        if n % 2:
            v = DD(np.concatenate([v.hi, np.zeros_like(v.hi[:1])]), np.concatenate([v.lo, np.zeros_like(v.lo[:1])]))
            n += 1
        v = v[: n // 2] + v[n // 2:]
    return v[0]
