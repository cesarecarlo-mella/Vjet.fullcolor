// iter_quad.cpp -- the exact (integer) iterated-integral system integrated in double or double-double precision.
//
//   dF = sum_k dlog W_k M_k F   (F: the functions of out/basis.wl, M_k integer, from out/deq.wl via quad_data.py)
//   along the complex-deformed straight path  x(t) = x0 + (t + 4 i dx t(1-t)) (X - x0),  y(t) likewise (dy),
//   from the start point of the data file (with its values F0 and the root r0 = Sqrt[x y z] there) to (X, Y);
//   r is integrated along (dr = r db / 2b) so that it continues correctly.
//   Bulirsch-Stoer (modified midpoint, polynomial extrapolation in h^2), templated on the real type.
//
// usage: iter_quad DATA.bin [quad|double] [tol] [dx dy [kmax]] < points   (one point "X Y" or "Xhi Xlo Yhi Ylo" per line)
// output per point: "# X Y steps evals seconds Re(r) Im(r)", then one line per function: Re.hi Re.lo Im.hi Im.lo
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "dd.hpp"

template <class R> struct Cx {
    R re, im;
    Cx() : re(0.0), im(0.0) {}
    Cx(R a) : re(a), im(0.0) {}
    Cx(R a, R b) : re(a), im(b) {}
};
template <class R> inline Cx<R> operator+(const Cx<R>& a, const Cx<R>& b) { return Cx<R>(a.re + b.re, a.im + b.im); }
template <class R> inline Cx<R> operator-(const Cx<R>& a, const Cx<R>& b) { return Cx<R>(a.re - b.re, a.im - b.im); }
template <class R> inline Cx<R> operator-(const Cx<R>& a) { return Cx<R>(-a.re, -a.im); }
template <class R> inline Cx<R> operator*(const Cx<R>& a, const Cx<R>& b) {
    return Cx<R>(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}
template <class R> inline Cx<R> operator*(const Cx<R>& a, const R& b) { return Cx<R>(a.re * b, a.im * b); }
inline Cx<dd> operator*(const Cx<dd>& a, double b) { return Cx<dd>(a.re * b, a.im * b); }
template <class R> inline Cx<R> operator/(const Cx<R>& a, const Cx<R>& b) {
    R d = b.re * b.re + b.im * b.im;
    return Cx<R>((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d);
}
template <class R> inline double cabs_d(const Cx<R>& a) { return std::hypot(to_double(a.re), to_double(a.im)); }

static void rd(std::ifstream& f, void* p, size_t n) { f.read((char*)p, n); if (!f) throw std::runtime_error("bad data file"); }

template <class R> struct System {
    typedef Cx<R> C;
    int n = 0, nl = 0;
    std::vector<int> letters;                           // letters with a non-zero matrix
    std::vector<std::vector<int>> ri, rj;               // per letter: sparse integer matrix (rows, cols, values)
    std::vector<std::vector<double>> rv;
    R x0, y0; C r0; std::vector<C> F0;
    C X, Y; double dlx = 0.1, dly = 0.2; long evals = 0;
    std::vector<C> tmp;

    explicit System(const char* fn) {
        std::ifstream f(fn, std::ios::binary);
        if (!f) throw std::runtime_error("cannot open data file");
        char mag[8]; rd(f, mag, 8);
        if (std::memcmp(mag, "ITERQ001", 8)) throw std::runtime_error("not an ITERQ001 file");
        int32_t h[2]; rd(f, h, sizeof h); n = h[0]; nl = h[1];
        ri.assign(nl + 1, {}); rj.assign(nl + 1, {}); rv.assign(nl + 1, {});
        for (int k = 1; k <= nl; ++k) {
            int32_t m; rd(f, &m, 4);
            ri[k].resize(m); rj[k].resize(m); rv[k].resize(m);
            if (m) { rd(f, ri[k].data(), 4 * m); rd(f, rj[k].data(), 4 * m); rd(f, rv[k].data(), 8 * m); letters.push_back(k); }
        }
        double p[4]; rd(f, p, sizeof p); x0 = R(p[0]) + R(p[1]); y0 = R(p[2]) + R(p[3]);
        double q[4]; rd(f, q, sizeof q); r0 = C(R(q[0]) + R(q[1]), R(q[2]) + R(q[3]));
        F0.resize(n);
        for (int i = 0; i < n; ++i) { double v[4]; rd(f, v, sizeof v); F0[i] = C(R(v[0]) + R(v[1]), R(v[2]) + R(v[3])); }
        tmp.resize(n);
    }

    // d/dt of the state (F_0..F_{n-1}, r)
    void rhs(const R& t, const std::vector<C>& s, std::vector<C>& ds) {
        ++evals;
        const C I(R(0.0), R(1.0));
        R one(1.0);
        C gx = C(t, R(4.0 * dlx) * t * (one - t)), gy = C(t, R(4.0 * dly) * t * (one - t));
        C dgx = C(one, R(4.0 * dlx) * (one - t * 2.0)), dgy = C(one, R(4.0 * dly) * (one - t * 2.0));
        C DX = X - C(x0), DY = Y - C(y0);
        C x = C(x0) + gx * DX, y = C(y0) + gy * DY, z = C(one) - x - y;
        C xp = dgx * DX, yp = dgy * DY, zp = -(xp + yp);
        C r = s[n];
        C b = x * y * z, bp = y * z * xp + x * z * yp + x * y * zp;
        std::vector<C> om(nl + 1);
        auto dlog = [&](const C& W, const C& Wx, const C& Wy) { return (Wx * xp + Wy * yp) / W; };
        C c1(one), c0(R(0.0)), c2(R(2.0));
        for (int k : letters) {
            switch (k) {
                case 1: om[k] = xp / x; break;
                case 2: om[k] = yp / y; break;
                case 3: om[k] = zp / z; break;
                case 4: om[k] = -xp / (c1 - x); break;
                case 5: om[k] = -yp / (c1 - y); break;
                case 6: om[k] = (xp + yp) / (x + y); break;
                case 7: case 8: {
                    C a = (k == 7) ? x * z : x * y, ap = (k == 7) ? xp * z + x * zp : xp * y + x * yp;
                    om[k] = I * (c2 * ap * b - a * bp) / (r * (a * a + b));
                    break;
                }
                case 9: om[k] = dlog(c1 - y - x * y, -y, -(c1 + x)); break;
                case 10: om[k] = dlog(c1 - x - x * y, -(c1 + y), -x); break;
                case 11: om[k] = dlog(x - y + y * y, c1, c2 * y - c1); break;
                case 12: om[k] = dlog(-x + x * x + y, c2 * x - c1, c1); break;
                case 13: om[k] = dlog((x - c1) * (x - c1) - y, c2 * (x - c1), -c1); break;
                case 14: om[k] = dlog((x + y) * (x + y) - y, c2 * (x + y), c2 * (x + y) - c1); break;
                case 15: om[k] = dlog((x + y) * (x + y) - x, c2 * (x + y) - c1, c2 * (x + y)); break;
                case 16: om[k] = dlog((c1 - y) * (c1 - y) - x, -c1, -c2 * (c1 - y)); break;
                case 17: om[k] = dlog(x + x * y + y * y, c1 + y, x + c2 * y); break;
                case 18: om[k] = dlog(x * x + y + x * y, c2 * x + y, c1 + x); break;
                case 19: om[k] = dlog((y - c1) * (y - c1) + x * y, y, c2 * (y - c1) + x); break;
                case 20: om[k] = dlog((x - c1) * (x - c1) + x * y, c2 * (x - c1) + y, x); break;
            }
        }
        for (int i = 0; i < n; ++i) ds[i] = c0;
        for (int k : letters) {
            const std::vector<int>& I_ = ri[k]; const std::vector<int>& J_ = rj[k]; const std::vector<double>& V_ = rv[k];
            for (int i = 0; i < n; ++i) tmp[i] = c0;
            for (size_t e = 0; e < I_.size(); ++e) tmp[I_[e]] = tmp[I_[e]] + s[J_[e]] * V_[e];
            const C w = om[k];
            for (int i = 0; i < n; ++i) ds[i] = ds[i] + w * tmp[i];
        }
        ds[n] = r * bp / (b * c2);
    }
};

template <class R> struct BS {
    typedef Cx<R> C;
    System<R>& S;
    double atol, rtol;
    int kmax;
    BS(System<R>& s, double tol, int km) : S(s), atol(tol), rtol(tol), kmax(km) {}

    void mmid(const std::vector<C>& y, const std::vector<C>& f0, const R& t, const R& H, int ns, std::vector<C>& out) {
        int N = (int)y.size();
        R h = H / R((double)ns);
        std::vector<C> ym = y, yn(N), f(N);
        for (int i = 0; i < N; ++i) yn[i] = y[i] + f0[i] * h;
        R x = t + h;
        S.rhs(x, yn, f);
        R h2 = h * 2.0;
        for (int m = 1; m < ns; ++m) {
            for (int i = 0; i < N; ++i) { C yt = ym[i] + f[i] * h2; ym[i] = yn[i]; yn[i] = yt; }
            x = x + h;
            S.rhs(x, yn, f);
        }
        out.resize(N);
        for (int i = 0; i < N; ++i) out[i] = (ym[i] + yn[i] + f[i] * h) * 0.5;
    }

    // integrate t = 0 .. 1; returns the number of accepted steps
    long run(std::vector<C>& y) {
        int N = (int)y.size();
        R t(0.0), one(1.0);
        double H = 0.05;
        long steps = 0;
        std::vector<C> f0(N);
        std::vector<int> nseq(kmax);
        for (int k = 0; k < kmax; ++k) nseq[k] = 2 * (k + 1);
        while (to_double(one - t) > 1e-30) {
            S.rhs(t, y, f0);
            bool ok = false;
            for (int tries = 0; tries < 60 && !ok; ++tries) {
                bool last = to_double(t) + H >= 1.0;
                R Hr = last ? one - t : R(H);
                double Hd = to_double(Hr);
                std::vector<std::vector<C>> prev, cur;
                for (int k = 0; k < kmax; ++k) {
                    cur.assign(k + 1, std::vector<C>());
                    mmid(y, f0, t, Hr, nseq[k], cur[0]);
                    for (int j = 1; j <= k; ++j) {                    // Neville in (h/n)^2
                        double a2 = double(nseq[k - j]) * nseq[k - j], b2 = double(nseq[k]) * nseq[k];
                        R fac = R(a2) / R(b2 - a2);                     // exact in R: a double factor limits dd to ~1e-24
                        cur[j].resize(N);
                        for (int i = 0; i < N; ++i) cur[j][i] = cur[j - 1][i] + (cur[j - 1][i] - prev[j - 1][i]) * fac;
                    }
                    if (k >= 2) {
                        double err = 0, ymax = 0;
                        for (int i = 0; i < N; ++i) ymax = std::max(ymax, cabs_d(y[i]));
                        for (int i = 0; i < N; ++i) {
                            double sc = atol * ymax + rtol * std::max(cabs_d(cur[k][i]), cabs_d(y[i]));   // relative to the largest function
                            err = std::max(err, cabs_d(cur[k][i] - cur[k - 1][i]) / sc);
                        }
                        if (err <= 1.0) {
                            y = cur[k]; t = last ? one : t + Hr; ++steps; ok = true;
                            double g = 0.94 * std::pow(0.65 / std::max(err, 1e-30), 1.0 / (2 * k + 1));
                            H = Hd * std::min(4.0, std::max(0.2, g));
                            if (k < kmax / 2) H *= 1.5;
                            break;
                        }
                    }
                    prev.swap(cur);
                }
                if (!ok) H = Hd * 0.25;
                if (H < 1e-12) throw std::runtime_error("step size underflow");
            }
            if (!ok) throw std::runtime_error("no convergence");
        }
        return steps;
    }
};

template <class R> int main_t(const char* fn, double tol, double dlx, double dly, int kmax) {
    typedef Cx<R> C;
    System<R> S(fn);
    S.dlx = dlx; S.dly = dly;
    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream is(line); std::vector<double> v; double a;
        while (is >> a) v.push_back(a);
        if (v.size() != 2 && v.size() != 4) continue;
        // "X Y", or "Xhi Xlo Yhi Ylo" for an end point in double-double (e.g. the image of a region point)
        double X = v[0], Y = v.size() == 4 ? v[2] : v[1];
        auto t0 = std::chrono::steady_clock::now();
        if (v.size() == 4) { S.X = C(R(v[0]) + R(v[1])); S.Y = C(R(v[2]) + R(v[3])); }
        else { S.X = C(R(X)); S.Y = C(R(Y)); }
        S.evals = 0;
        std::vector<C> y(S.n + 1);
        for (int i = 0; i < S.n; ++i) y[i] = S.F0[i];
        y[S.n] = S.r0;
        BS<R> bs(S, tol, kmax);
        long steps = bs.run(y);
        double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::printf("# %.17g %.17g %ld %ld %.6f %.17g %.17g\n", X, Y, steps, S.evals, dt,
                    to_double(y[S.n].re), to_double(y[S.n].im));          // the root r at the end point
        for (int i = 0; i < S.n; ++i) {
            dd re(y[i].re), im(y[i].im);
            std::printf("%.17g %.17g %.17g %.17g\n", re.hi, re.lo, im.hi, im.lo);
        }
        std::fflush(stdout);
    }
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) { std::cerr << "usage: iter_quad DATA.bin [quad|double] [tol] [dx dy] < points\n"; return 1; }
    std::string mode = argc > 2 ? argv[2] : "quad";
    double tol = argc > 3 ? std::atof(argv[3]) : (mode == "quad" ? 1e-28 : 1e-13);
    double dlx = argc > 5 ? std::atof(argv[4]) : 0.1, dly = argc > 5 ? std::atof(argv[5]) : 0.2;
    int kmax = argc > 6 ? std::atoi(argv[6]) : (mode == "quad" ? 16 : 10);
    try {
        if (mode == "quad") return main_t<dd>(argv[1], tol, dlx, dly, kmax);
        return main_t<double>(argv[1], tol, dlx, dly, kmax);
    } catch (std::exception& e) { std::cerr << "error: " << e.what() << "\n"; return 2; }
}
