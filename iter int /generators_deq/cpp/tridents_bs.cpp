// tridents_bs.cpp -- the 1318 trident functions (eps^0..eps^6) in double precision with the method of
// Czakon & Tancredi (arXiv:2606.30354): ONE first-order system containing all eps orders and the square root
// r = Sqrt[x y z] (dr = r db / 2b), integrated with boost::odeint's Bulirsch-Stoer along the complex-deformed
// straight path z(t) = P0 + (t + 4 i delta t (1-t)) (P - P0).
//
// The system is the redundant weight<=5 set W of 600 family elements (notes/redundant_W.tex), written in an
// eps-graded orthonormal basis:  W_w = Q_w y_w  with dim y_w = d_w = 1, 3, 10, 32, 105, 330 (w = 0..5):
//     dy_w/dt = sum_k f_k(t) A[k,w] y_{w-1}                 (A[k,w]: d_w x d_{w-1}, dense or low rank)
//     dJ[k,w]/dt = f_k(t) y_w                               (per-letter integrals, w = 0..5)
//     dr/dt = r b'/(2 b)
// and the tridents are read off at the end point:  T_w = T_w(P0) + sum_k R[k,w] J[k,w-1].
//
//   ./tridents_bs DATA.bin [eps_abs eps_rel delta_x delta_y] < points > values   (defaults 1e-10 1e-10 0.1 0.2)
// points: "x y" per line;  output per point: "# x y steps evals time_s" then 7 x 1318 lines "re im"
// (order w = 0..6, trident 1..1318).  Set QUIET=1 to print only the "#" line.
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include <complex>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <chrono>
#include <cstring>
#include <cstdlib>
#include <tuple>

using cd = std::complex<double>;
using state = std::vector<cd>;
namespace ode = boost::numeric::odeint;

struct Block {                      // out += f * M y,  M = U V (low rank) or dense; real/imag stored separately
    int type = 0, m = 0, n = 0, r = 0;
    std::vector<double> ar, ai, ur, ui, vr, vi;     // row-major
    mutable std::vector<double> tr, ti, yr, yi;
    static inline void mv(const double* Ar, const double* Ai, int m, int n, const double* xr, const double* xi,
                          double* outr, double* outi) {          // out = A x  (complex, split storage), 4 rows at a time
        int i = 0;
        for (; i + 4 <= m; i += 4) {
            const double *r0 = Ar + (size_t)i * n, *r1 = r0 + n, *r2 = r1 + n, *r3 = r2 + n;
            const double *i0 = Ai + (size_t)i * n, *i1 = i0 + n, *i2 = i1 + n, *i3 = i2 + n;
            double s0r = 0, s0i = 0, s1r = 0, s1i = 0, s2r = 0, s2i = 0, s3r = 0, s3i = 0;
            for (int j = 0; j < n; ++j) {
                const double a = xr[j], b = xi[j];
                s0r += r0[j] * a - i0[j] * b; s0i += r0[j] * b + i0[j] * a;
                s1r += r1[j] * a - i1[j] * b; s1i += r1[j] * b + i1[j] * a;
                s2r += r2[j] * a - i2[j] * b; s2i += r2[j] * b + i2[j] * a;
                s3r += r3[j] * a - i3[j] * b; s3i += r3[j] * b + i3[j] * a;
            }
            outr[i] = s0r; outi[i] = s0i; outr[i + 1] = s1r; outi[i + 1] = s1i;
            outr[i + 2] = s2r; outi[i + 2] = s2i; outr[i + 3] = s3r; outi[i + 3] = s3i;
        }
        for (; i < m; ++i) {
            const double* pr = Ar + (size_t)i * n; const double* pi = Ai + (size_t)i * n;
            double sr = 0, si = 0;
            for (int j = 0; j < n; ++j) { sr += pr[j] * xr[j] - pi[j] * xi[j]; si += pr[j] * xi[j] + pi[j] * xr[j]; }
            outr[i] = sr; outi[i] = si;
        }
    }
    void apply(const cd f, const cd* y, cd* out) const {
        if (type == 0) return;
        yr.resize(n); yi.resize(n);
        for (int j = 0; j < n; ++j) { yr[j] = y[j].real(); yi[j] = y[j].imag(); }
        tr.resize(m); ti.resize(m);
        if (type == 1) mv(ar.data(), ai.data(), m, n, yr.data(), yi.data(), tr.data(), ti.data());
        else {
            std::vector<double> qr(r), qi(r);
            mv(vr.data(), vi.data(), r, n, yr.data(), yi.data(), qr.data(), qi.data());
            mv(ur.data(), ui.data(), m, r, qr.data(), qi.data(), tr.data(), ti.data());
        }
        for (int i = 0; i < m; ++i) out[i] += f * cd(tr[i], ti[i]);
    }
    // projection V y for low-rank read-out blocks
    void project(const cd f, const cd* y, cd* out) const {
        yr.resize(n); yi.resize(n);
        for (int j = 0; j < n; ++j) { yr[j] = y[j].real(); yi[j] = y[j].imag(); }
        std::vector<double> qr(r), qi(r);
        mv(vr.data(), vi.data(), r, n, yr.data(), yi.data(), qr.data(), qi.data());
        for (int i = 0; i < r; ++i) out[i] = f * cd(qr[i], qi[i]);
    }
};

static void rd(std::ifstream& f, void* p, size_t n) { f.read((char*)p, n); if (!f) throw std::runtime_error("bad data file"); }
static Block read_block(std::ifstream& f) {
    Block b; int h[4]; rd(f, h, sizeof h); b.type = h[0]; b.m = h[1]; b.n = h[2]; b.r = h[3];
    auto get = [&](size_t cnt, std::vector<double>& re, std::vector<double>& im) {
        std::vector<cd> t(cnt); rd(f, t.data(), cnt * sizeof(cd)); re.resize(cnt); im.resize(cnt);
        for (size_t i = 0; i < cnt; ++i) { re[i] = t[i].real(); im[i] = t[i].imag(); } };
    if (b.type == 1) get((size_t)b.m * b.n, b.ar, b.ai);
    if (b.type == 2) { get((size_t)b.m * b.r, b.ur, b.ui); get((size_t)b.r * b.n, b.vr, b.vi); }
    return b;
}

// plain letters as polynomials: {i, j, c} for x^i y^j
struct Mono { int i, j; double c; };
static const std::vector<std::pair<int, std::vector<Mono>>> PLAIN = {
    {1, {{1,0,1}}}, {2, {{0,1,1}}}, {3, {{0,0,1},{1,0,-1},{0,1,-1}}}, {4, {{0,0,1},{1,0,-1}}}, {5, {{0,0,1},{0,1,-1}}},
    {6, {{1,0,1},{0,1,1}}}, {9, {{0,0,1},{0,1,-1},{1,1,-1}}}, {10, {{0,0,1},{1,0,-1},{1,1,-1}}},
    {11, {{1,0,1},{0,1,-1},{0,2,1}}}, {12, {{1,0,-1},{2,0,1},{0,1,1}}}, {13, {{0,0,1},{1,0,-2},{2,0,1},{0,1,-1}}},
    {14, {{0,1,-1},{2,0,1},{1,1,2},{0,2,1}}}, {15, {{1,0,-1},{2,0,1},{1,1,2},{0,2,1}}},
    {16, {{1,0,-1},{0,0,1},{0,1,-2},{0,2,1}}}, {17, {{1,0,1},{1,1,1},{0,2,1}}}, {18, {{2,0,1},{0,1,1},{1,1,1}}},
    {19, {{0,0,1},{0,1,-2},{0,2,1},{1,1,1}}}, {20, {{0,0,1},{1,0,-2},{2,0,1},{1,1,1}}}};

struct System {
    int NL, NW, NT; std::vector<int> d;                // d[w], w = 0..5
    double x0, y0; cd r0;
    std::vector<std::vector<cd>> y0v; std::vector<cd> T0;       // y0v[w], T0[(w)*NT + i]
    std::vector<std::vector<Block>> A, R;                         // A[k][w-1] (w=1..5), R[k][w-1] (w=1..6)
    std::vector<int> offy, offJ, dimJ; int ir, N;                       // state layout
    // path
    double Dx = 0, Dy = 0, dlx = 0.1, dly = 0.2;
    mutable long evals = 0;

    explicit System(const char* fn) {
        std::ifstream f(fn, std::ios::binary); char mag[8]; rd(f, mag, 8);
        if (std::memcmp(mag, "TRIDBS01", 8)) throw std::runtime_error("not a TRIDBS01 file");
        int h[3]; rd(f, h, sizeof h); NL = h[0]; NW = h[1]; NT = h[2];
        d.resize(NW); rd(f, d.data(), NW * sizeof(int));
        double p[2]; rd(f, p, sizeof p); x0 = p[0]; y0 = p[1]; rd(f, &r0, sizeof r0);
        y0v.resize(NW); for (int w = 0; w < NW; ++w) { y0v[w].resize(d[w]); rd(f, y0v[w].data(), d[w] * sizeof(cd)); }
        T0.resize(7 * NT); rd(f, T0.data(), T0.size() * sizeof(cd));
        A.assign(NL, {}); R.assign(NL, {});
        for (int k = 0; k < NL; ++k) for (int w = 1; w <= 5; ++w) A[k].push_back(read_block(f));
        for (int k = 0; k < NL; ++k) for (int w = 1; w <= 6; ++w) R[k].push_back(read_block(f));
        // layout: y_1..y_5, J[k][w] (w = 0..5), r
        int o = 0; offy.assign(NW, -1);
        for (int w = 1; w <= 5; ++w) { offy[w] = o; o += d[w]; }
        offJ.assign(NL * 6, 0); dimJ.assign(NL * 6, 0);
        for (int k = 0; k < NL; ++k) for (int w = 0; w <= 5; ++w) {
            const Block& B = R[k][w];                     // read-out block of T_{w+1}: acts on J[k][w]
            dimJ[k * 6 + w] = B.type == 0 ? 0 : (B.type == 2 ? B.r : d[w]);
            offJ[k * 6 + w] = o; o += dimJ[k * 6 + w];
        }
        ir = o; N = o + 1;
        build_stacks();
    }

    // ---- stacked matrices: one contiguous matrix-vector product per order ----
    struct Stack {
        int rows = 0, n = 0; std::vector<double> ar, ai;
        std::vector<std::tuple<int, int, int, const Block*>> parts;   // (letter, row offset, rows, block for U or null)
        mutable std::vector<double> tr, ti, xr, xi;
        void add(int k, const std::vector<double>& re, const std::vector<double>& im, int r, int nn, const Block* U) {
            n = nn; parts.emplace_back(k, rows, r, U); rows += r;
            ar.insert(ar.end(), re.begin(), re.end()); ai.insert(ai.end(), im.begin(), im.end());
        }
        void mul(const cd* y) const {
            xr.resize(n); xi.resize(n); tr.resize(rows); ti.resize(rows);
            for (int j = 0; j < n; ++j) { xr[j] = y[j].real(); xi[j] = y[j].imag(); }
            Block::mv(ar.data(), ai.data(), rows, n, xr.data(), xi.data(), tr.data(), ti.data());
        }
    };
    std::vector<Stack> SA, SJ;                     // SA[w] (w = 1..5) for dy_w; SJ[w] (w = 0..5) for J projections
    void build_stacks() {
        SA.assign(6, {}); SJ.assign(6, {});
        for (int w = 1; w <= 5; ++w)
            for (int k = 0; k < NL; ++k) {
                const Block& B = A[k][w - 1];
                if (B.type == 1) SA[w].add(k, B.ar, B.ai, B.m, B.n, nullptr);
                else if (B.type == 2) SA[w].add(k, B.vr, B.vi, B.r, B.n, &B);
            }
        for (int w = 0; w <= 5; ++w)
            for (int k = 0; k < NL; ++k) { const Block& B = R[k][w]; if (B.type == 2) SJ[w].add(k, B.vr, B.vi, B.r, B.n, nullptr); }
    }

    void letters(double t, cd r, cd* fk, cd& dr_over_r) const {
        const cd I(0, 1);
        const cd gx = t + 4.0 * I * dlx * t * (1 - t), gy = t + 4.0 * I * dly * t * (1 - t);
        const cd x = x0 + gx * Dx, y = y0 + gy * Dy;
        const cd dx = (1.0 + 4.0 * I * dlx * (1 - 2 * t)) * Dx, dy = (1.0 + 4.0 * I * dly * (1 - 2 * t)) * Dy;
        cd xp[3] = {1, x, x * x}, yp[3] = {1, y, y * y};
        for (auto& [k, P] : PLAIN) {
            cd v = 0, dv = 0;
            for (auto& m : P) {
                v += m.c * xp[m.i] * yp[m.j];
                if (m.i) dv += m.c * double(m.i) * xp[m.i - 1] * yp[m.j] * dx;
                if (m.j) dv += m.c * double(m.j) * xp[m.i] * yp[m.j - 1] * dy;
            }
            fk[k - 1] = dv / v;
        }
        const cd z = 1.0 - x - y, dz = -dx - dy;
        const cd b = x * y * z, db = dx * y * z + x * dy * z + x * y * dz;
        for (int k = 7; k <= 8; ++k) {
            cd a, da;
            if (k == 7) { a = x * z; da = dx * z + x * dz; } else { a = x * y; da = dx * y + x * dy; }
            fk[k - 1] = I * (2.0 * b * da - a * db) / (r * (a * a + b));
        }
        dr_over_r = db / (2.0 * b);
    }

    void operator()(const state& s, state& ds, double t) const {
        ++evals;
        cd fk[20], drr;
        letters(t, s[ir], fk, drr);
        std::fill(ds.begin(), ds.end(), cd(0));
        std::vector<double> qr, qi, ur, ui;
        for (int w = 1; w <= 5; ++w) {
            const cd* yprev = (w == 1) ? y0v[0].data() : &s[offy[w - 1]];
            cd* out = &ds[offy[w]];
            const Stack& St = SA[w]; St.mul(yprev);
            for (auto& [k, off, r, U] : St.parts) {
                const cd f = fk[k];
                if (!U) for (int i = 0; i < r; ++i) out[i] += f * cd(St.tr[off + i], St.ti[off + i]);
                else {                                   // low rank: out += f U (V y)
                    ur.resize(U->m); ui.resize(U->m);
                    Block::mv(U->ur.data(), U->ui.data(), U->m, r, &St.tr[off], &St.ti[off], ur.data(), ui.data());
                    for (int i = 0; i < U->m; ++i) out[i] += f * cd(ur[i], ui[i]);
                }
            }
        }
        for (int w = 0; w <= 5; ++w) {
            const cd* yw = (w == 0) ? y0v[0].data() : &s[offy[w]];
            const Stack& St = SJ[w]; if (St.rows) St.mul(yw);
            for (int k = 0; k < NL; ++k) {
                const int dj = dimJ[k * 6 + w]; if (!dj) continue;
                if (R[k][w].type == 1) { cd* out = &ds[offJ[k * 6 + w]]; const cd f = fk[k]; for (int i = 0; i < dj; ++i) out[i] = f * yw[i]; }
            }
            for (auto& [k, off, r, U] : St.parts) {
                cd* out = &ds[offJ[k * 6 + w]]; const cd f = fk[k];
                for (int i = 0; i < r; ++i) out[i] = f * cd(St.tr[off + i], St.ti[off + i]);
            }
        }
        ds[ir] = s[ir] * drr;
    }

    state initial() const {
        state s(N, cd(0));
        for (int w = 1; w <= 5; ++w) std::copy(y0v[w].begin(), y0v[w].end(), s.begin() + offy[w]);
        s[ir] = r0; return s;
    }

    std::vector<cd> readout(const state& s) const {
        std::vector<cd> T(T0);
        for (int w = 1; w <= 6; ++w)
            for (int k = 0; k < NL; ++k) {
                const Block& B = R[k][w - 1]; const cd* J = &s[offJ[k * 6 + (w - 1)]]; cd* out = &T[w * NT];
                if (B.type == 1) B.apply(cd(1), J, out);
                else if (B.type == 2) {               // out += U J
                    std::vector<double> jr(B.r), ji(B.r), tr(B.m), ti(B.m);
                    for (int j = 0; j < B.r; ++j) { jr[j] = J[j].real(); ji[j] = J[j].imag(); }
                    Block::mv(B.ur.data(), B.ui.data(), B.m, B.r, jr.data(), ji.data(), tr.data(), ti.data());
                    for (int i = 0; i < B.m; ++i) out[i] += cd(tr[i], ti[i]);
                }
            }
        return T;
    }
};

int main(int argc, char** argv) {
    if (argc < 2) { std::cerr << "usage: tridents_bs DATA.bin [eps_abs eps_rel delta_x delta_y] < points\n"; return 1; }
    System S(argv[1]);
    const double eabs = argc > 2 ? std::atof(argv[2]) : 1e-10, erel = argc > 3 ? std::atof(argv[3]) : 1e-10;
    if (argc > 5) { S.dlx = std::atof(argv[4]); S.dly = std::atof(argv[5]); }
    const bool quiet = std::getenv("QUIET") != nullptr;
    std::cerr << "state size " << S.N << " complex components\n";
    if (std::getenv("BENCH")) {
        S.Dx = 0.74; S.Dy = 0.05; state s = S.initial(), ds(S.N);
        auto t0 = std::chrono::steady_clock::now();
        for (int i = 0; i < 2000; ++i) S(s, ds, 0.37);
        double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        long mac = 0; for (auto& St : S.SA) mac += (long)St.rows * St.n; for (auto& St : S.SJ) mac += (long)St.rows * St.n;
        std::cerr << "RHS: " << dt / 2000 * 1e3 << " ms per evaluation; stacked MACs " << mac << "\n"; return 0;
    }
    std::string line;
    std::cout.precision(16);
    while (std::getline(std::cin, line)) {
        std::istringstream is(line); double x, y;
        if (!(is >> x >> y)) continue;
        auto t0 = std::chrono::steady_clock::now();
        S.Dx = x - S.x0; S.Dy = y - S.y0; S.evals = 0;
        state s = S.initial();
        ode::bulirsch_stoer<state, double, state, double> stepper(eabs, erel, 1.0, 0.0);
        size_t steps = ode::integrate_adaptive(stepper, std::ref(S), s, 0.0, 1.0, 0.05);
        auto T = S.readout(s);
        double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::cout << "# " << x << " " << y << " " << steps << " " << S.evals << " " << dt << "\n";
        if (!quiet) for (auto& v : T) std::cout << v.real() << " " << v.imag() << "\n";
        std::cout.flush();
    }
    return 0;
}
