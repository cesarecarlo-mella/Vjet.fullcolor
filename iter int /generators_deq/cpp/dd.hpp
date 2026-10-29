// dd.hpp -- double-double arithmetic (about 32 significant digits), after Hida, Li, Bailey (QD).
// Error-free transformations: compile WITHOUT -ffast-math / -fassociative-math, and with -ffp-contract=off.
#pragma once
#include <cmath>
#include <cstdio>

struct dd {
    double hi, lo;
    dd() : hi(0), lo(0) {}
    dd(double h) : hi(h), lo(0) {}
    dd(double h, double l) : hi(h), lo(l) {}
};

inline dd quick_two_sum(double a, double b) { double s = a + b; return dd(s, b - (s - a)); }
inline dd two_sum(double a, double b) { double s = a + b, bb = s - a; return dd(s, (a - (s - bb)) + (b - bb)); }
inline dd two_prod(double a, double b) { double p = a * b; return dd(p, std::fma(a, b, -p)); }

inline dd operator+(const dd& a, const dd& b) {
    dd s = two_sum(a.hi, b.hi), t = two_sum(a.lo, b.lo);
    s.lo += t.hi; s = quick_two_sum(s.hi, s.lo); s.lo += t.lo;
    return quick_two_sum(s.hi, s.lo);
}
inline dd operator-(const dd& a) { return dd(-a.hi, -a.lo); }
inline dd operator-(const dd& a, const dd& b) { return a + (-b); }
inline dd operator*(const dd& a, const dd& b) {
    dd p = two_prod(a.hi, b.hi);
    p.lo += a.hi * b.lo + a.lo * b.hi;
    return quick_two_sum(p.hi, p.lo);
}
inline dd operator*(const dd& a, double b) { dd p = two_prod(a.hi, b); p.lo += a.lo * b; return quick_two_sum(p.hi, p.lo); }
inline dd operator*(double b, const dd& a) { return a * b; }
inline dd operator/(const dd& a, const dd& b) {
    double q1 = a.hi / b.hi; dd r = a - b * q1;
    double q2 = r.hi / b.hi; r = r - b * q2;
    double q3 = r.hi / b.hi;
    dd q = quick_two_sum(q1, q2);
    return q + dd(q3);
}
inline dd operator/(const dd& a, double b) { return a / dd(b); }
inline dd& operator+=(dd& a, const dd& b) { a = a + b; return a; }
inline dd& operator-=(dd& a, const dd& b) { a = a - b; return a; }
inline dd& operator*=(dd& a, const dd& b) { a = a * b; return a; }
inline bool operator<(const dd& a, const dd& b) { return a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo); }
inline dd sqrt(const dd& a) {
    if (a.hi <= 0) return dd(0);
    double x = 1.0 / std::sqrt(a.hi);
    double ax = a.hi * x;
    return two_sum(ax, (a - two_prod(ax, ax)).hi * (x * 0.5));
}
inline double to_double(const dd& a) { return a.hi + a.lo; }
inline double to_double(double a) { return a; }
inline double fabs_d(const dd& a) { return std::fabs(a.hi); }
inline double fabs_d(double a) { return std::fabs(a); }
