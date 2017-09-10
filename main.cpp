#include <iostream>
#include "poly.h"
#include "matrix.h"

#include <sstream>
#include <deque>
#include <set>
#include <algorithm>

// Berlekamp algo, square free case
std::set<poly> factor_sq_free(const poly &f) {
    assert(f.double_factor().is_unit());

    auto B = f.powers_mod();
    B.sub_unit_diag();
    const auto &piv = B.to_rref();
    const auto &hv = B.rref_nullspace(piv);

    std::deque<poly> hq;
    for (const auto &p : hv)
        hq.push_back(p);
    assert(hq.front().is_unit());
    hq.pop_front();

    std::set<poly> in, out;
    in.insert(f);

    while (!hq.empty()) {
        const poly &h0 = hq.front();
        poly h1(h0);
        h1.xor_bit(0, 1);
        out.clear();
        for (const auto &p : in) {
            const poly d0 = gcd(p, h0);
            const poly d1 = gcd(p, h1);
            if (!d0.is_unit())
                out.insert(d0);
            if (!d1.is_unit())
                out.insert(d1);
        }
        hq.pop_front();
        std::swap(in, out);
    }
    return in;
}

// generic case
std::multiset<poly> factor(const poly &f) {
    const poly d = f.double_factor();
    if (d.is_unit()) {
        const auto &facs = factor_sq_free(f);
        std::multiset<poly> ret;
        for (const auto &f : facs)
            ret.insert(f);
        return ret;
    }
    if (d == f) {
        // f' = 0, f only has even powers
        poly g(f.bits());
        for (int i = 0; i <= f.degree(); i += 2)
            g.xor_bit(i / 2, f[i]);
        std::multiset<poly> gf = factor(g);
        std::multiset<poly> ret;
        for (const auto &p : gf) {
            ret.insert(p);
            ret.insert(p);
        }
        return ret;
    }
    std::multiset<poly> ret = factor(d);
    std::set<poly> add = factor_sq_free((f / d).quot);
    for (const auto &p : add)
        ret.insert(p);
    return ret;
}

int main() {

    int bits;
    std::cin >> bits;
    poly f(2*bits, std::cin);

    std::cout << "Factors of " << f.pretty() << ":" << std::endl;

    auto ff = factor(f);

    std::vector<poly> factors;
    std::vector<int> powers;

    for (const auto &p : ff) {
        if (factors.empty() || !(factors.back() == p))
            factors.push_back(p);
    }
    for (size_t i = 0; i < factors.size(); i++)
        powers.push_back(ff.count(factors[i]));

    for (size_t i = 0; i < factors.size(); i++)
        std::cout << "(" << factors[i].pretty() << ") ^ " << powers[i] << std::endl;

    return 0;
}
