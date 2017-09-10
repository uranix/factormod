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

void rec(const std::vector<poly> &fact, const std::vector<int> &pows, const int bits,
    std::vector<int> &p, std::vector<std::string> &ans, size_t k = 0)
{
    if (k == pows.size()) {
        int d1 = 0;
        int d2 = 0;
        for (size_t i = 0; i < k; i++) {
            int fd = fact[i].degree();
            d1 += p[i] * fd;
            d2 += (pows[i] - p[i]) * fd;
        }
        if (d1 < bits && d2 < bits) {
            poly p1(bits, {1});
            poly p2(bits, {1});
            for (size_t i = 0; i < k; i++) {
                int j = 0;
                for (; j < p[i]; j++)
                    p1 = p1 * fact[i];
                for (; j < pows[i]; j++)
                    p2 = p2 * fact[i];
            }
            ans.push_back(p1.to_string() + " " + p2.to_string());
        }
        return;
    }
    for (p[k] = 0; p[k] <= pows[k]; p[k]++)
        rec(fact, pows, bits, p, ans, k+1);
}

int main() {

    int bits;
    std::cin >> bits;
    poly f(2*bits, std::cin);

    //std::cout << f.pretty() << std::endl;

    auto ff = factor(f);

    std::vector<poly> factors;
    std::vector<int> powers;

    for (const auto &p : ff) {
        if (factors.empty() || !(factors.back() == p))
            factors.push_back(p);
    }
    for (size_t i = 0; i < factors.size(); i++)
        powers.push_back(ff.count(factors[i]));

    std::vector<int> pows(factors.size());
    std::vector<std::string> ans;
    rec(factors, powers, bits, pows, ans);
    std::sort(ans.begin(), ans.end());

    for (const auto &s : ans)
        std::cout << s << std::endl;

    return 0;
}
