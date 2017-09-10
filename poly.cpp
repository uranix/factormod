#include "poly.h"
#include "matrix.h"

#include <iomanip>
#include <sstream>

std::string poly::to_string() const {
    std::stringstream ss;

    for(size_t i = 0; i < x.size(); i++) {
        if (i > 0)
            ss << ' ';
        ss << std::setfill('0') << std::setw(8) << std::hex << x[i];
    }

    return ss.str();
}

poly::poly(size_t bits, std::istream &ss) {
    x.resize(bits / 32);
    for (size_t i = 0; i < bits / 32; i++)
        ss >> std::hex >> x[i];
    _degree = -1;
}

std::string poly::pretty() const {
    const auto &rthis = *this;
    if (degree() == 0)
        return std::to_string(rthis[0]);

    std::stringstream ss;
    for (int i = degree(); i >= 2; i--)
        if (rthis[i])
            ss << "x^" << i << " + ";
    if (rthis[1])
        ss << "x + ";
    if (rthis[0])
        ss << "1 + ";

    auto s = ss.str();
    return s.substr(0, s.size() - 3);
}

poly::div poly::operator/(const poly &f) const {
    div ans(f.bits());
    ans.rem = *this;
    if (f.degree() == 0) {
        if (f[0] == 0)
            throw std::runtime_error("Division by is_zero");
        ans.rem = poly(f.bits());
        ans.quot = *this;
        return ans;
    }

    auto &r = ans.rem;
    auto &q = ans.quot;
    for (int i = r.degree(); i >= f.degree(); i--) {
        if (r[i]) {
            int k = i - f.degree();
            q.xor_bit(k, 1);
            for (int j = 0; j <= f.degree(); j++)
                r.xor_bit(j + k, f[j]);
        }
    }

    return ans;
}

const poly gcd(const poly &f, const poly &g) {
    poly a = f.degree() > g.degree() ? f : g;
    poly b = f.degree() > g.degree() ? g : f;

    do {
        if (b.is_zero())
            return a;
        auto c = a / b;
        a = b;
        b = c.rem;
    } while (true);
}

const poly poly::double_factor() const {
    auto fp = this->derivative();
    return gcd(*this, fp);
}

const matrix poly::powers_mod() const {
    const auto &f = *this;
    matrix B(f.degree(), f.degree(), f.bits());
    poly q(f.bits(), {1});
    poly x2(f.bits(), {0, 0, 1});
    for (int i = 0; i < f.degree(); i++) {
        for (int j = 0; j < f.degree(); j++)
            B[j].xor_bit(i, q[j]);
        q = q * x2;
        auto c = q / f;
        q = c.rem;
    }
    return B;
}
