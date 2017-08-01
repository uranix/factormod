#pragma once

#include <cstdint>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>

#include <cassert>

class poly {
    std::vector<uint32_t> x;
    mutable int _degree;
public:
    poly(size_t size) {
        assert((size % 32) == 0);
        x.resize(size / 32);
        _degree = 0;
    }

    int bits() const {
        return 32 * x.size();
    }

    poly(size_t size, const std::vector<int> &lsb) {
        x.resize(size / 32);
        for (size_t i = 0; i < lsb.size(); i++)
            if (lsb[i])
                x[i / 32] ^= 1 << (i % 32);
        _degree = -1;
    }

    int coeff(int i) const {
        assert(i < bits());
        return (x[i/32] & (1ul << (i%32))) ? 1 : 0;
    }

    void xor_coeff(int i, int v) {
        assert(i < bits());
        v = v ? 1 : 0;
        x[i/32] ^= v << (i%32);
        _degree = -1;
    }

    int degree() const {
        if (_degree >= 0)
            return _degree;
        _degree = bits();
        do {
            _degree--;
        } while ((_degree > 0) && (coeff(_degree) == 0));
        return _degree;
    }

    bool zero() const {
        return (degree() == 0) && (coeff(0) == 0);
    }

    std::string to_string() const {
        std::stringstream ss;

        for(size_t i = 0; i < x.size(); i++) {
            if (i > 0)
                ss << ' ';
            ss << std::setfill('0') << std::setw(8) << std::hex << x[i];
        }

        return ss.str();
    }

    std::string pretty() const {
        if (degree() == 0)
            return std::to_string(coeff(0));

        std::stringstream ss;
        for (int i = degree(); i >= 2; i--)
            if (coeff(i))
                ss << "x^" << i << " + ";
        if (coeff(1))
            ss << "x + ";
        if (coeff(0))
            ss << "1 + ";

        auto s = ss.str();
        return s.substr(0, s.size() - 3);
    }

    const poly operator+(const poly &o) const {
        assert(bits() == o.bits());
        poly sum(bits());
        for (size_t i = 0; i < x.size(); i++)
            sum.x[i] = x[i] ^ o.x[i];
        sum._degree = -1;
        return sum;
    }

    poly &operator+=(const poly &o) {
        assert(bits() == o.bits());
        for (size_t i = 0; i < x.size(); i++)
            x[i] ^= o.x[i];
        _degree = -1;
        return *this;
    }

    const poly operator-(const poly &o) const {
        return (*this) + o;
    }

    poly &operator-=(const poly &o) {
        return (*this) += o;
    }

    bool operator==(const poly &o) const {
        assert(bits() == o.bits());
        for (size_t i = 0; i < x.size(); i++)
            if (x[i] != o.x[i])
                return false;
        return true;
    }

    const poly operator*(const poly &o) const {
        if (degree() + o.degree() >= bits())
            throw std::overflow_error(
                    "Attempt to multiply deg(f) = " + std::to_string(degree()) +
                    " with deg(g) = " + std::to_string(o.degree()) + " while using " +
                    std::to_string(bits()) + " bits");

        poly ans(bits());
        for (int i = 0; i <= degree(); i++)
            for (int j = 0; j <= o.degree(); j++)
                ans.xor_coeff(i+j, coeff(i) * o.coeff(j));
        return ans;
    }

    const poly derivative() const {
        poly ans(bits());
        for (int i = 1; i < bits(); i += 2)
            ans.xor_coeff(i-1, coeff(i));
        return ans;
    }

    struct div;

    div operator/(const poly &f) const;
};

struct poly::div {
    poly quot;
    poly rem;
    div(int bits) : quot(bits), rem(bits) { }
};

inline poly::div poly::operator/(const poly &f) const {
    div ans(f.bits());
    ans.rem = *this;
    if (f.degree() == 0) {
        if (f.coeff(0) == 0)
            throw std::runtime_error("Division by zero");
        ans.rem = poly(f.bits());
        ans.quot = *this;
        return ans;
    }

    auto &r = ans.rem;
    auto &q = ans.quot;
    for (int i = r.degree(); i >= f.degree(); i--) {
        if (r.coeff(i)) {
            int k = i - f.degree();
            q.xor_coeff(k, 1);
            for (int j = 0; j <= f.degree(); j++)
                r.xor_coeff(j + k, f.coeff(j));
        }
    }

    return ans;
}

inline const poly gcd(const poly &f, const poly &g) {
    poly a = f.degree() > g.degree() ? f : g;
    poly b = f.degree() > g.degree() ? g : f;

    do {
        if (b.zero())
            return a;
        auto c = a / b;
        a = b;
        b = c.rem;
    } while (true);
}
