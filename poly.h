#pragma once

#include <cstdint>
#include <string>
#include <istream>
#include <stdexcept>
#include <vector>

#include <cassert>

struct matrix;

class poly {
    std::vector<uint32_t> x;
    mutable int _degree;
public:
    int bits() const {
        return 32 * x.size();
    }

    poly(size_t bits) {
        assert((bits % 32) == 0);
        x.resize(bits / 32);
        _degree = 0;
    }

    poly(size_t bits, const std::vector<int> &lsb) {
        x.resize(bits / 32);
        for (size_t i = 0; i < lsb.size(); i++)
            if (lsb[i])
                x[i / 32] ^= 1 << (i % 32);
        _degree = -1;
    }

    poly(size_t bits, std::istream &ss);

    int operator[](int i) const {
        assert(i < bits());
        return (x[i/32] & (1ul << (i%32))) ? 1 : 0;
    }

    void xor_bit(int i, int v) {
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
        } while ((_degree > 0) && ((*this)[_degree] == 0));
        return _degree;
    }

    bool is_zero() const {
        return (degree() == 0) && ((*this)[0] == 0);
    }

    bool is_unit() const {
        return (degree() == 0) && ((*this)[0] == 1);
    }

    bool operator<(const poly &b) const {
        const auto &a = *this;
        if (a.degree() < b.degree())
            return true;
        if (a.degree() > b.degree())
            return false;
        for (int i = a.degree(); i >= 0; i--) {
            if (a[i] < b[i])
                return true;
            if (a[i] > b[i])
                return false;
        }
        return false;
    }

    std::string to_string() const;
    std::string pretty() const;

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
                ans.xor_bit(i+j, (*this)[i] * o[j]);
        return ans;
    }

    const poly derivative() const {
        poly ans(bits());
        for (int i = 1; i < bits(); i += 2)
            ans.xor_bit(i-1, (*this)[i]);
        return ans;
    }

    struct div;

    div operator/(const poly &f) const;

    /// gcd(f, f')
    const poly double_factor() const;

    /// x^{2i} = sum_j b_{ji} x^j (mod f)
    const matrix powers_mod() const;
};

struct poly::div {
    poly quot;
    poly rem;
    div(int bits) : quot(bits), rem(bits) { }
};

const poly gcd(const poly &f, const poly &g);
