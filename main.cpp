#include <iostream>
#include "poly.h"

void test() {
    auto f = poly(32, {1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1});
    auto x2 = poly(32, {0, 0, 1});
    auto q = poly(32, {1});
    std::cout << "f = " << f.pretty() << std::endl;
    auto fp = f.derivative();
    std::cout << "f' = " << fp.pretty() << std::endl;
    std::cout << "gcd(f, f') = " << gcd(f, fp).pretty() << std::endl;
    for (int i = 0; i < f.degree(); i++) {
        std::cout << "x^" << (2*i) << " === " << q.pretty() << " (mod f)" << std::endl;
        q = q * x2;
        auto c = q / f;
        q = c.rem;
    }
}

int main() {
    poly f(32, {1, 0, 1});
    auto g = f*f;
    auto h = g*f;
    auto h2 = f*g;
    auto u = h + h2;

    std::cout << "f = " << f.pretty() << std::endl;
    std::cout << "g = f*f = " << g.pretty() << std::endl;
    std::cout << "h = f*g = " << h.pretty() << std::endl;
    std::cout << "h2 = g*f = " << h2.pretty() << std::endl;
    std::cout << "u = h + h2 = " << u.pretty() << std::endl;

    f = poly(32, {1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1});
    g = poly(32, {1, 1, 0, 1, 1, 0, 0, 0, 1});
    std::cout << "f = " << f.pretty() << std::endl;
    std::cout << "g = " << g.pretty() << std::endl;
    auto qr = f / g;
    std::cout << "q = " << qr.quot.pretty() << std::endl;
    std::cout << "r = " << qr.rem.pretty() << std::endl;
    auto rf = qr.quot*g+qr.rem;
    std::cout << "q*g + r = " << rf.pretty() << std::endl;
    std::cout << "gcd(h*h*h, h*f) = " << gcd(h*h*h, h*f).pretty() << std::endl;
    std::cout << "f' = " << f.derivative().pretty() << std::endl;

    test();

    return 0;
}
