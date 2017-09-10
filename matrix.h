#pragma once

#include "poly.h"
struct matrix {
    std::vector<poly> rows;

    int ncols;
    matrix(int m, int n, int lda) : ncols(n) {
        for (int i = 0; i < m; i++)
            rows.emplace_back(lda);
    }
    poly &operator[](int i) {
        return rows[i];
    }
    const poly &operator[](int i) const {
        return rows[i];
    }
    void sub_unit_diag() {
        for (int i = 0; i < static_cast<int>(rows.size()); i++)
            rows[i].xor_bit(i, 1);
    }
    std::vector<int> to_rref();
    std::vector<poly> rref_nullspace(const std::vector<int> &piv) const;

    void print() const;
};
