#include "matrix.h"
#include <iostream>

std::vector<int> matrix::to_rref() {
    std::vector<int> piv;
    for (int j = 0; j < ncols; j++) {
        //this->print();
        int pivrow = -1;
        for (size_t i = piv.size(); i < rows.size(); i++) {
            if (rows[i][j]) {
                pivrow = i;
                break;
            }
        }
        //std::cout << "pivrow = " << pivrow << std::endl;
        if (pivrow == -1)
            continue;
        int i = piv.size();
        if (i != pivrow)
            std::swap(rows[i], rows[pivrow]);
        piv.push_back(j);
        pivrow = i;
        for (size_t i = 0; i < rows.size(); i++)
            if (static_cast<int>(i) != pivrow) {
                if (rows[i][j])
                    rows[i] += rows[pivrow];
            }
    }
    return piv;
}

void matrix::print() const {
    const auto &B = *this;
    for (int i = 0; i < static_cast<int>(rows.size()); i++) {
        for (int j = 0; j < ncols; j++)
            std::cout << B[i][j];
        std::cout << std::endl;
    }
}

std::vector<poly> matrix::rref_nullspace(const std::vector<int> &piv) const {
    std::vector<poly> ret;
    size_t k = 0;
    for (int j = 0; j < ncols; j++) {
        if (k < piv.size() && j == piv[k]) {
            k++;
            continue;
        }
        poly p(rows[0].bits());
        p.xor_bit(j, 1);
        for (size_t i = 0; i < piv.size(); i++)
            if (rows[i][j])
                p.xor_bit(piv[i], 1);
        ret.push_back(p);
    }
    return ret;
}
