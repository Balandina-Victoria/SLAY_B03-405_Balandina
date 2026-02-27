#include "CSRMatrix.hpp"

#include <map>
#include <stdexcept>

CSRMatrix::CSRMatrix(size_t m, size_t n, const std::vector<DOK> &x) : m(m), n(n) {
    std::map<std::pair<size_t, size_t>, double> t;
    r.resize(m + 1);
    for (const DOK &d: x) {
        if (d.i >= m || d.j >= n) {
            throw std::invalid_argument("Matrix index out of bounds");
        }
        t[std::make_pair(d.i, d.j)] += d.a;
    }
    for (auto [p, a]: t) {
        auto [i, j] = p;
        if (a == 0) {
            continue;
        }
        v.push_back(a);
        c.push_back(j);
        r[i + 1]++;
    }
    for (size_t i = 1; i <= m; i++) {
        r[i] += r[i - 1];
    }
}

double CSRMatrix::at(size_t i, size_t j) {
    if (i >= m || j >= n) {
        throw std::runtime_error("Coordinates not in matrix");
    }
    size_t ql = r[i];
    size_t qr = r[i + 1];
    if (ql == qr) {
        return 0;
    }
    while (qr - ql > 1) {
        size_t qm = (qr + ql) / 2;
        if (c[qm] > j) {
            qr = qm;
        } else {
            ql = qm;
        }
    }
    if (c[ql] != j) {
        return 0;
    }
    return v[ql];
}

std::vector<double> CSRMatrix::operator*(const std::vector<double> &x) {
    if (x.size() != n) {
        throw std::runtime_error("Can't multiply to vector wrong size");
    }
    std::vector<double> y(m);
    for (size_t i = 0; i < m; i++) {
        for (size_t j = r[i]; j < r[i + 1]; j++) {
            y[i] += v[j] * x[c[j]];
        }
    }
    return y;
}
