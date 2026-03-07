#include "DenseMatrix.hpp"

double DenseMatrix::at(size_t i, size_t j) {
    if (i >= m || j >= n) {
        throw std::runtime_error("Coordinates not in matrix");
    }
    return data[i * n + j];
}

std::vector<double> DenseMatrix::operator*(const std::vector<double> &x) {
    if (x.size() != n) {
        throw std::runtime_error("Can't multiply to vector wrong size");
    }
    std::vector<double> y(m);
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < n; j++) {
            y[i] += data[i * n + j] * x[j];
        }
    }
    return y;
}
