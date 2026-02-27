#pragma once
#include <stdexcept>
#include <vector>

class DenseMatrix {
public:
    DenseMatrix(size_t m, size_t n, const std::vector<double> &x) : m(m), n(n) {
        if (x.size() != m * n) {
            throw std::runtime_error("Wrong size of data array");
        }
        data = x;
    }

    double at(size_t i, size_t j);

    std::vector<double> operator*(const std::vector<double> &x);

private:
    size_t m;
    size_t n;
    std::vector<double> data;
};
