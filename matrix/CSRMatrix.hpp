#pragma once
#include <vector>

struct DOK {
    size_t i;
    size_t j;
    double a;
};

class CSRMatrix {
public:
    CSRMatrix(size_t m, size_t n, const std::vector<DOK> &x);

    double at(size_t i, size_t j);

    std::vector<double> operator*(const std::vector<double> &x);

private:
    size_t m;
    size_t n;
    std::vector<double> v;
    std::vector<size_t> c;
    std::vector<size_t> r;
};
