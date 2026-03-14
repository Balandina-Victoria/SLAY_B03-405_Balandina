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

    std::size_t rows() const {
        return m;
    }

    std::size_t cols() const {
        return n;
    }

    const std::vector<double> &get_v() const {
        return v;
    }

    const std::vector<size_t> &get_c() const {
        return c;
    }

    const std::vector<size_t> &get_r() const {
        return r;
    }

    std::vector<double> get_diag() const;

    double at(size_t i, size_t j);

    std::vector<double> operator*(const std::vector<double> &x) const;

private:
    size_t m;
    size_t n;
    std::vector<double> v;
    std::vector<size_t> c;
    std::vector<size_t> r;
};
