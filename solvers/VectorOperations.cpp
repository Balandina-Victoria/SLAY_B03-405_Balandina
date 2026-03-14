#include "VectorOperations.hpp"

#include <stdexcept>

std::vector<double> operator*(const std::vector<double> &a, double b) {
    std::vector<double> res = a;
    for (double &x: res) {
        x *= b;
    }
    return res;
}

std::vector<double> operator*(double b, const std::vector<double> &a) {
    return a * b;
}

std::vector<double> &operator*=(std::vector<double> &a, double b) {
    for (double &x: a) {
        x *= b;
    }
    return a;
}

std::vector<double> operator+(const std::vector<double> &a, const std::vector<double> &b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vector sizes do not match");
    }
    std::vector<double> res = a;
    for (size_t i = 0; i < res.size(); i++) {
        res[i] += b[i];
    }
    return res;
}

std::vector<double> operator-(const std::vector<double> &a, const std::vector<double> &b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vector sizes do not match");
    }
    std::vector<double> res = a;
    for (size_t i = 0; i < res.size(); i++) {
        res[i] -= b[i];
    }
    return res;
}

std::vector<double> &operator+=(std::vector<double> &a, const std::vector<double> &b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vector sizes do not match");
    }
    for (size_t i = 0; i < a.size(); i++) {
        a[i] += b[i];
    }
    return a;
}

std::vector<double> &operator-=(std::vector<double> &a, const std::vector<double> &b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vector sizes do not match");
    }
    for (size_t i = 0; i < a.size(); i++) {
        a[i] -= b[i];
    }
    return a;
}

double operator*(const std::vector<double> &a, const std::vector<double> &b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vector sizes do not match");
    }
    double res = 0;
    for (size_t i = 0; i < a.size(); i++) {
        res += a[i] * b[i];
    }
    return res;
}

std::vector<double> operator%(const std::vector<double> &a, const std::vector<double> &b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vector sizes do not match");
    }
    std::vector<double> res = a;
    for (size_t i = 0; i < a.size(); i++) {
        res[i] *= b[i];
    }
    return res;
}
