#pragma  once

#include <chrono>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <cstddef>

#include "CSRMatrix.hpp"
#include "VectorOperations.hpp"

class Solver {
public:
    struct SolveResult {
        std::vector<double> x;
        int iters = 0;
        double time = 0;
        double res = 0;
    };

    Solver(const CSRMatrix &A, const std::vector<double> &b, const std::vector<double> &x0, double eps,
           int max_it) : A_(A), b_(b), x0_(x0), eps_(eps), max_it_(max_it) {
        if (A_.rows() != A_.cols()) {
            throw std::runtime_error("A must be square");
        }
        if (b_.size() != A_.rows()) {
            throw std::runtime_error("b size mismatch");
        }
        if (x0_.size() != A_.rows()) {
            throw std::runtime_error("x0 size mismatch");
        }
        if (eps_ <= 0) {
            throw std::runtime_error("eps must be > 0");
        }
        if (max_it_ <= 0) {
            throw std::runtime_error("max_it must be > 0");
        }
    }

    virtual ~Solver() = default;

    virtual SolveResult solve() = 0;

protected:
    const CSRMatrix &A_;
    const std::vector<double> &b_;
    std::vector<double> x0_;
    double eps_;
    int max_it_;

    static double l2_norm(const std::vector<double> &v) {
        return std::sqrt(v * v);
    }

    double residual(const std::vector<double> &x) const {
        return l2_norm(A_ * x - b_);
    }
};

class MPISolver : public Solver {
public:
    MPISolver(const CSRMatrix &A, const std::vector<double> &b, const std::vector<double> &x0, double eps, int max_it,
              double tau) : Solver(A, b, x0, eps, max_it), tau_(tau) {
        if (tau_ <= 0) {
            throw std::runtime_error("tau must be > 0");
        }
    }

    SolveResult solve() override {
        auto t0 = std::chrono::steady_clock::now();

        SolveResult res;
        res.x = x0_;

        for (res.iters = 0; res.iters < max_it_; res.iters++) {
            auto r = A_ * res.x - b_;
            res.x -= tau_ * r;
            res.res = residual(res.x);
            if (res.res < eps_) {
                break;
            }
        }

        auto t1 = std::chrono::steady_clock::now();
        res.time = std::chrono::duration<double, std::milli>(t1 - t0).count();
        res.res = residual(res.x);
        return res;
    }

private:
    double tau_;
};


class JacobiSolver : public Solver {
public:
    JacobiSolver(const CSRMatrix &A, const std::vector<double> &b, const std::vector<double> &x0, double eps,
                 int max_it) : Solver(A, b, x0, eps, max_it) {
        diag_inv_ = A.get_diag();
        for (auto &d: diag_inv_) {
            if (d == 0) {
                throw std::runtime_error("Diagonal element cant be zero for Jacobi");
            }
            d = 1 / d;
        }
    }

    SolveResult solve() override {
        auto t0 = std::chrono::steady_clock::now();

        SolveResult res;
        res.x = x0_;

        for (res.iters = 0; res.iters < max_it_; res.iters++) {
            res.x += (b_ - A_ * res.x) % diag_inv_;

            res.res = residual(res.x);
            if (res.res < eps_) {
                break;
            }
        }

        auto t1 = std::chrono::steady_clock::now();
        res.time = std::chrono::duration<double, std::milli>(t1 - t0).count();
        return res;
    }

private:
    std::vector<double> diag_inv_;
};

class GaussSeidelSolver : public Solver {
public:
    GaussSeidelSolver(const CSRMatrix &A, const std::vector<double> &b, const std::vector<double> &x0, double eps,
                      int max_it) : Solver(A, b, x0, eps, max_it) {
        diag_inv_ = A.get_diag();
        for (auto &d: diag_inv_) {
            if (d == 0) {
                throw std::runtime_error("Diagonal element cant be zero for GaussSeidel");
            }
            d = 1 / d;
        }
    }

    SolveResult solve() override {
        auto t0 = std::chrono::steady_clock::now();

        SolveResult res;
        res.x = x0_;

        const auto &r = A_.get_r();
        const auto &c = A_.get_c();
        const auto &v = A_.get_v();

        const size_t n = A_.rows();

        for (res.iters = 0; res.iters < max_it_; ++res.iters) {
            for (size_t i = 0; i < n; ++i) {
                double sum = 0.0;
                for (size_t k = r[i]; k < r[i + 1]; ++k) {
                    size_t j = c[k];
                    if (j != i) {
                        sum += v[k] * res.x[j];
                    }
                }
                res.x[i] = (b_[i] - sum) * diag_inv_[i];
            }

            res.res = residual(res.x);
            if (res.res < eps_) break;
        }

        auto t1 = std::chrono::steady_clock::now();
        res.time = std::chrono::duration<double, std::milli>(t1 - t0).count();
        return res;
    }

private:
    std::vector<double> diag_inv_;
};
