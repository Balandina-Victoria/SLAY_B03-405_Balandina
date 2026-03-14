#include "CSRMatrix.hpp"
#include "VectorOperations.hpp"
#include "solvers.hpp"
#include <iostream>
#include <random>

std::mt19937 gen(1121141);

double gn() {
    return gen() % 10 + 1;
}

CSRMatrix get_matrix(size_t n, double dense) {
    std::vector<DOK> doks;
    size_t items = n * n * dense;
    items = std::max(n, items) - n;
    std::vector<double> sm(n);
    while (items--) {
        double x = gn();
        size_t i = gen() % n;
        size_t j = gen() % n;
        sm[i] += x;
        doks.push_back({i, j, x});
    }
    for (size_t i = 0; i < n; i++) {
        doks.push_back({i, i, sm[i] + gn() * 2});
    }
    return CSRMatrix(n, n, doks);
}

void run(size_t n, double dense, double eps, int max_it, double tau) {
    CSRMatrix A = get_matrix(n, dense);
    std::vector<double> x(n);
    for (size_t i = 0; i < n; i++) {
        x[i] = gn();
    }
    std::vector<double> b = A * x;
    std::vector<double> x0(n);
    MPISolver mpi(A, b, x0, eps, max_it, tau);
    JacobiSolver jacobi(A, b, x0, eps, max_it);
    GaussSeidelSolver gauss_seidel(A, b, x0, eps, max_it);
    auto r_mpi = mpi.solve();
    auto r_jac = jacobi.solve();
    auto r_gs = gauss_seidel.solve();


    std::cout << "n: " << n << ", dense: " << dense << ", eps: " << eps << ", max_it: " << max_it << ", tau:" << tau <<
            "\n";
    std::cout << "MPI: " << r_mpi.iters << " iters," << r_mpi.time << " ms," << r_mpi.res << " res\n";
    std::cout << "Jacobi: " << r_jac.iters << " iters," << r_jac.time << " ms," << r_jac.res << " res\n";
    std::cout << "GaussSeidel: " << r_gs.iters << " iters," << r_gs.time << " ms," << r_gs.res << " res\n";
    std::cout << std::endl << std::endl;
}

int main() {
    run(1000, 0.005, 1e-8, 1000, 0.01);
    run(10000, 0.0005, 1e-8, 1000, 0.01);
    run(10000, 0.005, 1e-8, 1000, 0.01);
    run(50000, 0.0001, 1e-8, 1000, 0.01);
    run(50000, 0.001, 1e-8, 1000, 0.01);
    run(70000, 0.0001, 1e-8, 1000, 0.01);
    run(70000, 0.0006, 1e-8, 1000, 0.01);
    return 0;
}
/*
n: 1000, dense: 0.005, eps: 1e-08, max_it: 1000, tau:0.01
MPI: 1000 iters,23.6571 ms,1.76649e-08 res
Jacobi: 62 iters,1.29459 ms,6.97819e-09 res
GaussSeidel: 17 iters,0.520815 ms,4.84306e-09 res


n: 10000, dense: 0.0005, eps: 1e-08, max_it: 1000, tau:0.01
MPI: 1000 iters,264.198 ms,1.26226e-07 res
Jacobi: 63 iters,16.5827 ms,7.98086e-09 res
GaussSeidel: 18 iters,5.99797 ms,7.6172e-09 res


n: 10000, dense: 0.005, eps: 1e-08, max_it: 1000, tau:0.01
MPI: 1000 iters,1220.46 ms,-nan res
Jacobi: 761 iters,893.93 ms,9.89216e-09 res
GaussSeidel: 19 iters,24.4574 ms,6.72282e-09 res


n: 50000, dense: 0.0001, eps: 1e-08, max_it: 1000, tau:0.01
MPI: 1000 iters,1483.95 ms,1.91984e-07 res
Jacobi: 65 iters,97.9184 ms,9.33254e-09 res
GaussSeidel: 19 iters,31.3698 ms,2.42994e-09 res


n: 50000, dense: 0.001, eps: 1e-08, max_it: 1000, tau:0.01
MPI: 1000 iters,16333.7 ms,-nan res
Jacobi: 781 iters,12784.8 ms,9.63295e-09 res
GaussSeidel: 20 iters,351.558 ms,4.16225e-09 res


n: 70000, dense: 0.0001, eps: 1e-08, max_it: 1000, tau:0.01
MPI: 1000 iters,2698.68 ms,7.95756e-08 res
Jacobi: 96 iters,263.673 ms,7.70654e-09 res
GaussSeidel: 18 iters,52.0686 ms,4.837e-09 res


n: 70000, dense: 0.0006, eps: 1e-08, max_it: 1000, tau:0.01
MPI: 1000 iters,24490.1 ms,-nan res
Jacobi: 659 iters,15715.5 ms,9.99926e-09 res
GaussSeidel: 20 iters,496.459 ms,3.30771e-09 res
 */
