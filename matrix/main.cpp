#include "CSRMatrix.hpp"
#include "DenseMatrix.hpp"
#include "VectorOperations.hpp"
#include <chrono>
#include <iostream>
#include <random>

int main() {
    size_t m = 5000;
    size_t n = 5000;

    std::vector<double> zeros = {0.5, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000004};

    std::mt19937 gen(37109231);
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    for (auto chance: zeros) {
        double sum_csr = 0;
        double sum_dense = 0;
        const int K = 2;
        for (int k = 0; k < K; k++) {
            std::vector<double> data(n * m);
            std::vector<DOK> dok;
            for (size_t i = 0; i < m; i++) {
                for (size_t j = 0; j < n; j++) {
                    if (dis(gen) <= chance) {
                        data[i * n + j] = gen();
                        dok.push_back({i, j, data[i * n + j]});
                    }
                }
            }

            std::vector<double> x(n);
            for (size_t i = 0; i < n; i++) {
                x[i] = gen();
            }
            CSRMatrix csr(m, n, dok);
            DenseMatrix dense(m, n, data);
            auto t0 = std::chrono::steady_clock::now();
            auto y1 = csr * x;
            auto t1 = std::chrono::steady_clock::now();
            auto y2 = dense * x;
            auto t2 = std::chrono::steady_clock::now();

            std::chrono::duration<double, std::milli> t_c = t1 - t0;
            std::chrono::duration<double, std::milli> t_d = t2 - t1;

            std::cout << "    Time for " << dok.size() << " not zero. CSR: " << t_c.count() << " ms, Dense: " << t_d.
                    count() <<
                    " ms"
                    << std::endl;
            sum_csr += t_c.count();
            sum_dense += t_d.count();
        }
        std::cout << "Average time for not zero chance " << chance << " is CSR: " << sum_csr / K << " ms, Dense: " <<
                sum_dense / K << " ms" << std::endl;
    }
    return 0;
}

/*
    Time for 12501073 not zero. CSR: 17.1585 ms, Dense: 25.2088 ms
    Time for 12496457 not zero. CSR: 16.643 ms, Dense: 25.5265 ms
Average time for not zero chance 0.5 is CSR: 16.9007 ms, Dense: 25.3676 ms
    Time for 2497674 not zero. CSR: 3.45444 ms, Dense: 24.6095 ms
    Time for 2497422 not zero. CSR: 3.51957 ms, Dense: 25.558 ms
Average time for not zero chance 0.1 is CSR: 3.487 ms, Dense: 25.0837 ms
    Time for 1250205 not zero. CSR: 1.84384 ms, Dense: 24.7531 ms
    Time for 1248339 not zero. CSR: 1.88868 ms, Dense: 25.8056 ms
Average time for not zero chance 0.05 is CSR: 1.86626 ms, Dense: 25.2793 ms
    Time for 249835 not zero. CSR: 0.404409 ms, Dense: 25.0091 ms
    Time for 249764 not zero. CSR: 0.370495 ms, Dense: 24.9509 ms
Average time for not zero chance 0.01 is CSR: 0.387452 ms, Dense: 24.98 ms
    Time for 24903 not zero. CSR: 0.067756 ms, Dense: 25.2844 ms
    Time for 25001 not zero. CSR: 0.0693 ms, Dense: 25.2552 ms
Average time for not zero chance 0.001 is CSR: 0.068528 ms, Dense: 25.2698 ms
    Time for 2567 not zero. CSR: 0.059801 ms, Dense: 28.0561 ms
    Time for 2522 not zero. CSR: 0.033152 ms, Dense: 25.9538 ms
Average time for not zero chance 0.0001 is CSR: 0.0464765 ms, Dense: 27.0049 ms
    Time for 244 not zero. CSR: 0.018814 ms, Dense: 27.3299 ms
    Time for 247 not zero. CSR: 0.019854 ms, Dense: 25.3164 ms
Average time for not zero chance 1e-05 is CSR: 0.019334 ms, Dense: 26.3232 ms
    Time for 23 not zero. CSR: 0.011144 ms, Dense: 25.4546 ms
    Time for 24 not zero. CSR: 0.012855 ms, Dense: 25.8015 ms
Average time for not zero chance 1e-06 is CSR: 0.0119995 ms, Dense: 25.6281 ms
    Time for 8 not zero. CSR: 0.009894 ms, Dense: 25.2052 ms
    Time for 3 not zero. CSR: 0.009927 ms, Dense: 25.5501 ms
Average time for not zero chance 4e-07 is CSR: 0.0099105 ms, Dense: 25.3776 ms
*/
