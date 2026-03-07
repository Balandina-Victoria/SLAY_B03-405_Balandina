#include "DenseMatrix.hpp"

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

const double EPS = 0.0000000001;

void apply_householder(vector<vector<double> > &A, vector<double> &b,
                       const vector<double> &v, double beta, int k) {
    int n = A.size();

    for (int j = k; j < n; ++j) {
        double s = 0.0;
        for (int i = k; i < n; ++i) {
            s += v[i - k] * A[i][j];
        }
        for (int i = k; i < n; ++i) {
            A[i][j] -= beta * v[i - k] * s;
        }
    }

    double s = 0.0;
    for (int i = k; i < n; ++i) {
        s += v[i - k] * b[i];
    }
    for (int i = k; i < n; ++i) {
        b[i] -= beta * v[i - k] * s;
    }
}

vector<double> solve_qr(vector<vector<double> > A, vector<double> b) {
    int n = A.size();

    for (int k = 0; k < n; ++k) {
        vector<double> x(n - k);
        for (int i = k; i < n; ++i) {
            x[i - k] = A[i][k];
        }

        double norm = 0.0;
        for (double value: x) {
            norm += value * value;
        }
        norm = sqrt(norm);

        if (norm < EPS) {
            continue;
        }

        vector<double> v = x;
        if (x[0] >= 0.0) {
            v[0] += norm;
        } else {
            v[0] -= norm;
        }

        double vv = 0.0;
        for (double value: v) {
            vv += value * value;
        }

        double beta = 2.0 / vv;

        apply_householder(A, b, v, beta, k);

        for (int i = k + 1; i < n; ++i) {
            A[i][k] = 0.0;
        }
    }

    vector<double> answer(n);

    for (int i = n - 1; i >= 0; --i) {
        double sum = b[i];
        for (int j = i + 1; j < n; ++j) {
            sum -= A[i][j] * answer[j];
        }

        if (abs(A[i][i]) < EPS) {
            throw runtime_error("Matrix is singular");
        }

        answer[i] = sum / A[i][i];
    }

    return answer;
}

int main() {
    int n;
    cin >> n;

    vector<vector<double> > A(n, vector<double>(n));
    vector<double> b(n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> A[i][j];
        }
    }

    for (int i = 0; i < n; ++i) {
        cin >> b[i];
    }

    vector<double> x = solve_qr(A, b);

    cout << fixed << setprecision(10);
    for (double value: x) {
        cout << value << ' ';
    }
    cout << '\n';

    return 0;
}
