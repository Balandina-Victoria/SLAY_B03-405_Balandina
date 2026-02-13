#include <iostream>
#include <vector>

using namespace std;

vector<double>
progonka(vector<double> under_diag, vector<double> main_diag, vector<double> over_diag, vector<double> out, size_t n) {
    if (under_diag.size() != n || over_diag.size() != n || main_diag.size() != n || out.size() != n) {
        throw runtime_error("Матрица не трёхдиагональная");
    }
    for(size_t i = 0; i < n; i++){
        if (abs(main_diag[i]) <= abs(under_diag[i]) + abs(over_diag[i])) {
            throw runtime_error("Матрица не удовлетворяет условию строго диагонального преобладания");
        }
    }
    vector<double> alpha(n);
    vector<double> betta(n);
    alpha[0] = -over_diag[0] / main_diag[0];
    betta[0] = out[0] / main_diag[0];
    for (size_t i = 1; i < n - 1; i++) {
        alpha[i] = -over_diag[i] / (under_diag[i] * alpha[i - 1] + main_diag[i]);
        betta[i] = (out[i] - under_diag[i] * betta[i - 1]) / (under_diag[i] * alpha[i - 1] + main_diag[i]);
    }
    betta[n - 1] =
            (out[n - 1] - under_diag[n - 1] * betta[n - 2]) / (under_diag[n - 1] * alpha[n - 2] + main_diag[n - 1]);
    alpha[n - 1] = 0;

    vector<double> x(n);
    x[n - 1] = betta[n - 1];
    for (size_t i = n - 1; i > 0; i--) {
        x[i - 1] = alpha[i - 1]*x[i] + betta[i - 1];
    }
    return x;

}

int main(){
    vector<double> under_diag = {0, 1, 2, 3};
    vector<double> main_diag = {3, 5, 4, 6};
    vector<double> over_diag = {-2, 3, -1, 0};
    vector<double> out = {2, 0, -1, 5};
    size_t n = 4;
    vector<double> answers = {double(68)/117, double(-5)/39, double(7)/351, double(289)/351};
    vector<double> x = progonka(under_diag, main_diag, over_diag, out, n);
    double eps = 0.000000001;
    for (size_t i = 0; i < n; i++){
        if (abs(x[i] - answers[i]) > eps){
            cout << "Wrong answer";
            return -1;
        }
    }
    cout << "Correct";
}