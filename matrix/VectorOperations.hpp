#pragma once
#include <vector>

std::vector<double> operator*(const std::vector<double> &a, double b);

std::vector<double> &operator*=(std::vector<double> &a, double b);

std::vector<double> operator+(const std::vector<double> &a, const std::vector<double> &b);

std::vector<double> &operator+=(std::vector<double> &a, const std::vector<double> &b);

double operator*(const std::vector<double> &a, const std::vector<double> &b);
