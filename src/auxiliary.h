#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#define EIGEN_DONT_PARALLELIZE
#include <Eigen/Dense>

std::vector<std::vector<double>> log_returns(const std::vector<std::vector<double>>& arr_data);

std::vector<double> rank(const std::vector<double>& arr_data);

std::vector<std::vector<double>> pobs(const std::vector<std::vector<double>>& arr_data);


std::vector<double> linear_least_squares(const std::vector<std::vector<double>>& A,
                                         const std::vector<double>& b);

bool containsNaN(const std::vector<std::vector<double>>& vec);

template <typename T>
std::vector<T> multiply_vector(const std::vector<std::vector<T>>& A, const std::vector<T>& b);
