#pragma once

#include <cmath>
#include <random>
#include <iostream>
#include <vector>
#include <stdexcept>

std::vector<std::vector<std::vector<double>>> normal_params(const std::vector<std::vector<double>>& data,
                                                            int window_len);

std::vector<std::vector<std::vector<double>>> marginals_params(const std::vector<std::vector<double>>& data,
                                                               int window_len,
                                                               const std::string& method);


std::vector<std::vector<double>> normal_rvs(const std::vector<std::vector<double>>& params,
                                            int N);

std::vector<std::vector<double>> rvs(const std::vector<std::vector<double>>& params,
                                     int N,
                                     const std::string& method);





