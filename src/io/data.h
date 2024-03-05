#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

template <typename T>
void print_arr(const std::vector<std::vector<T>>& data, const int n = 0);

std::vector<std::vector<double>> read_file(const std::string& filename, const char sep = ',');

template <typename T>
void to_csv(const std::vector<std::vector<T>>& data, const std::string& path, const char sep = ',');
