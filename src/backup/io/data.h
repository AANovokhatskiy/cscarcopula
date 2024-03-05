#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

template <typename T>
void print_arr(const std::vector<std::vector<T>>& data, int n = 0)
{
    int rows_to_print = (n == 0) ? data.size() : std::min(n, static_cast<int>(data.size()));

    for (int i = 0; i < rows_to_print; ++i)
    {
        for (const auto& cell : data[i])
        {
            std::cout << cell << " ";
        }
        std::cout << std::endl;
    }
}


std::vector<std::vector<double>> read_file(const std::string& filename)
{
    std::ifstream file(filename);
    std::vector<std::vector<double>> data;

    if (file) {
        std::string line;
        bool first_line = true;
        while (std::getline(file, line)) {
            if (first_line) {
                first_line = false;
                continue;  // Skip the first line (header)
            }
            std::vector<double> row;
            std::stringstream ss(line);
            std::string cell;
            bool first_cell = true;
            while (std::getline(ss, cell, ',')) {
                if (first_cell) {
                    first_cell = false;
                    continue;  // Skip the first cell (date)
                }
                row.push_back(std::stod(cell));
            }
            data.push_back(row);
        }
    } else {
        std::cerr << "Error: Unable to open the file." << std::endl;
    }

    return data;
}

