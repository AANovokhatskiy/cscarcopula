#include "data.h"

template <typename T>
void print_arr(const std::vector<std::vector<T>>& data, const int n)
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

template void print_arr<double>(const std::vector<std::vector<double>>& data, const int n);


std::vector<std::vector<double>> read_file(const std::string& filename, const char sep)
{
    std::ifstream file(filename);
    std::vector<std::vector<double>> data;

    if (file)
    {
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
            while (std::getline(ss, cell, sep))
            {
                if (first_cell) {
                    first_cell = false;
                    continue;  // Skip the first cell (date)
                }
                row.push_back(std::stod(cell));
            }
            data.push_back(row);
        }
    }
    else
    {
        std::cerr << "Error: Unable to open the file." << std::endl;
    }

    return data;
}

template <typename T>
void to_csv(const std::vector<std::vector<T>>& data, const std::string& path, const char sep)
{
    std::ofstream file(path);

    if (!file.is_open())
    {
        throw std::runtime_error("Could not open file");
    }

    for (const auto& row : data) {
        for (size_t i = 0; i < row.size(); ++i)
        {
            file << row[i];

            if (i != row.size() - 1) {
                file << sep;
            }
        }

        file << "\n";
    }

    file.close();
}

template void to_csv(const std::vector<std::vector<double>>& data, const std::string& path, const char sep);

