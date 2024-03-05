#pragma once

#include <cmath>
#include <random>
#include <iostream>
#include <vector>
#include <stdexcept>

std::vector<std::vector<std::vector<double>>> normal_params(const std::vector<std::vector<double>>& data,
                                                            int window_len)
{
    int T = data.size();
    int dim = data[0].size();
    std::vector<std::vector<std::vector<double>>> res(T, std::vector<std::vector<double>>(dim, std::vector<double>(2)));

    for (int i = 0; i < T - window_len + 1; ++i)
    {
        for (int j = 0; j < dim; ++j)
        {
            double mean = 0.0;
            double std_dev = 0.0;

            for (int k = i; k < i + window_len; ++k)
            {
                mean += data[k][j];
            }
            mean /= window_len;

            for (int k = i; k < i + window_len; ++k)
            {
                std_dev += (data[k][j] - mean) * (data[k][j] - mean);
            }
            std_dev = std::sqrt(std_dev / window_len);

            res[i + window_len - 1][j][0] = mean;
            res[i + window_len - 1][j][1] = std_dev;
        }
    }

    return res;
}

std::vector<std::vector<std::vector<double>>> marginals_params(const std::vector<std::vector<double>>& data,
                                                               int window_len,
                                                               const std::string& method)
{
    std::cout << "calc marginals_params" << std::endl;
    std::vector<std::vector<std::vector<double>>> res;

    if (method == "normal")
    {
        res = normal_params(data, window_len);
    }
    else
    {
        throw std::invalid_argument("given method " + method + " is not implemented. Available methods: normal");
    }

    return res;
}


std::vector<std::vector<double>> normal_rvs(const std::vector<std::vector<double>>& params,
                                            int N)
{
    int dim = params.size();
    std::vector<std::vector<double>> res(N, std::vector<double>(dim));

    std::random_device rd;
    std::mt19937 gen(rd());

    for (int i = 0; i < dim; ++i)
    {
        std::normal_distribution<> dist(params[i][0], params[i][1]);
        for (int j = 0; j < N; ++j)
        {
            res[j][i] = dist(gen);
        }
    }

    return res;
}

std::vector<std::vector<double>> rvs(const std::vector<std::vector<double>>& params,
                                     int N,
                                     const std::string& method)
{
    std::vector<std::vector<double>> res;
    if (method == "normal")
    {
        res = normal_rvs(params, N);
    }
    else
    {
        throw std::invalid_argument("given method " + method + " is not implemented. Available methods: normal");
    }
    return res;
}





