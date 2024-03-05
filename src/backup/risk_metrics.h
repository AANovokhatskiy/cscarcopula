#pragma once

#include <vector>
#include <cmath>
#include "marginals.h"
#include "copula.h"
#include <omp.h>

double F_cvar_q(const std::vector<double>& q,
                const double gamma,
                std::vector<double>& loss,
                std::vector<double>& copula_pdf_data
                )
{
    int n = copula_pdf_data.size();
    double mean = 0.0;
    for (int k = 0; k < n; ++k)
    {
        mean += copula_pdf_data[k] * std::max(loss[k] - q[0], 0.0);
    }
    mean = mean / n;
    double F = q[0] + 1 / (1 - gamma) * mean;
    return F;
}

std::vector<std::vector<double>>latent_process_params(const std::vector<std::vector<double>>& pobs,
                                                      const std::string& method,
                                                      const std::vector<std::vector<double>>& dwt = {},
                                                      const int window_len = 0)
{
    int numThreads = 8;
    //omp_set_num_threads(numThreads);
    auto num = omp_get_max_threads();
    std::cout<<num<<std::endl;
    int T = pobs.size();
    int iters = T - window_len + 1;
    std::vector<std::vector<double>> res(T, std::vector<double>(4, 0.0));
    //#pragma omp for
    for (int i = 0; i < iters; ++i)
    {
        int idx = (i + window_len - 1);
        std::vector<std::vector<double>> sub_pobs(pobs.begin() + i, pobs.begin() + i + window_len);

        if (method == "mle")
        {
            std::vector<double> res_i = fit(sub_pobs, method);
            res[idx][0] = res_i[0];
            res[idx][1] = res_i[1];
        }
        else
        {
            std::vector<std::vector<double>> sub_dwt(dwt.begin() + i, dwt.begin() + i + window_len);
            std::vector<double> res_i = fit(sub_pobs, method, sub_dwt);
            std::cout<<i<<" ";
            for (auto& cell : res_i)
            {
                std::cout<<cell<<" ";
            }
            std::cout<<std::endl;
            res[idx] = res_i;
        }
    }
    return res;
}
