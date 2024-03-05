#include "risk_metrics.h"

#include <omp.h>

double F_cvar_q(const std::vector<double>& q,
                const std::vector<double>& copula_pdf_data,
                const std::vector<double>& loss,
                const double gamma
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
                                                      const std::vector<std::vector<double>>& dwt,
                                                      const int window_len)
{
    std::cout << "calculation of latent process parameters" << std::endl;
    int T = pobs.size();
    if (T != dwt.size() and dwt.size() > 0)
    {
         throw std::invalid_argument("Size of data do not match to size of dwt. Data size = " + std::to_string(T) + ", dwt size = " + std::to_string(dwt.size()) + ".");
    }
    int numThreads = 8;
    omp_set_num_threads(numThreads);
    //auto num = omp_get_max_threads();
    int iters = T - window_len + 1;
    std::vector<std::vector<double>> res(T, std::vector<double>(4, 0.0));

    #pragma omp parallel
    {
    #pragma omp for
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
            // std::cout<<i<<" ";
            // for (auto& cell : res_i)
            // {
            //     std::cout<<cell<<" ";
            // }
            // std::cout<<std::endl;
            res[idx] = res_i;
        }
    }
    }
    return res;
}


std::vector<double> latent_process_sampler(const std::vector<double>& alpha,
                                           const std::string& copula_method,
                                           const int T,
                                           const int MC_iterations)
{
    std::vector<double> result;
    if (copula_method == "mle")
    {
        result = {alpha[0]};
    }
    else if (copula_method == "scar-p-ou")
    {
        result = p_sampler_no_hist_ou(alpha, T, MC_iterations);
    }
    else if (copula_method == "scar-m-ou")
    {
        result = p_sampler_no_hist_ou(alpha, T, MC_iterations);
    }
    else
    {
        throw std::invalid_argument("Given method " + copula_method + " is not implemented. Available methods: mle, scar-p-ou, scar-m-ou");
    }
    return result;
}


double fit_objective_cvar(const std::vector<double>& x, std::vector<double>& grad, void* f_data)
{
    fit_cvar_args* args = static_cast<fit_cvar_args*>(f_data);

    const std::vector<double>& copula_pdf_data = args->copula_pdf_data;
    const std::vector<double>& loss = args->loss;
    double gamma = args->gamma;

    double res = F_cvar_q(x, copula_pdf_data, loss, gamma);
    if (!grad.empty())
    {
        double h = 0.0001;
        std::vector<double> x1 = x;
        x1[0] = x[0] + h;
        grad[0] = (F_cvar_q(x1, copula_pdf_data, loss, gamma) - res) / h;
    }
    return res;
}


std::vector<double> fit_cvar(const std::vector<double>& weight,
                             const std::vector<double>& lp_params,
                             const std::vector<std::vector<double>>& marginals,
                             const std::string& copula_method,
                             const std::string& marginals_method,
                             const int T,
                             const int MC_iterations,
                             const double gamma)
{
    int dim = marginals.size();
    std::vector<std::vector<double>> returns_sample = rvs(marginals, MC_iterations, marginals_method);
    std::vector<std::vector<double>> pseudo_obs = pobs(returns_sample);
    std::vector<double> loss(MC_iterations, 0.0);
    for (int i = 0; i < MC_iterations; ++i)
    {
        double temp = 0.0;
        for (int j = 0; j < dim; ++j)
        {
            temp += exp(returns_sample[i][j]) * weight[j];
        }
        loss[i] = 1 - temp;
    }

    for (auto& innerVec : returns_sample)
    {
        innerVec.clear();
    }
    returns_sample.clear();

    std::vector<double> latent_process_sample = latent_process_sampler(lp_params, copula_method, T, MC_iterations);
    std::vector<double> copula_pdf_data;
    if (copula_method == "mle")
    {
        copula_pdf_data = pdf(pseudo_obs, transform(latent_process_sample[0]));
    }
    else
    {
        copula_pdf_data = pdf(pseudo_obs, transform(latent_process_sample));
    }

    for (auto& innerVec : pseudo_obs)
    {
        innerVec.clear();
    }
    pseudo_obs.clear();

    fit_cvar_args args{copula_pdf_data, loss, gamma};
    double minf = 0.0;
    std::vector<double> x;
    nlopt::result result;
    x = {0.0};
    double accuracy = 1e-5;
    nlopt::opt optimizer(nlopt::LD_SLSQP, x.size());
    optimizer.set_maxeval(10000);

    optimizer.set_min_objective(fit_objective_cvar, static_cast<void*>(&args));

    optimizer.set_ftol_rel(accuracy);

    result = optimizer.optimize(x, minf);

    std::vector<double> z(2, 0.0);
    z[0] = minf;
    z[1] = x[0];
    return z;
}

std::vector<std::vector<double>> cvar(const std::vector<std::vector<double>>& data,
                                      const std::vector<double>& weight,
                                      const std::string& copula_method,
                                      const std::string& marginals_method,
                                      const int latent_process_tr,
                                      const int MC_iterations,
                                      const int window_len,
                                      const double gamma
                                      )
{
    //std::vector<std::vector<double>> res = {{1.0, 2.0}, {3.0, 4.0}};
    int dim = data[0].size();
    if (dim != weight.size())
    {
         throw std::invalid_argument("Dimension of data does not match the dimension of weights. Data dim = "
         + std::to_string(dim) + ", Weights dim = " + std::to_string(weight.size()) + ".");
    }
    std::vector<std::vector<double>> pobs_data = pobs(data);
    int T = data.size();
    double dt = 1.0/window_len;
    double sqrt_dt = sqrt(dt);
    std::vector<std::vector<double>> dwt;
    if (copula_method != "mle")
    {
        dwt = calculate_crns(T, latent_process_tr);
        for(auto& row : dwt)
        {
            for(auto& el : row)
            {
                el *= sqrt_dt;
            }
        }
    }

    std::vector<std::vector<std::vector<double>>> marginals = marginals_params(data, window_len, marginals_method);
    std::vector<std::vector<double>> lp_params = latent_process_params(pobs_data, copula_method, dwt, window_len);

    for (auto& innerVec : dwt)
    {
        innerVec.clear();
    }
    dwt.clear();

    int iters = T - window_len + 1;
    std::vector<std::vector<double>> res(T, std::vector<double>(2, 0.0));

    std::cout << "calculation of risk metrics" << std::endl;

    int numThreads = 4;
    omp_set_num_threads(numThreads);
    #pragma omp parallel
    {
    #pragma omp for
    for (int i = 0; i < iters; ++i)
    {
        //std::cout<<i<<std::endl;
        int idx = (i + window_len - 1);
        std::vector<std::vector<double>> marginals_i = marginals[idx];
        std::vector<double> lp_params_i(lp_params[idx].begin() + 1, lp_params[idx].begin() + 4);
        std::vector<double> res_i;

        res_i = fit_cvar(weight,
                         lp_params_i,
                         marginals_i,
                         copula_method,
                         marginals_method,
                         window_len,
                         MC_iterations,
                         gamma);

        res[idx] = res_i;
    }
    }
    return res;
}



