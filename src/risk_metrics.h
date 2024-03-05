#pragma once

#include <vector>
#include <cmath>
#include "marginals.h"
#include "copula.h"

struct fit_cvar_args
{
    const std::vector<double>& copula_pdf_data;
    const std::vector<double>& loss;
    const double gamma;
};

double F_cvar_q(const std::vector<double>& q,
                const double gamma,
                std::vector<double>& loss,
                std::vector<double>& copula_pdf_data
                );

std::vector<std::vector<double>>latent_process_params(const std::vector<std::vector<double>>& pobs,
                                                      const std::string& method,
                                                      const std::vector<std::vector<double>>& dwt = {},
                                                      const int window_len = 0);

std::vector<double> latent_process_sampler(const std::vector<double>& alpha,
                                           const std::vector<std::vector<double>>& dwt,
                                           const std::string& copula_method,
                                           const int MC_iterations);

double fit_objective_cvar(const std::vector<double>& x, std::vector<double>& grad, void* f_data);

std::vector<double> fit_cvar(const std::vector<double>& weight,
                             const std::vector<double>& lp_params,
                             const std::vector<std::vector<double>>& marginals,
                             const std::string& copula_method,
                             const std::string& marginals_method,
                             const int T,
                             const int MC_iterations,
                             const double gamma);

std::vector<std::vector<double>> cvar(const std::vector<std::vector<double>>& data,
                                      const std::vector<double>& weight,
                                      const std::string& copula_method,
                                      const std::string& marginals_method,
                                      const int latent_process_tr,
                                      const int MC_iterations,
                                      const int window_len,
                                      const double gamma
                                      );
