#pragma once

#include <stdexcept>
#include <cmath>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include "pdf/generated_pdf.h"
#include "sampler/scar-ou.h"
#include "sampler/mle.h"
//#include <eigen3/Eigen/Dense>
#include <Eigen/Dense>
#include <nlopt.hpp>

struct fit_copula_args
{
    const std::vector<std::vector<double>>& pobs;
    const std::vector<std::vector<double>>& dwt;
    const int m_iters;
};

std::vector<std::vector<double>> calculate_crns(const int T, const int latent_process_tr);

double log_likelihood(const std::vector<std::vector<double>>& u, const double r);

double log_likelihood(const std::vector<std::vector<double>>& u, const std::vector<double>& r);

double fit_objective_mle(const std::vector<double>& x, std::vector<double>& grad, void* f_data);

double fit_objective_scar_p_ou(const std::vector<double>& x, std::vector<double>& grad, void* f_data);

double fit_objective_scar_m_ou(const std::vector<double>& x, std::vector<double>& grad, void* f_data);

std::vector<double> fit(const std::vector<std::vector<double>>& pobs,
                        const std::string& method,
                        const std::vector<std::vector<double>>& dwt = {},
                        const bool print_result = false);

