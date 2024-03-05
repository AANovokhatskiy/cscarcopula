#pragma once

#include <vector>
#include <cmath>
#include <algorithm>
#include <random>

#include "../pdf/generated_pdf.h"
#include "../auxiliary.h"

std::vector<std::vector<double>> p_sampler_ou(const std::vector<double>& alpha,
                                                     const std::vector<std::vector<double>>& dwt);


std::vector<double> p_sampler_no_hist_ou(const std::vector<double>& alpha,
                                         const int T,
                                         const int latent_process_tr
                                         );


double get_avg_p_log_likelihood(const std::vector<std::vector<double>>& data,
                                       const std::vector<std::vector<double>>& lambda_data
                                       );


double p_mlog_likelihood_ou(const std::vector<double>& alpha,
                                   const std::vector<std::vector<double>>& data,
                                   const std::vector<std::vector<double>>& dwt,
                                   bool print_path = false
                                   );


std::vector<double> polynom_fit(const std::vector<double>& x,
                                const std::vector<double>& y,
                                int dim,
                                bool fit_intercept = true);


std::vector<double> polynom_values(const std::vector<double>& x,
                                   const std::vector<double>& coef,
                                   bool intercept = true);


std::vector<double> polynom_values_correction(const std::vector<double>& t_data,
                                              const std::vector<double>& coef,
                                              const std::vector<double>& alpha,
                                              bool intercept);


std::vector<double> moving_average(const std::vector<double>& a,
                                   int n = 3);


std::vector<double> correction(const std::vector<double>& t_data,
                               const std::vector<double>& x_data,
                               const std::vector<double>& alpha);


std::vector<std::vector<double>> m_sampler_ou(const std::vector<double>& alpha,
                                      const std::vector<double>& a1t,
                                      const std::vector<double>& a2t,
                                      const std::vector<std::vector<double>>& dwt);


std::vector<double> log_norm_ou(const std::vector<double>& alpha,
                          const std::vector<double>& a1,
                          const std::vector<double>& a2,
                          double t,
                          const std::vector<double>& x0);


double log_norm_ou(const std::vector<double>& alpha,
                   double a1,
                   double a2,
                   double t,
                   double x0);


double m_mlog_likelihood_ou(const std::vector<double>& alpha,
                                const std::vector<std::vector<double>>& data,
                                const std::vector<std::vector<double>>& dwt,
                                int m_iters,
                                bool print_path = false);

