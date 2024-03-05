#include "scar-ou.h"

std::vector<std::vector<double>> p_sampler_ou(const std::vector<double>& alpha,
                                              const std::vector<std::vector<double>>& dwt)
{
        double alpha1 = alpha[0];
        double alpha2 = alpha[1];
        double alpha3 = alpha[2];
        size_t T = dwt.size();
        size_t n = dwt[0].size();
        double dt = 1.0/T;
        std::vector<std::vector<double>> xt(T, std::vector<double>(n, 0.0));
        double mu = -alpha1 / alpha2;

        for (size_t i = 0; i < n; ++i)
        {
            xt[0][i] = mu;
        }
        for (size_t k = 1; k < T; ++k)
        {
            for (size_t j = 0; j < n; ++j)
            {
                xt[k][j] = xt[k - 1][j] + (alpha1 + alpha2 * xt[k - 1][j]) * dt + alpha3 * dwt[k][j];

            }
        }
        return xt;
}

std::vector<double> p_sampler_no_hist_ou(const std::vector<double>& alpha,
                                         const int T,
                                         const int latent_process_tr)
{
        double alpha1 = alpha[0];
        double alpha2 = alpha[1];
        double alpha3 = alpha[2];
        double dt = 1.0/T;
        double sqrt_dt = sqrt(dt);
        double mu = -alpha1 / alpha2;

        std::vector<double> xt(latent_process_tr, 0.0);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> dis(0.0, 1.0);

        for (size_t j = 0; j < latent_process_tr; ++j)
        {
            double x_km1 = mu;
            double x_k = 0.0;
            for (size_t k = 1; k < T; ++k)
            {
                x_k = x_km1 + (alpha1 + alpha2 * x_km1) * dt + alpha3 * sqrt_dt * dis(gen);
                x_km1 = x_k;
            }
            xt[j] = x_k;
        }
        return xt;
}



double get_avg_p_log_likelihood(const std::vector<std::vector<double>>& data,
                                       const std::vector<std::vector<double>>& lambda_data
                                       )
{
    double avg_likelihood = 0.0;
    size_t T = data.size();
    size_t latent_process_tr = lambda_data[0].size();
    std::vector<double> copula_log_data(latent_process_tr);
    for (size_t k = 0; k < latent_process_tr; ++k)
    {
        double res = 0.0;
        for (size_t j = 0; j < T; ++j)
        {
            res += log(pdf(data[j], transform(lambda_data[j][k])));
        }
        copula_log_data[k] = res;
    }
    double xc = *std::max_element(copula_log_data.begin(), copula_log_data.end());

    for (size_t k = 0; k < latent_process_tr; ++k)
    {
        avg_likelihood += exp(copula_log_data[k] - xc);
    }

    avg_likelihood = avg_likelihood / latent_process_tr;
    avg_likelihood = log(avg_likelihood) + xc;
    return avg_likelihood;
}

double p_mlog_likelihood_ou(const std::vector<double>& alpha,
                                   const std::vector<std::vector<double>>& data,
                                   const std::vector<std::vector<double>>& dwt,
                                   bool print_path
                                   )
{
    std::vector<std::vector<double>> lambda_data = p_sampler_ou(alpha, dwt);
    double avg_log_likelihood = get_avg_p_log_likelihood(data, lambda_data);
    double res = -avg_log_likelihood;
    if (print_path)
        {
            std::cout<<alpha[0]<<" "<<alpha[1]<<" "<<alpha[2]<< " "<< res<< std::endl;
        }
    return res;
}


std::vector<double> polynom_fit(const std::vector<double>& x,
                                const std::vector<double>& y,
                                int dim,
                                bool fit_intercept)
{
    int fi = fit_intercept ? 1 : 0;
    std::vector<std::vector<double>> A(x.size(), std::vector<double>(dim + fi, 0.0));

    for (size_t i = 0; i < x.size(); ++i)
    {
        for (int j = 0; j < dim; ++j)
        {
            A[i][j + fi] = std::pow(x[i], j + 1);
        }
        if (fit_intercept)
        {
            A[i][0] = 1.0;
        }
    }

    std::vector<double> res = linear_least_squares(A, y);

    return res;
}

std::vector<double> polynom_values(const std::vector<double>& x,
                                   const std::vector<double>& coef,
                                   bool intercept)
{
    int dim = coef.size();
    std::vector<double> res(x.size(), 0.0);
    int fi = intercept ? 1 : 0;

    for (size_t i = 0; i < x.size(); ++i)
    {
        double temp = 0.0;
        for (int j = 0; j < dim; ++j)
        {
            temp += coef[j] * std::pow(x[i], 1 - fi + j);
        }
        res[i] = temp;
    }

    return res;
}


std::vector<double> polynom_values_correction(const std::vector<double>& t_data,
                                              const std::vector<double>& coef,
                                              const std::vector<double>& alpha,
                                              bool intercept)
{
    std::vector<double> result = polynom_values(t_data, coef, intercept);
    double alpha1 = alpha[0];
    double alpha2 = alpha[1];
    double alpha3 = alpha[2];

    size_t dim = t_data.size();

    for (size_t i = 0; i < dim; ++i)
    {
        double sigma2 = alpha3 * alpha3 / (-2 * alpha2) * (1 - std::exp(2 * alpha2 * t_data[i])) + 0.0001;
        double max_res = 1 / (2 * sigma2) - 1;
        double exp_res = std::exp(-0.05 * (max_res - result[i]));
        result[i] = 1 / (1 + exp_res) * result[i];
    }

    return result;
}

std::vector<double> moving_average(const std::vector<double>& a,
                                   int n)
{
    std::vector<double> ret(a.size(), 0.0);
    std::vector<double> cumsum(a.size(), 0.0);

    cumsum[0] = a[0];
    for (size_t i = 1; i < a.size(); ++i)
    {
        cumsum[i] = cumsum[i - 1] + a[i];
    }

    ret[n - 1] = cumsum[n - 1] / n;
    for (size_t i = n; i < a.size(); ++i)
    {
        ret[i] = (cumsum[i] - cumsum[i - n]) / n;
    }

    std::vector<double> lin_arr1, lin_arr2;
    for (int i = 0; i < n / 2; ++i)
    {
        lin_arr1.push_back(a[0] + (ret[0] - a[0]) * (i + 1) / (n / 2));
        lin_arr2.push_back(ret.back() + (a.back() - ret.back()) * (i + 1) / (n / 2));
    }

    ret.insert(ret.begin(), lin_arr1.begin(), lin_arr1.end());
    ret.insert(ret.end(), lin_arr2.begin(), lin_arr2.end());

    return ret;
}

std::vector<double> correction(const std::vector<double>& t_data,
                               const std::vector<double>& x_data,
                               const std::vector<double>& alpha)
{
    double alpha1 = alpha[0];
    double alpha2 = alpha[1];
    double alpha3 = alpha[2];

    std::vector<double> sigma2(t_data.size(), 0.0);
    for (size_t i = 0; i < t_data.size(); ++i) {
        sigma2[i] = alpha3 * alpha3 / (-2 * alpha2) * (1 - std::exp(2 * alpha2 * t_data[i])) + 0.0001;
    }

    std::vector<double> max_res(t_data.size(), 0.0);
    for (size_t i = 0; i < t_data.size(); ++i) {
        max_res[i] = 1 / (2 * sigma2[i]) - 1;
    }

    std::vector<double> exp_res(t_data.size(), 0.0);
    for (size_t i = 0; i < t_data.size(); ++i) {
        exp_res[i] = std::exp(-0.05 * (max_res[i] - x_data[i]));
    }

    std::vector<double> result(t_data.size(), 0.0);
    for (size_t i = 0; i < t_data.size(); ++i) {
        result[i] = 1 / (1 + exp_res[i]) * x_data[i];
    }

    return result;
}

std::vector<std::vector<double>> m_sampler_ou(const std::vector<double>& alpha,
                                      const std::vector<double>& a1t,
                                      const std::vector<double>& a2t,
                                      const std::vector<std::vector<double>>& dwt)
{
    double alpha1 = alpha[0];
    double alpha2 = alpha[1];
    double alpha3 = alpha[2];
    int T = dwt.size();
    int dim = dwt[0].size();
    double dt = 1.0 / T;
    std::vector<std::vector<double>> xt(T, std::vector<double>(dim, 0.0));
    double mu = -alpha1 / alpha2;
    double theta = -alpha2;
    double nu = alpha3;
    double x0 = mu;
    xt[0] = std::vector<double>(dim, x0);

    for (int i = 1; i < T; ++i) {
        double a1 = a1t[i];
        double a2 = a2t[i];
        double a1dt = (a1t[i] - a1t[i - 1]) / dt;
        double a2dt = (a2t[i] - a2t[i - 1]) / dt;
        double t = static_cast<double>(i) / T;
        double xs = (x0 - mu) * std::exp(alpha2 * t);
        double xsdt = -theta * xs;
        double sigma2 = alpha3 * alpha3 / (-2 * alpha2) * (1 - std::exp(2 * alpha2 * t));
        double sigma2dt = nu * nu - 2 * theta * sigma2;
        double p = (1 - 2 * a2 * sigma2);
        double sigma2w = sigma2 / p;
        double sigma2wdt = (sigma2dt + 2 * sigma2 * sigma2 * a2dt) / (p * p);
        double xsw = (xs + a1 * sigma2) / p;
        double xswdt = (xsdt + a1 * sigma2dt + a1dt * sigma2) / p + 2 * xsw * (a2dt * sigma2 + a2 * sigma2dt) / p;

        for (int j = 0; j < dim; ++j) {
            double B = nu;
            double A = xswdt - (xt[i - 1][j] - mu - xsw) * (B * B - sigma2wdt) / (2 * sigma2w);
            xt[i][j] = xt[i - 1][j] + A * dt + B * dwt[i][j];
        }
    }

    return xt;
}

std::vector<double> log_norm_ou(const std::vector<double>& alpha,
                          const std::vector<double>& a1,
                          const std::vector<double>& a2,
                          double t,
                          const std::vector<double>& x0)
{
    double alpha1 = alpha[0];
    double alpha2 = alpha[1];
    double alpha3 = alpha[2];
    double mu = -alpha1 / alpha2;
    double sigma2 = alpha3 * alpha3 / (-2 * alpha2) * (1 - std::exp(2 * alpha2 * t));

    std::vector<double> xs(x0.size());
    std::vector<double> res(x0.size());

    for (size_t i = 0; i < x0.size(); ++i)
    {
        xs[i] = (x0[i] - mu) * std::exp(alpha2 * t);
        res[i] = (2 * xs[i] * (a1[i] + a2[i] * xs[i]) + a1[i] * a1[i] * sigma2) / (2 - 4 * a2[i] * sigma2) - 0.5 * std::log(1 - 2 * a2[i] * sigma2);
    }

    return res;
}


double log_norm_ou(const std::vector<double>& alpha,
                   double a1,
                   double a2,
                   double t,
                   double x0)
{
    double alpha1 = alpha[0];
    double alpha2 = alpha[1];
    double alpha3 = alpha[2];
    double mu = -alpha1 / alpha2;
    double xs = (x0 - mu) * std::exp(alpha2 * t);
    double sigma2 = alpha3 * alpha3 / (-2 * alpha2) * (1 - std::exp(2 * alpha2 * t));
    double res = (2 * xs * (a1 + a2 * xs) + a1 * a1 * sigma2) / (2 - 4 * a2 * sigma2) - 0.5 * std::log(1 - 2 * a2 * sigma2);
    return res;
}


double m_mlog_likelihood_ou(const std::vector<double>& alpha,
                                const std::vector<std::vector<double>>& data,
                                const std::vector<std::vector<double>>& dwt,
                                int m_iters,
                                bool print_path)
{
    int T = data.size();
    int latent_process_tr = dwt[0].size();
    std::vector<std::vector<double>> norm_log_data(T, std::vector<double>(latent_process_tr, 0.0));
    std::vector<std::vector<double>> lambda_data(T, std::vector<double>(latent_process_tr, 0.0));
    double dt = 1.0 / T;
    std::vector<double> t_data(T);
    std::vector<std::vector<double>> a_data(T, std::vector<double>(3, 0.0));
    double alpha1 = alpha[0];
    double alpha2 = alpha[1];
    double alpha3 = alpha[2];
    double mu = -alpha1 / alpha2;
    std::vector<double> a1t(T, 0.0);
    std::vector<double> a2t(T, 0.0);

    for (int i = 0; i < T; ++i) {
        t_data[i] = static_cast<double>(i) / T;
    }

    for (int j = 0; j < m_iters; ++j)
    {
        if (j == 0)
        {
            lambda_data = p_sampler_ou(alpha, dwt);
        }
        else
        {
            lambda_data = m_sampler_ou(alpha, a1t, a2t, dwt);
            if (containsNaN(lambda_data))
            {
                double res = 10000;
                if (print_path)
                {
                    std::cout << alpha1<< " "<< alpha2 <<" "<< alpha3 << " m sampler nan " << res << std::endl;
                }
                return res;
            }
        }
        for (int i = 1; i < T; ++i)
        {
            std::vector<std::vector<double>> A(latent_process_tr, std::vector<double>(3, 0.0));
            std::vector<double> b(latent_process_tr);
            double sigma2 = alpha3 * alpha3 / (-2 * alpha2) * (1 - std::exp(2 * alpha2 * t_data[i]));

            for (int k = 0; k < latent_process_tr; ++k)
            {
                double copula_log_data = std::log(pdf(data[i], transform(lambda_data[i][k])));
                A[k][0] = 1.0;
                A[k][1] = lambda_data[i][k] - mu;
                A[k][2] = (lambda_data[i][k] - mu) * (lambda_data[i][k] - mu);
                b[k] = copula_log_data + norm_log_data[i - 1][k];
            }

            try {
                a_data[i] = linear_least_squares(A, b);
                a_data[i][2] = std::min(a_data[i][2], 1 / (2 * sigma2) - 10);
                a_data[i] = {std::max(std::min(a_data[i][0], 30.0), -30.0),
                             std::max(std::min(a_data[i][1], 30.0), -30.0),
                             std::max(std::min(a_data[i][2], 30.0), -30.0)};
            }
            catch (const std::exception&)
            {
                double res = 10000;
                if (print_path)
                {
                    std::cout << alpha1<< " "<< alpha2 <<" "<< alpha3 << " ls problem fail " << res << " " << i << std::endl;
                }
                return res;
            }
            for (int k = 0; k < latent_process_tr; ++k)
            {
                norm_log_data[i][k] = log_norm_ou(alpha, a_data[i][1], a_data[i][2], dt, lambda_data[i - 1][k]);

            }
        }
        std::vector<double> a_data_a1(T);
        std::vector<double> a_data_a2(T);
        for (int i = 0; i < T; ++i)
        {
            a_data_a1[i] = a_data[i][1];
            a_data_a2[i] = a_data[i][2];
        }
        int dim = 2;
        bool intercept = false;
        std::vector<double> a1_params = polynom_fit(t_data, a_data_a1, dim, intercept);
        std::vector<double> a2_params = polynom_fit(t_data, a_data_a2, dim, intercept);
        a1t = polynom_values(t_data, a1_params, intercept);
        a2t = polynom_values_correction(t_data, a2_params, alpha, intercept);

        //int n = T / 10;
        //a1t = moving_average(a_data_a1, n);
        //a2t = correction(t_data, moving_average(a_data_a2, n), alpha);
    }

    std::vector<double> log_likelihood(latent_process_tr, 0.0);
    lambda_data = m_sampler_ou(alpha, a1t, a2t, dwt);
    //lambda_data = p_sampler_ou(alpha, dwt);
    for (int i = 1; i < T; ++i)
    {
        double a1 = a1t[i];
        double a2 = a2t[i];

        for (int k = 0; k < latent_process_tr; ++k)
        {
            norm_log_data[i][k] = log_norm_ou(alpha, a1, a2, dt, lambda_data[i - 1][k]);
        }
    }
    for (int k = 0; k < latent_process_tr; ++k)
    {
        double temp = 0.0;
        for (int i = 0; i < T; ++i)
        {
            double copula_log_data = log(pdf(data[i], transform(lambda_data[i][k])));
            double g = a1t[i] * (lambda_data[i][k] - mu) + a2t[i] * (lambda_data[i][k] - mu) * (lambda_data[i][k] - mu);
            temp += copula_log_data + norm_log_data[i][k] - g;
        }
        log_likelihood[k] = temp;
    }
    double xc = *std::max_element(log_likelihood.begin(), log_likelihood.end());

    double avg_likelihood = 0.0;
    for (int k = 0; k < latent_process_tr; ++k)
    {
        avg_likelihood += exp(log_likelihood[k] - xc);
    }
    avg_likelihood = avg_likelihood / latent_process_tr;

    double res = std::log(avg_likelihood) + xc;
    res = -res;
    if (print_path) {
        std::cout << alpha1<< " "<< alpha2 <<" "<< alpha3  << " " << res << std::endl;
    }
    return res;

}




