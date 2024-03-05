#include "copula.h"

std::vector<std::vector<double>> calculate_crns(const int T, const int latent_process_tr)
{
    std::vector<std::vector<double>> crns(T, std::vector<double>(latent_process_tr));
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dis(0.0, 1.0);
    for (size_t i = 0; i < T; ++i)
    {
        for (size_t j = 0; j < latent_process_tr; ++j)
        {
            crns[i][j] = dis(gen);
        }
    }
    return crns;
}


double log_likelihood(const std::vector<std::vector<double>>& u, const double r)
{
    std::vector<double> pdf_data = pdf(u, r);
    size_t n = pdf_data.size();
    double log_lik = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        log_lik += log(pdf_data[i]);
    }
    return log_lik;
}


double log_likelihood(const std::vector<std::vector<double>>& u, const std::vector<double>& r)
{
    std::vector<double> pdf_data = pdf(u, r);
    size_t n = pdf_data.size();
    double log_lik = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        log_lik += log(pdf_data[i]);
    }
    return log_lik;
}


double fit_objective_mle(const std::vector<double>& x, std::vector<double>& grad, void* f_data)
{
    // for (const auto& val : x)
    // {
    //     std::cout << val << " ";
    // }
    fit_copula_args* args = static_cast<fit_copula_args*>(f_data);

    // std::vector<std::vector<double>>& pobs = args->pobs;
    const std::vector<std::vector<double>>& pobs = args->pobs;
    double mll = mlog_likelihood_mle(x[0], pobs);
    if (!grad.empty())
        {
            double h = 0.01;
            grad[0] = (mlog_likelihood_mle(x[0] + h, pobs) - mll) / h;
        }
    //std::cout<<mll<< std::endl;
    return mll;
}


double fit_objective_scar_p_ou(const std::vector<double>& x, std::vector<double>& grad, void* f_data)
{
    // for (const auto& val : x)
    // {
    //     std::cout << val << " ";
    // }
    fit_copula_args* args = static_cast<fit_copula_args*>(f_data);

    const std::vector<std::vector<double>>& pobs = args->pobs;
    const std::vector<std::vector<double>>& dwt = args->dwt;

    if (!grad.empty())
        {
            double h = 0.01;
            int dim = x.size();

            for (int k = 0; k < dim; ++k)
            {
                std::vector<double> x1(dim);
                x1 = x;
                x1[k] = x1[k] + h;
                grad[k] = (p_mlog_likelihood_ou(x1, pobs, dwt) - p_mlog_likelihood_ou(x, pobs, dwt)) / h;
            }
        }
    double mll = p_mlog_likelihood_ou(x, pobs, dwt);
    // std::cout<<mll<< std::endl;
    return mll;
}


double fit_objective_scar_m_ou(const std::vector<double>& x, std::vector<double>& grad, void* f_data)
{
    // for (const auto& val : x)
    // {
    //     std::cout << val << " ";
    // }
    fit_copula_args* args = static_cast<fit_copula_args*>(f_data);

    const std::vector<std::vector<double>>& pobs = args->pobs;
    const std::vector<std::vector<double>>& dwt = args->dwt;
    int m_iters = args->m_iters;

    if (!grad.empty())
        {
            double h = 0.01;
            int dim = x.size();

            for (int k = 0; k < dim; ++k)
            {
                std::vector<double> x1(dim);
                x1 = x;
                x1[k] = x1[k] + h;
                grad[k] = (m_mlog_likelihood_ou(x1, pobs, dwt, m_iters) - m_mlog_likelihood_ou(x, pobs, dwt, m_iters)) / h;
                //grad[k] = (p_mlog_likelihood_ou(x1, pobs, dwt) - p_mlog_likelihood_ou(x, pobs, dwt)) / h;
            }
        }
    double mll = m_mlog_likelihood_ou(x, pobs, dwt, m_iters);
    // std::cout<<mll<< std::endl;
    return mll;
}


std::vector<double> fit(const std::vector<std::vector<double>>& pobs,
                        const std::string& method,
                        const std::vector<std::vector<double>>& dwt,
                        const bool print_result)
{
    int T = pobs.size();
    double dt = 1.0/T;
    double sqrt_dt = sqrt(dt);
    int m_iters = 5;
    int latent_process_tr = dwt.size();

    fit_copula_args args{pobs, dwt, m_iters};

    double minf = 0.0;
    std::vector<double> x;
    nlopt::result result;

    if (method == "mle")
    {
        x = {0.5};
        double accuracy = 1e-5;
        nlopt::opt optimizer(nlopt::LD_SLSQP, x.size());

        optimizer.set_min_objective(fit_objective_mle, static_cast<void*>(&args));

        optimizer.set_ftol_rel(accuracy);

        result = optimizer.optimize(x, minf);
    }
    else if (method == "scar-p-ou")
    {
        double accuracy = 1e-3;
        x = {0.05, 0.95, 0.05};
        std::vector<double> lb = {-10.0,-3.999,0.001};
        std::vector<double> ub = {10.0, 3.999, 0.999};

        nlopt::opt optimizer(nlopt::LD_SLSQP, x.size());

        optimizer.set_min_objective(fit_objective_scar_p_ou, static_cast<void*>(&args));

        optimizer.set_lower_bounds(lb);
        optimizer.set_upper_bounds(ub);

        optimizer.set_ftol_rel(accuracy);

        result = optimizer.optimize(x, minf);
    }
    else if (method == "scar-m-ou")
    {
        double accuracy = 1e-3;
        x = {0.05, 0.95, 0.05};
        std::vector<double> lb = {-10.0,-3.999,0.001};
        std::vector<double> ub = {10.0, 3.999, 0.999};

        nlopt::opt optimizer(nlopt::LD_SLSQP, x.size());

        optimizer.set_min_objective(fit_objective_scar_m_ou, static_cast<void*>(&args));

        optimizer.set_lower_bounds(lb);
        optimizer.set_upper_bounds(ub);

        optimizer.set_ftol_rel(accuracy);

        result = optimizer.optimize(x, minf);
    }
    else
    {
        throw std::invalid_argument("Given method " + method + " is not implemented. Available methods: mle, scar-p-ou, scar-m-ou");
    }

    if (print_result)
    {
        if (result < 0)
        {
            std::cout << "Optimization failed. Error code: " << result << std::endl;
        }
        else
        {
            std::cout << "Optimization successful." << std::endl;
            std::cout << "Minimum value: " << minf << std::endl;
            std::cout << "Solution: ";
            for (const auto& val : x) {
                std::cout << val << " ";
            }
            std::cout << std::endl;
        }
    }
    // construct vector z = (minf, x)
    std::vector<double> z;
    z.reserve(x.size() + 1);
    z.push_back(minf);
    z.insert(z.end(), x.begin(), x.end());
    return z;

    //return x;
}

