#include "mle.h"

double mlog_likelihood_mle(const double r,
                           const std::vector<std::vector<double>>& data
                           )
{
    double log_likelihood = 0.0;
    int T = data.size();
    double transformed_r = transform(r);
    for (int i = 0; i < T; ++i)
    {
        log_likelihood += log(pdf(data[i], transformed_r));
    }
    return -log_likelihood;
}
