#include <iostream>
#include <vector>
#include <iomanip>
#include <string>
#include <tuple>

#include "src/io/data.h"
#include "src/auxiliary.h"
#include "src/copula.h"
#include "src/risk_metrics.h"

int main() {

    std::vector<std::vector<double>> csv_data = read_file("../data/test.csv");
    std::vector<std::vector<double>> returns = log_returns(csv_data);
    //std::vector<std::vector<double>> pobs_data = pobs(returns);

    int latent_process_tr = 1000;
    int window_len = 250;
    double gamma = 0.95;
    int dim = returns[0].size();
    int N_mc = static_cast<int>(pow(10.0, 7));

    std::vector<double> weight(dim, 1.0/dim);
    std::vector<std::vector<double>> cvar_data = cvar(returns, weight, "scar-m-ou", "normal", latent_process_tr, N_mc, window_len, gamma);

    to_csv(cvar_data, "../risk_data/cvar_data.csv");
    return 0;
}
