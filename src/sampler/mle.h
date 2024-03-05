#pragma once

#include <vector>
#include <cmath>
#include "../pdf/generated_pdf.h"

double mlog_likelihood_mle(const double r,
                           const std::vector<std::vector<double>>& data
                           );
