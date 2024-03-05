#include "auxiliary.h"

std::vector<std::vector<double>> log_returns(const std::vector<std::vector<double>>& arr_data)
{
    size_t n = arr_data.size();
    size_t dim = arr_data[0].size();
    std::vector<std::vector<double>>  log_return(n - 1, std::vector<double>(dim));
    for (size_t i = 0; i < n - 1; ++i)
    {
        for (size_t j = 0; j < dim; ++j)
        {
            log_return[i][j] = log(arr_data[i + 1][j] / arr_data[i][j]);
        }
    }
    return log_return;
}

std::vector<double> rank(const std::vector<double>& arr_data)
{
    size_t n = arr_data.size();
    std::vector<size_t> order(n);
    std::vector<double> ranks(n);

    // Create an index array
    for (size_t i = 0; i < n; ++i) {
        order[i] = i;
    }

    // Sort the index array based on the values in arr_data
    std::sort(order.begin(), order.end(),
              [&](size_t i, size_t j) { return arr_data[i] < arr_data[j]; });

    // Compute the ranks
    for (size_t i = 0; i < n; ++i) {
        ranks[order[i]] = static_cast<double>(i + 1) / (n + 1);
    }

    return ranks;
}

std::vector<std::vector<double>> pobs(const std::vector<std::vector<double>>& arr_data)
{
    size_t n = arr_data.size();
    size_t dim = arr_data[0].size();
    std::vector<std::vector<double>> res(n, std::vector<double>(dim));

    for (size_t k = 0; k < dim; ++k) {
        std::vector<double> column_data(n);
        for (size_t i = 0; i < n; ++i) {
            column_data[i] = arr_data[i][k];
        }
        std::vector<double> column_ranks = rank(column_data);
        for (size_t i = 0; i < n; ++i) {
            res[i][k] = column_ranks[i];
        }
    }

    return res;
}


// Linear least squares function using Eigen library
std::vector<double> linear_least_squares(const std::vector<std::vector<double>>& A,
                                         const std::vector<double>& b)
{

  // Convert input to Eigen types
  Eigen::MatrixXd eigenA(A.size(), A[0].size());
  Eigen::VectorXd eigenB(b.size());

  for (int i = 0; i < A.size(); ++i) {
    for (int j = 0; j < A[0].size(); ++j) {
      eigenA(i, j) = A[i][j];
    }
    eigenB(i) = b[i];
  }

  // Solve the linear least squares problem using Eigen's BDCSVD class
  Eigen::VectorXd x = eigenA.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(eigenB);

  // Convert the result back to std::vector
  std::vector<double> result(x.size());
  for (int i = 0; i < x.size(); ++i) {
    result[i] = x(i);
  }

  return result;
}


bool containsNaN(const std::vector<std::vector<double>>& vec)
{
    for (const auto& innerVec : vec) {
        for (const auto& value : innerVec) {
            if (std::isnan(value)) {
                return true;
            }
        }
    }
    return false;
}


template <typename T>
std::vector<T> multiply_vector(const std::vector<std::vector<T>>& A, const std::vector<T>& b)
{
    int dim11 = A.size();
    int dim12 = A[0].size();
    int dim2 = b.size();
    if (dim12 != dim2)
    {
         throw std::invalid_argument("Matricies shapes is incompatible. Got: (" + std::to_string(dim11) + ","
         + std::to_string(dim12) + ") and  (" + std::to_string(dim2) + ", ).");
    }
    std::vector<double> result(dim11, 0.0);
    for (int i = 0; i < dim11; ++i)
    {
        double temp = 0.0;
        for (int j = 0; j < dim12; ++j)
        {
            temp += A[i][j] * b[j];
        }
        result[i] = temp;
    }

    return result;
}

template std::vector<double> multiply_vector(const std::vector<std::vector<double>>& A, const std::vector<double>& b);
