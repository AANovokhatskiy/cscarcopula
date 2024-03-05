//Gumbel
//9
 
#pragma once
 
#include <vector>
#include <cmath>
 
std::vector<double> pdf(const std::vector<std::vector<double>>& u, const double r);
 
std::vector<double> pdf(const std::vector<std::vector<double>>& u, const std::vector<double>& r);
 
double pdf(const std::vector<double>& u, const double r);
 
double transform(const double r);
 
std::vector<double> transform(const std::vector<double>& r);
 
