This project is made for calculation risk metrics VaR and CVaR for multidimensional time series using Archimedian copulas and stohastic processes as copula parameter.
Launch order:
1. Choose copula type and dimension of data and generate pdf function using
(currently no parameters are available. Calculated data/test.csv dataset with dim = 9 and Gumbel Copula)
./generate_pdf.sh
2. compile cpp code using
./compile.sh
3. run the code using
./run.sh

The program will calculate risk metrics in running window of lenght set in main.cpp (you can set here desirable parameters) and saves the result in risk_data/cvar_data.csv file.
First column is CVaR, Second VaR.

Available copula methods:
1. mle
2. scar-p-ou
3. scar-m-ou

Available marginals methods:
1. normal

Please make sure that you have installed libraries nlopt and eigen. The paths of this libs are set in CMakeLists.txt.
