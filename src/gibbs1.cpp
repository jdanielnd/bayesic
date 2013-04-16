#include <Rcpp.h>

RcppExport SEXP gibbs1(SEXP ite, SEXP dx, SEXP dy) {  
  try {
    
    int dite = Rcpp::as<int>(ite);
    double a_star, b_sum, b_star, sig_0, mu_0sum, mu_0, sig_1sum, sig_1, mu_1sum, mu_1, sigma;
    double beta_0 = 0;
    double beta_1 = 0;
    double a = 2.1;
    double b = 1.1;
    double tau_0 = 5;
    double tau_1 = 5;
    Rcpp::NumericVector x(dx), y(dy);
    Rcpp::NumericMatrix mat(dite,3);
    int n = x.size();
    
    for(int i=0; i<dite; ++i) {
      a_star = a + n/2;
      b_sum = 0;
      for(int j=0; i<n; ++i) {
        b_sum += pow((y[i] - (beta_0 + beta_1*x[i])),2);
      }
      b_star = b + b_sum/2;
      
      sigma = 1/Rcpp::as<double>(Rcpp::rgamma(1, a_star, b_star));
      
      sig_0 = (sigma*tau_0)/(n*tau_0+sigma);
      for(int j=0; i<n; ++i) {
        mu_0sum += (y[i]*x[i] - x[i]*beta_0);
      }
      mu_0 = sig_0*mu_0sum/sigma;
      
      beta_0 = Rcpp::as<double>(Rcpp::rnorm(1,mu_0,sqrt(sig_0)));
      
      for(int j=0; i<n; ++i) {
        sig_1sum += pow(x[i],2);
      }
      sig_1 = (sigma*tau_1)/(tau_1*sig_1sum+sigma);
      for(int j=0; i<n; ++i) {
        mu_1sum += y[i]*x[i] - x[i]*beta_0;
      }
      mu_1 = sig_1*mu_1sum/sigma;
      
      beta_1 = Rcpp::as<double>(Rcpp::rnorm(1, mu_1, sqrt(sig_1)));
      
      mat(i,0) = beta_0;
      mat(i,1) = beta_1;
      mat(i,2) = sigma;
    }
    
    return Rcpp::wrap(mat);
    
  } catch(...) {
    ::Rf_error("c++ error");
  }
}
