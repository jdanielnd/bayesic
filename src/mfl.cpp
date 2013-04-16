#include <RcppArmadillo.h>

RcppExport SEXP mfl(SEXP ite, SEXP sx) {  

	int dite = Rcpp::as<int>(ite);
	double omega = 10, a = 2.1, b = 1.1, gamma_a = 1, gamma_b = 1;

	Rcpp::NumericMatrix rx(sx);
	int m = rx.nrow(), n = rx.ncol();
	arma::mat x(rx.begin(), m, n, false);
	arma::mat alpha_mat(dite,m);
	arma::mat lambda_mat(dite,n);
	arma::mat v_lambda_mat(dite,n);
	arma::mat sigma_mat(dite,m);
	arma::mat z_mat(dite,m);
	arma::mat p_mat(dite,m);
	arma::mat dinv(m,m);
	double A, B, gamma_1, gamma_2, p = 0.5, v_alpha, m_alpha, v_lambda, m_lambda;
	Rcpp::NumericVector sigma(m);
	arma::vec alpha(m);
	arma::rowvec lambda = Rcpp::as<arma::rowvec>(Rcpp::rnorm(n,0,0.1));
	Rcpp::NumericVector p_star(m, 0.5);
	// criar vetor de class integer
	// Rcpp::NumericVector z(m);
	std::vector<int> z(m,0);

	Rcpp::RNGScope scope;

try{

	for(int k=0; k<dite; k++) {
	  
	  for(int i=0; i<m; i++) {
	  	A = a + n/2;
	  	// double b_1=0, b_2=0, b_3=0;
	  	// for(int j=0; j<n; j++) {
	  	// 	b_1 = b_1 + pow(x(i,j),2);
	  	// 	b_2 = b_2 + lambda[j] * x(i,j);
	  	// 	b_3 = b_3 + pow(lambda[j],2);
	  	// }
	    double b_1 = arma::as_scalar(x.row(i)*x.row(i).t());
	    double b_2 = arma::as_scalar(lambda*x.row(i).t());
	    double b_3 = arma::as_scalar(lambda*lambda.t());
	    B = b + 0.5*(b_1 - 2*b_2*alpha[i] +alpha[i]*b_3*alpha[i]);
	    if (A<0) {throw std::runtime_error("A cannot be negative");}
	    if (B<0) {throw std::runtime_error("B cannot be negative");}
	    sigma[i] = 1/::Rf_rgamma(A, 1/B);

	    // verificar valores diferentes de 0 e 1
	    if (::Rf_runif(0,1) < p_star[i]) {
	    	z[i] = 1;
	    } else {
	    	z[i] = 0;
	    }
	    gamma_1 = gamma_a + z[i];
	    gamma_2 = gamma_b + 1 - z[i];
	   	if (gamma_1<0) {throw std::runtime_error("gamma_1 cannot be negative");}
	   	if (gamma_2<0) {throw std::runtime_error("gamma_2 cannot be negative");}
	    p = ::Rf_rbeta(gamma_1, gamma_2);
	    
	    Rcpp::NumericVector aux_1(Rcpp::wrap(arma::pow(lambda,2)));
	    v_alpha = 1/((1/omega) + (Rcpp::sum(aux_1)/sigma[i]));
	    Rcpp::NumericVector aux_2(Rcpp::wrap(x.row(i)*lambda.t()));
	    m_alpha = v_alpha * (Rcpp::sum(aux_2)/sigma[i]);
		if (v_alpha<0) {throw std::runtime_error("v_alpha cannot be negative");}

	    if(z[i]==1) {
	      alpha[i] = ::Rf_rnorm(m_alpha,sqrt(v_alpha));
	    } else {
	      alpha[i] = 0;
	    }
	    double p_1 = Rcpp::dnorm(Rcpp::NumericVector::create(1),m_alpha,sqrt(v_alpha))[0];
	    double p_2 = Rcpp::dnorm(Rcpp::NumericVector::create(1),0,sqrt(omega))[0];
	    p_star[i] = p/(p + p_1/(p_2*(1-p)));

	    alpha_mat(k,i) = alpha[i];
	    sigma_mat(k,i) = sigma[i];
	    z_mat(k,i) = z[i];
	    p_mat(k,i) = p_star[i];
	  	dinv(i,i) = 1/sigma[i];
	  	
	    if (isnan(sigma[i])) {throw std::runtime_error("sigma should not be NaN");}
	    if (isnan(alpha[i])) {throw std::runtime_error("alpha should not be NaN");}
	    if (isnan(p_star[i])) {throw std::runtime_error("p_star should not be NaN");}
	  }

	  for(int j=0; j<n; j++) {

	  	// lambda NaN
	  	// a conta abaixo pode ser substituida por sum(alpha[i]^2*sigma[i]) verificar sigma, sigma2
	  	v_lambda = 1/arma::as_scalar(alpha.t()*dinv*alpha + 1);
	  	// double vl = 0;
	  	// for(int i; i<m; i++) {
	  	// 	vl = vl + pow(alpha[i],2)/sigma[i];
	  	// }
	  	// v_lambda = 1/(vl + 1);
	  	// v_lambda_mat(k,j) = v_lambda;
	  	// sum((alpha[i]*x[i,j])/sigma2[i])
	  	m_lambda = v_lambda*arma::as_scalar(alpha.t()*dinv*x.col(j));
	  	
	  	// double ml = 0;
	  	// for(int i; i<m; i++){
	  	// 	ml = ml + (alpha[i]*x(i,j))/sigma[i];
	  	// }
	  	// m_lambda = v_lambda*ml;

	  	if (v_lambda<0) {throw std::runtime_error("v_lambda cannot be negative");}
	  	lambda[j] = ::Rf_rnorm(m_lambda,sqrt(v_lambda));
	  	lambda_mat(k,j) = lambda[j];
	  	if (isnan(lambda[j])) {throw std::runtime_error("lambda should not be NaN");}
	  }


	}

	return Rcpp::List::create(
        Rcpp::Named("alpha")  = Rcpp::wrap(alpha_mat),
        Rcpp::Named("lambda") = Rcpp::wrap(lambda_mat),
        Rcpp::Named("sigma")  = Rcpp::wrap(sigma_mat),
        Rcpp::Named("z")      = Rcpp::wrap(z_mat),
        Rcpp::Named("p")      = Rcpp::wrap(p_mat),
        Rcpp::Named("v_lambda") = Rcpp::wrap(v_lambda_mat)
    ) ;

    } catch(std::runtime_error e) {
		return Rcpp::List::create(
	        Rcpp::Named("alpha")  = Rcpp::wrap(alpha_mat),
	        Rcpp::Named("lambda") = Rcpp::wrap(lambda_mat),
	        Rcpp::Named("sigma")  = Rcpp::wrap(sigma_mat),
	        Rcpp::Named("z")      = Rcpp::wrap(z_mat),
	        Rcpp::Named("p")      = Rcpp::wrap(p_mat),
	        Rcpp::Named("v_lambda") = Rcpp::wrap(v_lambda_mat)
	    );
       ::Rf_error(e.what());
    }
}