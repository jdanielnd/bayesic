#include <RcppArmadillo.h>

RcppExport SEXP pca(SEXP mats) {  
	try { 
  
		Rcpp::NumericMatrix matr(mats);
		int n = matr.nrow(), k = matr.ncol();
		arma::mat mat(matr.begin(), n, k, false);

		arma::mat pca = arma::princomp(mat);
		return Rcpp::wrap(pca);
  
	}  catch(...) {
		::Rf_error("C++ Error");
	}
}

arma::mat read_matrix(SEXP mats) {
	try {

		Rcpp::NumericMatrix matr(mats);
		int n = matr.nrow(), k = matr.ncol();
		arma::mat mat(matr.begin(), n, k, false);
		
		return mat;

	} catch(...) {
		::Rf_error("read_matrix error");
	}
}

arma::mat divide_mean(arma::mat mat) {
	try {
		int n = mat.n_rows, k = mat.n_cols;

		for(int i = 0; i < k; i++) {
			// Calculate columns mean
			double mean, total = 0;
			for(int j=0; j < n; j++) {
				total += mat(j,i);
			}
			mean = total/n;
			// Divide elements by column mean
			for(int j=0; j < n; j++) {
				mat(j,i) = mat(j,i)/mean;
			}
		}

		return mat;
	} catch(...) {
		::Rf_error("divide_mean error");
	}

}

arma::mat log_trans(arma::mat mat) {

	try{ 
		int n = mat.n_rows, k = mat.n_cols;
		for(int i = 0; i < k; i++) {
			for(int j=0; j < n; j++) {
				mat(j,i) = log(mat(j,i));
			}
		}

		return mat;
	} catch(...) {
		::Rf_error("log_trans error");
	}

}

arma::mat stand(arma::mat mat) {

	try{
		int n = mat.n_rows, k = mat.n_cols;

		for(int j = 0; j < n; j++) {
			// Calculate mean
			double mean, sd, total = 0, sum = 0;
			for(int i = 0; i < k; i++) {
				total += mat(j,i);
			}
			mean = total/n;
			// Subtract mean and calculate sd
			for(int i = 0; i < k; i++) {
				sum += (mat(j,i) - mean)*(mat(j,i) - mean);
				mat(j,i) = mat(j,i) - mean;
			}
			sd = sqrt(1/(k*sum));
			// Divide elements by sd
			for(int i = 0; i < k; i++) {
				mat(j,i) = mat(j,i)/sd;
			}
		}

		return mat;
	} catch(...) {
		::Rf_error("stand error");
	}
}

arma::mat pca_subt(arma::mat mat, arma::mat pca) {

	try {
		int n = mat.n_rows, k = mat.n_cols;


		for(int j = 0; j < n; j++) {
			mat.row(j) -= mat.row(j)*pca;
		}	
		return mat;
	} catch(...) {
		::Rf_error("pca_subt error");
	}

}

RcppExport SEXP treat(SEXP mats, SEXP pcas) {

	// Rcpp::NumericMatrix matr(mats);
	// int n = matr.nrow(), k = matr.ncol();
	// arma::mat mat(matr.begin(), n, k, false);
	arma::mat mat = read_matrix(mats);
	// mat.print("mat");

	// Rcpp::NumericMatrix pcar(pcas);
	// int pcan = pcar.nrow(), pcak = pcar.ncol();
	// arma::mat pca(pcar.begin(), pcan, pcak, false);
	arma::mat pca = read_matrix(pcas);
	// mat.print("pca");

	// Step 1
	// for(int i = 0; i < k; i++) {
	// 	// Calculate columns mean
	// 	double mean, total = 0;
	// 	for(int j=0; j < n; j++) {
	// 		total += mat(j,i);
	// 	}
	// 	mean = total/n;
	// 	// Divide elements by column mean
	// 	for(int j=0; j < n; j++) {
	// 		mat(j,i) = mat(j,i)/mean;
	// 	}
	// }
	// mat = divide_mean(mat);
	// mat.print("divide_mean");

	// Log transformation
	// for(int i = 0; i < k; i++) {
	// 	for(int j=0; j < n; j++) {
	// 		mat(j,i) = log(mat(j,i));
	// 	}
	// }
	mat = log_trans(mat);
	// mat.print("log_trans");

	// Standardization of rows
	// for(int j = 0; j < n; j++) {
	// 	// Calculate mean
	// 	double mean, sd, total = 0, sum = 0;
	// 	for(int i = 0; i < k; i++) {
	// 		total += mat(j,i);
	// 	}
	// 	mean = total/n;
	// 	// Subtract mean and calculate sd
	// 	for(int i = 0; i < k; i++) {
	// 		sum += (mat(j,i) - mean)*(mat(j,i) - mean);
	// 		mat(j,i) = mat(j,i) - mean;
	// 	}
	// 	sd = sqrt(1/(k*sum));
	// 	// Divide elements by sd
	// 	for(int i = 0; i < k; i++) {
	// 		mat(j,i) = mat(j,i)/sd;
	// 	}
	// }
	mat = stand(mat);
	// mat.print("stand");

	// for(int j = 0; j < n; j++) {
	// 	mat.row(j) -= mat.row(j)*pca;
	// }	
	mat = pca_subt(mat, pca);
	// mat.print("pca_subt");

	return Rcpp::wrap(mat);

}