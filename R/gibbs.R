## MODELO GIBBS - REGRESSAO LINEAR Y = A + BX - ULTIMA REVISAO: 23/01/2013
##
gibbsR <- function(ite, beta_0 = 0, beta_1 = 0, a=2.1, b=1.1, tau_0=5, tau_1=5, x, y) {

  mat <- matrix(nrow=ite,ncol=3)
  
  n <- length(x)
  
  for(i in 1:ite) {
    # Setting Precision Posteriori parameters a and b
    a_star <- a + n/2
    b_star <- b + sum((y - (beta_0 + beta_1*x))^2)/2
    
    # Generating precision from Gamma
    pre <- rgamma(1,a_star, rate=b_star)
    
    # Sigma2 from the inverse of precision
    sig <- 1/pre
    
    # Setting Beta0 Posteriori parameters mu_0 and sig_0
    sig_0 <- (sig*tau_0)/(n*tau_0+sig)
    mu_0 <- sig_0*(sum(y - x*beta_1))/sig
    
    # Generating Beta0 from Normal
    beta_0 <- rnorm(1, mu_0, sqrt(sig_0))
    
    # Setting Beta1 Posteriori parameters mu_1 and sig_1
    sig_1 <- (sig*tau_1)/(tau_1*sum(x^2)+sig)
    mu_1 <- sig_1*(sum(y*x - x*beta_0))/sig
    
    # Generating Beta1 from Normal
    beta_1 <- rnorm(1, mu_1, sqrt(sig_1))
    
    mat[i,] <- c(beta_0, beta_1, sig)
  }
  colnames(mat) <- c("beta_0","beta_1","sigma2")
  mat
}

gibbsC <- function(ite, x, y) {
  .Call("gibbs1", ite, x, y, PACKAGE="bayesic")
}
