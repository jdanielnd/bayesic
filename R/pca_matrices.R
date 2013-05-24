pca_matrices <- function(sample_size, dir="/home/tecto/UFMG/Xk_Matrices") {
  all_mat_names <- list.files(dir, recursive=TRUE)
  if(sample_size > length(all_mat_names)) stop("sample_size is too big")
  sample_mat_names <- sample(all_mat_names, sample_size)
  lmats <- lapply(sample_mat_names, function(x) {
    read.table(paste(dir,x,sep="/"))
  })
  big_mat <- do.call(rbind, lmats)
  big_mat_star <- exp(big_mat)
  
  cm <- apply(big_mat_star, 2, mean)
  w <- big_mat_star
  for(i in 1:nrow(w)) {
    w[i,] <- w[i,]/cm
  }
  
  lw <- log(w)
  lcm <- log(cm)
  lw <- normalize(lw,byrow=T)

  pcas <- pcaC(lw)
  pc1 <- pcas[,1]
}