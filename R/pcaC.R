pcaC <- function(mats) {
	.Call("pca", mats, PACKAGE="bayesic")
}

treat <- function(mats, pcas) {
	.Call("treat", mats, pcas, PACKAGE="bayesic")
}

pca <- function(mat.list, nmat, nsample) {
	pc <- matrix(ncol=nsample,nrow=ncol(mat.list[[1]]))
	for(i in 1:nsample) {
		aux <- sample(length(mat.list), nmat)
		mats <- mat.list[aux]
		mat <- do.call(rbind, mats)
		pc[,i] <- pcaC(mat)[,1]
	}
	pc1 <- apply(pc,1,mean)
}

treat.mat <- function(mat, pc1) {
	lcm <- apply(mat, 2, mean)
	mat.star <- mat
	for(i in 1:nrow(mat.star)) {
		mat.star[i,] <- mat.star[i,] - lcm
	}
	mat.star <- normalize(mat.star,byrow=T)

	res <- apply(mat.star, 1, function(y) {
		y - y%*%pc1%*%pc1
	})
	t(res)
}

treat.data <- function(mats, sample.size = NULL, pc1 = NULL) {

	if(pc1 == NULL) {
		## calculo do pca
		aux <- sample(length(mats),sample.size)
		mats.sample <- mats[aux]
		big.mat <- do.call(rbind,mats.sample)
		big.mat.star <- exp(big.mat)
		
		# w <- scale(big.mat.star, scale=FALSE)
		cm <- apply(big.mat.star, 2, mean)
		w <- big.mat.star
		for(i in 1:nrow(w)) {
			w[i,] <- w[i,]/cm
		}
		
		lw <- log(w)
		lcm <- log(cm)
		lw <- normalize(lw,byrow=T)

		pcas <- pcaC(lw)
		pc1 <- pcas[,1]
	}

	lapply(mats, function(x) {
		# x.star <- scale(x,center=lcm,scale=FALSE)
		lcm <- apply(x, 2, mean)
		x.star <- x
		for(i in 1:nrow(x.star)) {
			x.star[i,] <- x.star[i,] - lcm
		}
		x.star <- normalize(x.star,byrow=T)

		res <- apply(x.star, 1, function(y) {
			y - y%*%pc1%*%pc1
		})

		t(res)
	})

}