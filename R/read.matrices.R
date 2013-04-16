read.matrices <- function(start=1, end=1000, base_dir = "/home/tecto/Dropbox/UFMG/IC/Xk_Matrices_1/") {
	lapply(seq(start,end), function(x) {
		as.matrix(read.table(paste(base_dir,"Xk_Probeset", sprintf("%05d",x) ,".txt",sep="")))
	})
}