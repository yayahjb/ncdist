library("Rcpp")
dyn.load("rcpp_s6dist.so")
dyn.load("rcpp_cs6dist.so")
dyn.load("rcpp_ncdist.so")
dyn.load("rcpp_d7dist.so")
dyn.load("rcpp_cs6dist_in_g6.so")

args <- commandArgs(trailingOnly=TRUE)
if (length(args)<14) {
  stop("At least 14  arguments must be supplied (2 cells).\n", call.=FALSE)
}
latsym1 <- args[1]
a1 <- as.numeric(args[2])
b1 <- as.numeric(args[3])
c1 <- as.numeric(args[4])
alpha1 <- as.numeric(args[5])
beta1 <- as.numeric(args[6])
gamma1 <- as.numeric(args[7])
latsym2 <- args[8]
a2 <- as.numeric(args[9])
b2 <- as.numeric(args[10])
c2 <- as.numeric(args[11])
alpha2 <- as.numeric(args[12])
beta2 <- as.numeric(args[13])
gamma2 <- as.numeric(args[14])

cat("cell1: [",latsym1,a1,b1,c1,alpha1,beta1,gamma1,"] ")
cat("cell2: [",latsym2,a2,b2,c2,alpha2,beta2,gamma2,"] ")

cat("s6dist: ",.Call("rcpp_s6dist", latsym1,a1,b1,c1,alpha1,beta1,gamma1, latsym2,a2,b2,c2,alpha2,beta2,gamma2)," ")
cat("cs6dist: ",.Call("rcpp_cs6dist", latsym1,a1,b1,c1,alpha1,beta1,gamma1, latsym2,a2,b2,c2,alpha2,beta2,gamma2)," ")
cat("cs6dist_in_g6: ",.Call("rcpp_cs6dist_in_g6", latsym1,a1,b1,c1,alpha1,beta1,gamma1, latsym2,a2,b2,c2,alpha2,beta2,gamma2)," ")
cat("ncdist: ",.Call("rcpp_ncdist", latsym1,a1,b1,c1,alpha1,beta1,gamma1, latsym2,a2,b2,c2,alpha2,beta2,gamma2)," ")
cat("d7dist: ",.Call("rcpp_d7dist", latsym1,a1,b1,c1,alpha1,beta1,gamma1, latsym2,a2,b2,c2,alpha2,beta2,gamma2),"\n")
