rm(list=ls())
library(mvPoisVA)
library(susiF.alpha)
library(susieR)
data(N3finemapping)
mysd=1
N =100

X <- N3finemapping$X[1:N,]
X <- X[, -which(apply(X,2,var)==0)]
genotype <-X[1:N,1:500]
n=2^7
sigma=0.5
mu=c(rep(0.3,n/4), rep(3, n/4), rep(10, n/4), rep(0.3, n/4))[-1]
x = rpois(n,exp(log(mu)+rnorm(n,sd=sigma)))
idx <- which( apply( genotype,2, var ) <1e-15)
genotype <- genotype [, -idx]


count.data  <- list()
L <-  1#actual number of effect

lf <-  list()
for(l in 1:L){
  lf[[l]] <-mu#functional effect for effect l
}

true_pos <- sample (size=1, 1:ncol(X))
X[,true_pos] <- X[,true_pos]-min( X[,true_pos] )

noise_data <- list()
for ( i in 1:N){

  predictor <-   0.1*(DJ.EX(n=2^7)[[4]]) *X[i,true_pos] +rnorm(length(lf[[1]]), sd=.5)

  noise_data[[i]]<- rpois(n, lambda=exp (predictor) )
}

Y <- do.call(rbind, noise_data)
out <-acc_Pois_fSuSiE(Y=Y,X=X)
out2 <-susiF(Y=log(Y+1),X=X)
out$susiF.obj$cs
true_pos
out2$cs


plot( 0.051*(DJ.EX(n=2^7)[[4]]))
lines (out2$fitted_func[[1]])
lines(  2*out$susiF.obj$fitted_func[[1]],col="green" )

lines(mu)
plot ( out$susiF.obj$lfsr_func[[1]])
