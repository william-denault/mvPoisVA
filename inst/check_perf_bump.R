rm(list = ls())
library(mvPoisVA)
library(fsusieR)
library(susieR)
library(ebnm)
data(N3finemapping)
X <- N3finemapping$X
mysd=0.2
N =50
genotype <-X[1:N,1:100]

idx <- which( apply( genotype,2, var ) <1e-15)
genotype <- genotype [, -idx]
lev_res =6
count.data  <- list()
L <-2# sample(1:2, size =1)#actual number of effect

lf <-  list()
 lf[[1]]<- rep(0.1, 2^6)
 lf[[1]][10:20] <-2

 lf[[2]]<- rep(0.1, 2^6)
 lf[[2]][50:60] <-2



 mu= lf[[2]]
 x = rpois(length(mu),exp(log(mu)+rnorm(n=length(mu),sd=mysd)))

 fit = pois_smooth_split(x,maxiter=30)
 plot(x,col='grey80')
 lines(exp(fit$Eb))
 lines( smashr::smash(haarfisz::hft(x))  , col="green")
 lines(lf[[2]])

data(N3finemapping)
X <- N3finemapping$X
genotype <-X[sample(1:nrow(X), size=N),1:100]

idx <- which( apply( genotype,2, var ) <1e-15)
if( length(idx)==0){
  X <-genotype

  Rtrue <- cor (genotype )
}else{
  genotype <- genotype [, -idx]
  X <-genotype

}
G<- genotype
X <- (X -0.99*min(X))/(0.5*max(X ))

G <-  (G -0.99*min(G ))/(0.5*max(G ))

tpos <- sample(1:ncol(genotype), replace = FALSE,size=2)
true_pos <- tpos
pos1 <- tpos[1]
pos2 <- tpos[2]
if( length(which(apply(G,2,var)==0))>0){
  G <- G[,-which(apply(G,2,var)==0)]
}
# G <- matrix( rnorm(nrow(genotype)*300), nrow = nrow(genotype))


predictor <-rep (0, length(lf[[1]] ))
count.data  <- list()

G[ , true_pos[1]] <-G[ , true_pos[1]] -min(G[ , true_pos[1]] )
G[ , true_pos[2]] <-G[ , true_pos[2]] -min(G[ , true_pos[2]] )
for ( i in 1:N)
{

  predictor <-rep (0, length(lf[[1]] ))

  for ( l in 1:L){
    predictor <-predictor + G[i, true_pos[l]]*lf[[l]]+0.3
  }
  predictor <- exp( predictor + rnorm(  length(lf[[1]]), sd=mysd))

  count.data [[i]] <-   rpois(n= length(lf[[1]]) ,
                              lambda =predictor  )

}
count.data <- do.call(rbind, count.data)


Y <- count.data

res3<-  acc_Pois_fSuSiE2 (Y=Y,X=X, L=3 , post_processing = "HMM")

res31 <- acc_Pois_fSuSiE2 (Y=Y,X=X, L=3 , post_processing = "HMM", cal_obj.fsusie = TRUE)
res3$susiF.obj$cs
res31$susiF.obj$cs


true_pos


res0 <-susiF (Y=log(Y+1),X=X, L=3, post_processing = "HMM")
res0$cs
res01 <-HF_susiF2  (Y=Y,X=X, L=3 , post_processing = "HMM")
res01$cs

Y_t  <-fit_latent_space(Y )



Y_t1  <-fit_latent_nugget(Y )

plot(Y_t $Y, log(Y+1) )
points(log(Y_t1 $Y), log(Y+1), col="green")
abline(a=0,b=1)


res02 <-susiF (Y=Y_t$Y,X=X, L=3, post_processing = "HMM")
res02$cs

plot ( exp( res3$susiF.obj$fitted_func[[1]]  )-1)
lines(res3$susiF.obj$lfsr_func[[1]])
lines ( exp( res31$susiF.obj$fitted_func[[1]]  )-1)
lines(res31$susiF.obj$lfsr_func[[1]])
plot ( exp( res3$susiF.obj$fitted_func[[2]]  )-1)
lines(res3$susiF.obj$lfsr_func[[2]])
lines ( exp( res31$susiF.obj$fitted_func[[2]]  )-1)
lines(res31$susiF.obj$lfsr_func[[2]])


plot ( exp( res3$susiF.obj$fitted_func[[2]]  )-1)
lines(res3$susiF.obj$lfsr_func[[2]])

plot(lf[[1]], ylim= c( min (lf[[1]])-1, max(lf[[1]])+1) )
lines(res01$fitted_func[[1]], col='blue')
lines ( exp(  res3$susiF.obj$fitted_func[[1]]/ sqrt(var(X[, true_pos[1]])   ))-1, col='red')
lines ( res0$fitted_func[[1]], col='green')
lines ( exp(res02$fitted_func[[1]])-1, col='lightblue')

lines (exp( res31$susiF.obj$fitted_func[[1]])-1, col='blue4')


plot(lf[[2]], ylim= c( min (lf[[2]])-1, max(lf[[2]])+1) )
lines(res01$fitted_func[[2]], col='blue')
lines ( exp(  res3$susiF.obj$fitted_func[[2]]/ sqrt(var(X[, true_pos[2]])   ))-1, col='red')
lines ( exp( res0$fitted_func[[2]])-1 , col='green')
lines (exp( res02$fitted_func[[2]])-1, col='lightblue')

lines (exp( res31$susiF.obj$fitted_func[[2]])-1, col='blue4')


