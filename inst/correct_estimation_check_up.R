rm(list = ls())
library(mvPoisVA)
library(fsusieR)
library(susieR)
library(ebnm)
'%!in%' <- function(x,y)!('%in%'(x,y))
data(N3finemapping)
X <- N3finemapping$X
mysd=1.1
N =150
lev_res=6
genotype <-X[1:N,1:100]
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




idx <- which( apply( genotype,2, var ) <1e-15)

if ( length(idx)>0){
  genotype <- genotype [, -idx]
}

count.data  <- list()
L <-2# sample(1:2, size =1)#actual number of effect

lf <-  list()
lf[[1]]<- 2*cos((1:2^lev_res) /(0.5*2^lev_res))# #rep(0.1, 2^lev_res)
#lf[[1]][10:20] <-2

lf[[2]]<-2*sin((1:2^lev_res) /(0.5*2^lev_res)) # rep(0.1, 2^lev_res)
#lf[[2]][50:60] <-2

true_pos <- sample(1:ncol(genotype), replace = FALSE,size=2)

plot (lf[[1]], type="l", main=paste ( "effect of SNP",true_pos[1]))

plot (lf[[2]], type="l", main=paste ( "effect of SNP",true_pos[2]))


pos1 <- true_pos[1]
pos2 <- true_pos[2]
if( length(which(apply(G,2,var)==0))>0){
  G <- G[,-which(apply(G,2,var)==0)]
}
# G <- matrix( rnorm(nrow(genotype)*300), nrow = nrow(genotype))
scaling= runif( n= N, min=0.8,max=1.3)

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
  predictor <- scaling[i] *exp( predictor  + rnorm(  length(lf[[1]]), sd=mysd))

  count.data [[i]] <-   rpois(n= length(lf[[1]]) ,
                              lambda =predictor  )

}
count.data <- do.call(rbind, count.data)


Y <- count.data

#'pois_mean_split',
#'ind_pois_mean_split',
#'ind_ebps',
#'ind_poisson_smoothing',
#'nugget'),

print( true_pos)
res01 <-acc_Pois_fSuSiE2 (Y=Y,X=X, L=3, post_processing = "TI"  ,scaling=scaling,
                          ebps_method='pois_mean_split')
res01$susiF.obj$cs


plot ( res01$susiF.obj$fitted_func[[2]]  , type="l", main=paste ( "effect of SNP",true_pos[2]), col= "green4")
lines(lf[[2]])

plot ( res01$susiF.obj$fitted_func[[1]]    , type="l", main=paste ( "effect of SNP",true_pos[2]), col= "green4")
lines(lf[[1]] )


res02 <-acc_Pois_fSuSiE2 (Y=Y,X=X, L=3, post_processing = "TI"  ,scaling=scaling,
                          ebps_method='ind_pois_mean_split')
res02$susiF.obj$cs





plot ( res02$susiF.obj$fitted_func[[2]]  , type="l", main=paste ( "effect of SNP",true_pos[2]), col= "green4")
lines(lf[[2]])

plot ( res02$susiF.obj$fitted_func[[1]]    , type="l", main=paste ( "effect of SNP",true_pos[2]), col= "green4")
lines(lf[[1]] )




res03 <-acc_Pois_fSuSiE2 (Y=Y,X=X, L=3, post_processing = "TI" ,  scaling=scaling,
                          ebps_method='ind_ebps')
res03$susiF.obj$cs


plot ( res03$susiF.obj$fitted_func[[2]]  , type="l", main=paste ( "effect of SNP",true_pos[2]), col= "green4")
lines(lf[[2]])

plot ( res03$susiF.obj$fitted_func[[1]]    , type="l", main=paste ( "effect of SNP",true_pos[2]), col= "green4")
lines(lf[[1]] )




res04 <-acc_Pois_fSuSiE2 (Y=Y,X=X, L=3, post_processing = "TI" ,scaling=scaling,
                          ebps_method='ind_poisson_smoothing')
res04$susiF.obj$cs

res00= fsusieR::susiF(Y=log1p( Y/scaling),X=X, L=3, post_processing = "TI" )
res00$cs

true_pos
plot( lf[[2]])
lines(res03$susiF.obj$fitted_func[[2]], col="red")
lines(res03$susiF.obj$cred_band[[2]][1,], col="red", lty =2)
lines(res03$susiF.obj$cred_band[[2]][2,], col="red", lty =2)


lines(res01$susiF.obj$fitted_func[[2]], col="blue")
lines(res01$susiF.obj$cred_band[[2]][1,], col="blue", lty =2)
lines(res01$susiF.obj$cred_band[[2]][2,], col="blue", lty =2)

lines(res04$susiF.obj$fitted_func[[2]], col="green")
lines(res04$susiF.obj$cred_band[[2]][1,], col="green", lty =2)
lines(res04$susiF.obj$cred_band[[2]][2,], col="green", lty =2)

lines(res00$fitted_func[[2]] )
lines(res00$ cred_band[[2]][1,] , lty =2)
lines(res00$ cred_band[[2]][2,] , lty =2)

plot( lf[[1]])
lines(res03$susiF.obj$fitted_func[[1]], col="red")
lines(res03$susiF.obj$cred_band[[1]][1,], col="red", lty =2)
lines(res03$susiF.obj$cred_band[[1]][2,], col="red", lty =2)


lines(res01$susiF.obj$fitted_func[[1]], col="blue")
lines(res01$susiF.obj$cred_band[[1]][1,], col="blue", lty =2)
lines(res01$susiF.obj$cred_band[[1]][2,], col="blue", lty =2)

lines(res04$susiF.obj$fitted_func[[1]], col="green")
lines(res04$susiF.obj$cred_band[[1]][1,], col="green", lty =2)
lines(res04$susiF.obj$cred_band[[1]][2,], col="green", lty =2)

lines(res00$fitted_func[[1]] )
lines(res00$ cred_band[[1]][1,] , lty =2)
lines(res00$ cred_band[[1]][2,] , lty =2)
sqrt(res01$susiF.obj$sigma2)
sqrt(res03$susiF.obj$sigma2)
sqrt(res04$susiF.obj$sigma2)
sqrt(res00$sigma2)
mysd
