rm(list = ls())
library(mvPoisVA)
library(susiF.alpha)
library(susieR)
library(ebnm)

'%!in%' <- function(x,y)!('%in%'(x,y))
data(N3finemapping)
X <- N3finemapping$X
mysd=2
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
  predictor <- exp( predictor+ rnorm(  length(lf[[1]]), sd=mysd))

  count.data [[i]] <-   rpois(n= length(lf[[1]]) ,
                              lambda =predictor  )

}
count.data <- do.call(rbind, count.data)


Y <- count.data

false_disc1 <- NULL
true_disc1  <- NULL
false_disc2 <- NULL
true_disc2  <- NULL
res01 <-acc_Pois_fSuSiE2 (Y=Y,X=X, L=3 , post_processing = "HMM")

res01 <-acc_Pois_fSuSiE2 (Y=Y,X=X, L=3 , post_processing = "HMM")
false_disc1 <- NULL
true_disc1  <- NULL
false_disc2 <- NULL
true_disc2  <- NULL
if ( true_pos[1] %in% res01$susiF.obj$cs[[1]]){

  no_zero1 <-  which(res01$susiF.obj$lfsr_func[[1]] <0.05
  )
  if (length(no_zero1)>0){
    false_disc1 <- length( which ( no_zero1 %!in% which(lf[[1]]>0)))

    true_disc1 <- length(which ( no_zero1 %in% which(lf[[1]]>0)))
  }

}
if ( length(res01$susiF.obj$cs)>1){
  if ( true_pos[1] %in% res01$susiF.obj$cs[[2]]){

    no_zero1 <-  which(res01$susiF.obj$lfsr_func[[1]] <0.05
    )
    if (length(no_zero1)>0){
      false_disc1 <- length( which ( no_zero1 %!in% which(lf[[1]]>0)))

      true_disc1 <- length(which ( no_zero1 %in% which(lf[[1]]>0)))
    }

  }
}



if ( true_pos[2] %in% res01$susiF.obj$cs[[1]]){

  no_zero2 <-  which(res01$susiF.obj$lfsr_func[[1]] <0.05
  )
  if (length(no_zer2o)>0){
    false_disc2 <- length( which ( no_zero2 %!in% which(lf[[1]]>0)))

    true_disc2 <- length(which ( no_zero2 %in% which(lf[[1]]>0)))
  }

}
if ( length(res01$susiF.obj$cs)>1){
  if ( true_pos[2] %in% res01$susiF.obj$cs[[2]]){

    no_zero2 <-  which(res01$susiF.obj$lfsr_func[[1]] <0.05
    )
    if (length(no_zero2)>0){
      false_disc2 <- length( which ( no_zero2 %!in% which(lf[[1]]>0)))

      true_disc2 <- length(which ( no_zero2 %in% which(lf[[1]]>0)))
    }

  }
}
