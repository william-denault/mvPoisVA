rm(list = ls())
library(mvPoisVA)
library(fsusieR)
library(susieR)
library(ebnm)


lf1_fit <-  list()
lf2_fit <-  list()

for (o in 1:400){
print(o)
  data(N3finemapping)
  X <- N3finemapping$X
  mysd=0.5
  N =200
  genotype <-X[1:N,1:100]

  idx <- which( apply( genotype,2, var ) <1e-15)
  genotype <- scale( genotype [, -idx])
  lev_res =6

  count.data  <- list()
  L <-2# sample(1:2, size =1)#actual number of effect

  lf <-  list()
  lf[[1]]<- rep(0.1, 2^6)
  lf[[1]][10:20] <-5

  lf[[2]]<- rep(0.1, 2^6)
  lf[[2]][50:60] <-10



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

  res3$susiF.obj$cs



  lf1_fit[[o]] <-  exp( res3$susiF.obj$fitted_func[[1]]  )-1
  lf2_fit[[o]] <-  exp( res3$susiF.obj$fitted_func[[2]]  )-1


}

lf1_out <- do.call(rbind, lf1_fit)

lf2_out <- do.call(rbind, lf2_fit)

lf_out = lf1_out + lf2_out


lf_out = log(lf_out+1)
plot(lf[[1]]+lf[[2]], ylim= c( -0.2, 12) , col="red")

lines(apply(lf_out, 2, mean) +0.1, col="blue")

for (o in 1: length(lf1_fit)){
  lines (lf1_fit[[o]], col="blue")
}

true_pos

plot(lf[[2]], ylim= c( -0.2, 5) , col="red")

for (o in 1: length(lf2_fit)){
  lines (lf2_fit[[o]], col="blue")
}

