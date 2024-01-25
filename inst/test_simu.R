rm(list=ls())
library(mvPoisVA)
library(susiF.alpha)
library(susieR)
data(N3finemapping)
X <- N3finemapping$X
mysd=0.1
N =30
genotype <-X[1:N,1:500]

idx <- which( apply( genotype,2, var ) <1e-15)
genotype <- genotype [, -idx]
library(gplots)#

if(file.exists("D:/Document/Serieux/Travail/Package/mvPoisVA/comp_fsusie.RData")){
  load("D:/Document/Serieux/Travail/Package/mvPoisVA/comp_fsusie.RData")

}else{
  res <-list()
}
lev_res =6
for (o  in (length(res)+1):10000) {


  count.data  <- list()
  L <-2# sample(1:2, size =1)#actual number of effect

  lf <-  list()
  for(l in 1:L){
    lf[[l]] <-log(abs(0.2*sim_intenisty(lev_res )$sim_intens) )#functional effect for effect l
  }




  data(N3finemapping)
  X <- N3finemapping$X
  genotype <-X[sample(1:nrow(X), size=N),]

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
  for ( i in 1:N)
  {

    predictor <-rep (0, length(lf[[1]] ))

    for ( l in 1:L){
      predictor <-predictor + G[i, true_pos[l]]*lf[[l]]
    }
    predictor <- exp(predictor+ rnorm(  length(lf[[1]]), sd=mysd))

    count.data [[i]] <-   rpois(n= length(lf[[1]]) ,
                                lambda =predictor  )

  }
  count.data <- do.call(rbind, count.data)


  Y <- count.data


  res3 <- acc_Pois_fSuSiE(Y=Y,X=X, L=3)
  res0 <-susiF (Y=log(Y+1),X=X, L=3)
  #res1 <- mv_Poisproc_reg (Y=Y,X=X, L=3)

  res2 <- Pois_fSuSiE (Y=Y,X=X, L=3)
  out <-  list( #mv_POIS = res1$susiF.obj$pip,
                #mv_POIS_cs= res1$susiF.obj$cs,
                susiF_pip= res0$pip,
                susiF_cs= res0$cs,
                g_susiF=res2$susiF.obj$pip,
                g_susiF_cs=res2$susiF.obj$cs,
                acc_g_susiF=res3$susiF.obj$pip,
                acc_g_susiF_cs=res3$susiF.obj$cs,
                true_pos=true_pos)
  res[[o]] <- out
  save(res, file = "D:/Document/Serieux/Travail/Package/mvPoisVA/comp_fsusie.RData")
}



tt <- res[[1]]
pip_acc <- c()
pip     <- c()
pip_acc <-  c(pip_acc, tt$acc_g_susiF)
pip  <-  c(pip , tt$g_susiF)

col <- rep("black", length(tt$g_susiF))
my_pch <- rep( 1, length(tt$g_susiF))

col[tt$true_pos] <- "red"
my_pch[tt$true_pos] <- 19
plot( pip, pip_acc,col=col, pch=my_pch)
for ( i in 1:38){
  tt <- res[[i]]
  col <- rep("black", length(tt$g_susiF))
  my_pch <- rep( 1, length(tt$g_susiF))

  col[tt$true_pos] <- "red"
  my_pch[tt$true_pos] <- 19
  points(   tt$acc_g_susiF,tt$g_susiF,col=col, pch=my_pch)
  pip_acc <-  c(pip_acc, tt$acc_g_susiF)
  pip  <-  c(pip , tt$g_susiF)

}

