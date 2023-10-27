rm(list=ls())
library(mvPoisVA)
library(susiF.alpha)
library(susieR)
data(N3finemapping)
X <- N3finemapping$X

N =50
genotype <-X[1:N,]

idx <- which( apply( genotype,2, var ) <1e-15)
genotype <- genotype [, -idx]
library(gplots)#

if(file.exists("/home/wdenault/benchmark_mvPois/sim/comparison_mvPois_fusie_gFSuSiE.RData")){
  load("/home/wdenault/benchmark_mvPois/comparison_mvPois_fusie_gFSuSiE.RData")

}else{
  res <-list()
}
lev_res =7
for (o  in (length(res)+1):10000) {


  beta0       <- 0
  beta1       <- 1
  beta2       <- 1
  count.data  <- list()
  L <-2# sample(1:2, size =1)#actual number of effect

  f1 <- abs(0.2*sim_intenisty(lev_res )$sim_intens[-1])#first effect
  f2 <- abs(0.2*sim_intenisty(lev_res )$sim_intens[-1])#second effect

  lf <-  list()
  for(l in 1:L){
    lf[[l]] <- simu_IBSS_per_level(lev_res=8)$sim_func #functional effect for effect l
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
  X <- X+5

  G <- genotype+5

  tpos <- sample(1:ncol(genotype), replace = FALSE,size=2)
  true_pos <- tpos
  pos1 <- tpos[1]
  pos2 <- tpos[2]
  if( length(which(apply(G,2,var)==0))>0){
    G <- G[,-which(apply(G,2,var)==0)]
  }
  # G <- matrix( rnorm(nrow(genotype)*300), nrow = nrow(genotype))


  beta0       <- 0
  beta1       <- 1
  beta2       <- 1
  count.data  <- list()
  for ( i in 1:N)
  {
    predictor <-beta1*G[i,pos1]*f1 + beta2*G[i,pos2]*f2

    count.data [[i]] <-   rpois(n= length(f1) ,
                                lambda =predictor  )

  }
  count.data <- do.call(rbind, count.data)


  Y <- count.data
  res0 <-susiF (Y=Y,X=X, L=3)
  res <- mv_Poisproc_reg (Y=Y,X=X, L=3)

  res2 <- Pois_fSuSiE (Y=Y,X=X, L=3)

  out <-  list( mv_POIS = res$susiF.obj$pip,
                mv_POIS_cs= res$susiF.obj$cs,
                susiF_pip= res0$pip,
                susiF_cs= res0$cs,
                g_susiF=res2$susiF.obj$pip,
                g_susiF_cs=res2$susiF.obj$cs,
                true_pos=true_pos)
  res[[o]] <- out
  save(res, file="/home/wdenault/benchmark_mvPois/comparison_mvPois_fusie_gFSuSiE.RData")

}
