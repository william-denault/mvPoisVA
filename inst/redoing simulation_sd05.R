
rm(list = ls())
library(mvPoisVA)
library(fsusieR)
library(susieR)
library(ebnm)
'%!in%' <- function(x,y)!('%in%'(x,y))
data(N3finemapping)
X <- N3finemapping$X
mysd=1
N =200

if(file.exists("C:/Document/Serieux/Travail/Package/mvPoisVA/simredoing_mv_sd_0.5.RData")){
  load("C:/Document/Serieux/Travail/Package/mvPoisVA/simredoing_mv_sd_0.5.RData")
}else{

  res_list=list()
}


for( o in (length(res_list)+1):300){
  set.seed(length(res_list)+1)

  genotype <-X
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






  idx <- which( apply( genotype,2, var ) <1e-15)

  if ( length(idx)>0){
    genotype <- genotype [, -idx]
  }
  lev_res =6
  count.data  <- list()
  L <- sample(1:5, size =1)#actual number of effect

  lf <-  list()
  for ( l in 1:L ){

    if ( l%%2==1){
      lf[[l]]<- cos((1:2^6) /(0.5*2^6))
    }else{

      lf[[l]]<- sin((1:2^6) /(0.5*2^6))
    }
    # #rep(0.1, 2^6)
    #lf[[1]][10:20] <-2

    # rep(0.1, 2^6)
    #lf[[2]][50:60] <-2
  }


  true_pos <- sample(1:ncol(genotype), replace = FALSE,size=L)


  if( length(which(apply(G,2,var)==0))>0){
    G <- G[,-which(apply(G,2,var)==0)]
  }
  # G <- matrix( rnorm(nrow(genotype)*300), nrow = nrow(genotype))


  predictor <-rep (0, length(lf[[1]] ))
  count.data  <- list()

  for ( l in 1:L){

    G[ , true_pos[l]] <-G[ , true_pos[l]] -min(G[ , true_pos[l]] )
  }

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
  res01 <-acc_Pois_fSuSiE2 (Y=Y,X=X, L=5, post_processing = "TI" ,
                            ebps_method='ind_ebps')
  m1= res01$susiF.obj


  m2 <-fsusieR::susiF (Y=log(Y+1),X=X, L=5, post_processing = "TI"  )

  cal_purity <- function(l_cs,X){
    tt <- list()
    for (k in 1:length(l_cs)){
      if(length(unlist(l_cs[[k]]))==1 ){
        tt[[k]] <- 1
      }else{
        x <-abs( cor(X[,unlist(l_cs[[k]]   ) ]))


        tt[[k]] <-  min( x[col(x) != row(x)])
      }
    }
    return( mean( unlist(tt )))
  }



  out <- c( length(m1$cs), #number of CS
            length(which(true_pos%in% do.call(c, m1$cs))), #number of effect found
            Reduce("+",sapply(1:length(m1$cs), function(k)
              ifelse( length(which(true_pos%in%m1$cs[[k]] ))==0, 1,0)
            )
            ),#number of CS without any effect
            cal_purity(m1$cs, X=as.matrix(G)),#mean purity
            mean(sapply( m1$cs, length)), #CS size
            length(m2$cs),
            length(which(true_pos%in% do.call(c, m2$cs))),
            Reduce("+",sapply(1:length(m2$cs), function(k)
              ifelse( length(which(true_pos%in%m2$cs[[k]] ))==0, 1,0)
            )
            ),#number of CS without any effect
            cal_purity(m2$cs, X=as.matrix(G)),#mean purity
            mean(sapply( m2$cs, length)), #CS size
            L)


  names(out)=c("ncs_pois",
               "true_cs_pois",
               "false_cs_pois",
               "purity_pois",
               "cs_len_pois",
               "ncs",
               "true_cs",
               "false_cs",
               "purity",
               "cs_len","L")

  res_list[[o]] <-unlist(out)

 # print(do.call(rbind, res_list))

  save(res_list, file="simredoing_mv_sd_1.RData")
}

