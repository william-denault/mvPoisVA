rm(list=ls())
library(ashr)
library(mvPoisVA)
library(fsusieR)
library(susieR)
'%!in%' <- function(x,y)!('%in%'(x,y))


lev_res =7
o=2
set.seed(o+2)
lstfiles <- list.files("/home/wdenault/scratch-midway2/enhancer_fine_mapping/formated_data/")


thin_base <- function (mat, designmat, coefmat, relative = TRUE, type = c("thin",
                                                                          "mult"))
{
  assertthat::assert_that(is.matrix(mat))
  assertthat::assert_that(is.matrix(designmat))
  assertthat::assert_that(is.matrix(coefmat))
  assertthat::assert_that(is.numeric(mat))
  assertthat::assert_that(is.numeric(designmat))
  assertthat::assert_that(is.numeric(coefmat))
  assertthat::are_equal(nrow(mat), nrow(coefmat))
  assertthat::are_equal(ncol(mat), nrow(designmat))
  assertthat::are_equal(ncol(designmat), ncol(coefmat))
  assertthat::assert_that(is.logical(relative))
  assertthat::are_equal(1L, length(relative))
  stopifnot(mat >= 0)
  type <- match.arg(type)
  meanmat <- tcrossprod(coefmat, designmat)
  maxvec <- apply(meanmat, 1, max)
  if (!relative) {
    if (any(maxvec > 0)) {
      stop(paste0("thin_base: tcrossprod(coefmat, designmat) produced positive entries\n",
                  "       and relative = FALSE. Either set relative = TRUE or change your\n",
                  "       coefficient and design matrices."))
    }
    qmat <- 2^meanmat
  }
  else {
    qmat <- 2^(meanmat - maxvec)
  }
  if (type == "thin") {
    newmat <- stats::rbinom(n = prod(dim(mat)), size = mat,
                            prob = qmat)
  }
  else if (type == "mult") {
    newmat <- round(qmat * mat)
  }
  dim(newmat) <- dim(mat)
  return(newmat)
}

load("C:/Document/Serieux/Travail/Package/susiF.alpha/inst/data_RCC/ZFAND2A_locus_7.RData")

Y <- data$Y$Y_f[[1]]




start <- sample(1:819, size=1)
Y <- Y[,(start:(start+127))]
library(susieR)

X <- N3finemapping$X



X <-X[sample(1:nrow(X), size=nrow(Y)),1:200]

X     <-  apply(X, 2, function(x){
  x[which(is.na(x))] <- median(x, na.rm=T)
  return(x)
})




genotype <-X

idx <- which( apply( genotype,2, var ) <1e-15)
if (length(idx)>0){

  genotype <- genotype [, -idx]
}
L=2


true_pos <- sample(1:ncol(genotype), size=L)
lev_res=7
lf <-  list()
for(l in 1:L){
  lf[[l]] <- rep (0, 2^lev_res)

  idx <- sample( size =1,1:( -20+ 2^lev_res))
  lf[[l]][(idx: (idx+10))] <-1
  #functional effect for effect l
}

B <- do.call(cbind, lf)
Ynew <- t(thin_base(mat = t(Y ), designmat = genotype[, true_pos], coefmat = B))

res2 <- Pois_fSuSiE (Y=Ynew,X=genotype, L=5, max.iter=10, cor_small =TRUE, cal_obj.fsusie = TRUE)
res4 <-acc_Pois_fSuSiE2(Y= Ynew ,X=genotype, L=5 , post_processing = "HMM", cal_obj.fsusie = TRUE)

res0 <-susiF (Y=log(Ynew+1),X=genotype, L=5 , cor_small = TRUE, post_processing="none")
res01 <-HF_susiF2  (Y=Ynew,X=genotype, L=5  )


res3 <- susiF  (Y=Ynew,X=genotype, L=5, cor_small = FALSE , post_processing="none" )

Y= Ynew
X=genotype
dim(Y)
dim(X)
#plot (res4$susiF.obj$fitted_func[[1]])
#lines (res4$susiF.obj$lfsr_func[[1]])
#lines (lf[[1]], col="green")
#lines (lf[[2]], col="green")




false_disc1 <- -9
true_disc1  <- -9
false_disc2 <- -9
true_disc2  <- -9
if ( true_pos[1] %in% res4$susiF.obj$cs[[1]]){

  no_zero1 <-  which(res4$susiF.obj$lfsr_func[[1]] <0.05
  )
  if (length(no_zero1)>0){
    false_disc1 <- length( which ( no_zero1 %!in% which(lf[[1]]>0)))

    true_disc1 <- length(which ( no_zero1 %in% which(lf[[1]]>0)))
  }

}
if ( length(res4$susiF.obj$cs)>1){
  if ( true_pos[1] %in% res4$susiF.obj$cs[[2]]){

    no_zero1 <-  which(res4$susiF.obj$lfsr_func[[2]] <0.05)
    if (length(no_zero1)>0){
      false_disc1 <- length( which ( no_zero1 %!in% which(lf[[1]]>0)))

      true_disc1 <- length(which ( no_zero1 %in% which(lf[[1]]>0)))
    }

  }
}



if ( true_pos[2] %in% res4$susiF.obj$cs[[1]]){

  no_zero2 <-  which(res4$susiF.obj$lfsr_func[[1]] <0.05 )
  if (length(no_zero2)>0){
    false_disc2 <- length( which ( no_zero2 %!in% which(lf[[2]]>0)))

    true_disc2 <- length(which ( no_zero2 %in% which(lf[[2]]>0)))
  }

}
if ( length(res4$susiF.obj$cs)>1){
  if ( true_pos[2] %in% res4$susiF.obj$cs[[2]]){

    no_zero2 <-  which(res4$susiF.obj$lfsr_func[[2]] <0.05
    )
    if (length(no_zero2)>0){
      false_disc2 <- length( which ( no_zero2 %!in% which(lf[[2]]>0)))

      true_disc2 <- length(which ( no_zero2 %in% which(lf[[2]]>0)))
    }

  }
}

perf = c( false_disc1  ,
          true_disc1   ,
          false_disc2   ,
          true_disc2   )
perf

out <- list( true_pos=true_pos,
             idxt= idxt,
             perf = perf,
             susiF_pip= res0$pip,
             susiF_cs= res0$cs,

             HF_susiF_pip= res01$pip,
             HF_susiF_cs= res01$cs,


             g_susiF=res2$susiF.obj$pip,
             g_susiF_cs=res2$susiF.obj$cs,

             acc_g_susiF=res3$pip,
             acc_g_susiF_cs=res3$cs,

             acc_g_susiF2=res4$susiF.obj$pip,
             acc_g_susiF2_cs=res4$susiF.obj$cs

)


