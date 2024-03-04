rm(list=ls())
library(ashr)
library(mvPoisVA)
library(fsusieR)
library(susieR)
#o=13
#o=46 pb with effect =0.4

#79
'%!in%' <- function(x,y)!('%in%'(x,y))

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

load("C:/Document/Serieux/Travail/Package/susiF.alpha/inst/data_RCC/IRF2_locus_3.RData")
Y <- data$Y$Y_f[[1]]




start <- 495
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
  lf[[l]][(idx: (idx+10))] <-10
  #functional effect for effect l
}

B <- do.call(cbind, lf)
Ynew <- t(thin_base(mat = t(Y ), designmat = genotype[, true_pos], coefmat = B))
plot(Ynew,Y )

plot(lf[[1]])
par (mfrow=c(1,2))
plot (Ynew[1,], type="l",
      xlim = c(min(which(lf[[1]]>0))-5,
               max(which(lf[[1]]>0))+5),
      ylim = c(-.5,30))
for ( i in 1:nrow(Ynew)){
  lines (Ynew[i,], col= ceiling(X[i,true_pos[1]]+2))
}
lines(lf[[1]], lwd=2)


plot (Y [1,], type="l",
      xlim = c(min(which(lf[[1]]>0))-5,
               max(which(lf[[1]]>0))+5),
      ylim = c(-.5,30))
for ( i in 1:nrow(Ynew)){
  lines (Y [i,], col= ceiling(X[i,true_pos[1]]+2))
}
lines(lf[[1]], lwd=2)





par (mfrow=c(1,1))


true_pos

res1 <- Pois_fSuSiE(Y= Ynew ,X=genotype, L=5 , cor_small = TRUE )
res1$cs
plot (res1$susiF.obj$fitted_func[[1]])

lines(lf[[2]], lwd=2)
lines(lf[[1]], lwd=2)


plot (res1$fitted_func[[2]])
lines(res1$lfsr_func[[2]])
lines(lf[[2]], lwd=2)
lines(lf[[1]], lwd=2)

res4 <-acc_Pois_fSuSiE2(Y= Ynew ,
                        X=genotype,
                        L=5 ,
                        post_processing = "HMM",
                        cal_obj.fsusie = TRUE)
res4$susiF.obj$cs

plot (res4$susiF.obj$fitted_func[[1]])
lines(lf[[2]], lwd=2)

lines(lf[[1]], lwd=2)
lines(res4$susiF.obj$lfsr_func[[1]], lwd=2, col="red")
abline(h=0.05)


plot (res4$susiF.obj$fitted_func[[2]])
lines(lf[[2]], lwd=2)

lines(lf[[1]], lwd=2)
lines(res4$susiF.obj$lfsr_func[[2]], lwd=2, col="red")
abline(h=0.05)

Y=log(Ynew+1)
Y=colScale(Y, scale=FALSE)
X=fsusieR::colScale(genotype)

susiF.obj <-res1



plot ( res1$fitted_func[[1]])

lines(lf[[1]], lwd=2)
lines(lf[[2]], lwd=2)
lines(res1$lfsr_func[[1]], lwd=2, col="red")
abline(h=0.05)

