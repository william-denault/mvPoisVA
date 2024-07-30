library(mvPoisVA)
library(susiF.alpha)
library(susieR)
library(ebnm)


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

'%!in%' <- function(x,y)!('%in%'(x,y))
data(N3finemapping)
X <- N3finemapping$X
mysd=2
L=2
N =50
genotype <-X[1:N,1:100]

idx <- which( apply( genotype,2, var ) <1e-15)
genotype <- genotype [, -idx]
lev_res =8
count.data  <- list()

Y <- matrix ( rpois(N* 2^(lev_res), lambda = 10), ncol = 2^(lev_res))
image (Y)
X <- N3finemapping$X
true_pos <- sample(1:ncol(genotype), L)

lf <-  list()
for(l in 1:L){
  lf[[l]] <- rep (0, 2^lev_res)

  idx <- sample( size =1,1:( -20+ 2^lev_res))
  lf[[l]][(idx: (idx+10))] <-1
  #functional effect for effect l
}


B <- do.call(cbind, lf)
Ynew <- t(thin_base(mat = t(Y), designmat = genotype[, true_pos], coefmat = B))
plot (Ynew,Y)

res01 <-acc_Pois_fSuSiE2(Y= Ynew ,X=genotype, L=3 , post_processing = "HMM")
true_pos
plot (res01$susiF.obj$fitted_func[[1]])

lines(res01$susiF.obj$lfsr_func[[1]], col="red")
lines(lf[[1]])
lines(lf[[2]])

plot (res01$susiF.obj$fitted_func[[2]])

lines(res01$susiF.obj$lfsr_func[[2]], col="red")
lines(lf[[2]])



false_disc1 <- -9
true_disc1  <- -9
false_disc2 <- -9
true_disc2  <- -9
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

    no_zero1 <-  which(res01$susiF.obj$lfsr_func[[2]] <0.05)
    if (length(no_zero1)>0){
      false_disc1 <- length( which ( no_zero1 %!in% which(lf[[1]]>0)))

      true_disc1 <- length(which ( no_zero1 %in% which(lf[[1]]>0)))
    }

  }
}



if ( true_pos[2] %in% res01$susiF.obj$cs[[1]]){

  no_zero2 <-  which(res01$susiF.obj$lfsr_func[[1]] <0.05 )
  if (length(no_zero2)>0){
    false_disc2 <- length( which ( no_zero2 %!in% which(lf[[2]]>0)))

    true_disc2 <- length(which ( no_zero2 %in% which(lf[[2]]>0)))
  }

}
if ( length(res01$susiF.obj$cs)>1){
  if ( true_pos[2] %in% res01$susiF.obj$cs[[2]]){

    no_zero2 <-  which(res01$susiF.obj$lfsr_func[[2]] <0.05
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

