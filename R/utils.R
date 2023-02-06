

#' @param x a vector of count of length of 2^S
#' @param indx_list a list generated susiF.alpha::gen_wavelet_indx
#' sort intensity in the same order than wavelthresh package, allow use of susiF.alpha routines
get_empirical_intensity <- function(x,indx_lst)
{

  J <- log2(length(x))
  Y_min <- rep(NA, length(x))
  Y_tot <-  rep(NA, length(x))
  for (s in 0 :(J-1)){
    nD = 2^(J - s)
    nDo2 = nD/2
    twonD = 2 * nD
    tt <-1
    for( l in 0:(2^s-1)){

    #  print(s)
      ind <- (l * nD + 1):((l + 1) * nD)
      #print(paste("ind ",ind))
      ind_len <-  length(ind)

      ind_l = ind[1:(ind_len/2)]
      #print(paste("ind_l ", ind_l))
      Y_min[indx_lst[[(s+1)]][tt]]  <-  sum(x[ind_l])
      Y_tot[indx_lst[[(s+1)]][tt]]  <-  sum(x[ind])


      tt <- tt+1



    }

  }
  sx <- sum(x)
  Y_min[length(x)] <- sx
  Y_tot[length(x)] <- sx

  return(list(Y_min= Y_min,
              Y_tot= Y_tot)
  )
}








logit <- function( p)
{
  if( p>1 & p<0)
  {
    print(" p must be between 0 and 1, introducing NA")
    out <-NA
  }
  else{
    out <- log( p/(1-p))
  }

  return(out)
}










#### From other package -----



log1pexp = function (x){
  indx <- .bincode(x, c(-Inf, -37, 18, 33.3, Inf), right = TRUE,
                   include.lowest = TRUE)
  kk <- which(indx == 1)
  if (length(kk)) {
    x[kk] <- exp(x[kk])
  }
  kk <- which(indx == 2)
  if (length(kk)) {
    x[kk] <- log1p(exp(x[kk]))
  }
  kk <- which(indx == 3)
  if (length(kk)) {
    x[kk] <- x[kk] + exp(-x[kk])
  }
  return(x)
}





















reflect_vec <- function (x)
{
  n = length(x)
  J = log2(n)
  if ((J%%1) == 0) {
    x = c(x, x[n:1])
    return(list(x = x, idx = 1:n))
  }
  else {
    n.ext = 2^ceiling(J)
    lnum = round((n.ext - n)/2)
    rnum = n.ext - n - lnum
    if (lnum == 0) {
      x.lmir = NULL
    }
    else {
      x.lmir = x[lnum:1]
    }
    if (rnum == 0) {
      x.rmir = NULL
    }
    else {
      x.rmir = x[n:(n - rnum + 1)]
    }
    x.ini = c(x.lmir, x, x.rmir)
    x.mir = x.ini[n.ext:1]
    x = c(x.ini, x.mir)
    return(list(x = x, idx = (lnum + 1):(lnum + n)))
  }
}
