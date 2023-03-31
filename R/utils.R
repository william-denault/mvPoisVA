

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


compute_inverse_effect <- function(lp,lq, l_lambda){




  ## inverse function of get_empirical_intensity


    J=log2(length(lp)+1)

      log_lambda_tot <-l_lambda




    out <- rep( log_lambda_tot , 2^J)

    for(s in (J ):1){

      nD = 2^(J-s+1)
      nDo2 = nD/2
      tt <-1
      for(l in 0:(2^(s-1)-1)){
        ind = (l*nD+1):((l+1)*nD) # all "sub index for coef s,l (here s=D)

        print(ind)
        ind_l <-  ind[1:nDo2] #all "sub index in the left for coef s,l (here s=D)
        ind_r <-  ind[(nDo2+1):nD] # all "sub index in the right for coef s,l (here s=D)
        out[ind_l] <- out[ind_l]+ lp[indx_lst[[(s )]][tt]]
        out[ind_r] <- out[ind_r]+ lq[indx_lst[[(s )]][tt]]
        tt <- tt+1
      }
    }
    return( )


}


## inverse function of get_empirical_intensity

reverse_intensity_transform =function(vec_int, indx_lst,
                                      is.logprob=TRUE,
                                      is.prob=FALSE,
                                      is.log_int=TRUE){

  J=log2(length(vec_int))

  if( is.prob){
    lp <- log(vec_int[- length(vec_int)])
    lq = log(1-pmin(exp(lp),1-1e-10))
  }
  if(is.logprob){

    lp <-  (vec_int[- length(vec_int)])
    lq =   log(1-pmin( exp(lp),1-1e-10))
  }
  if( is.log_int){
    log_lambda_tot <- vec_int[  length(vec_int)]
  }else{
    print( "assuming input is total intensity ")
    log_lambda_tot <- log(vec_int[  length(vec_int)])
  }



  out <- rep( log_lambda_tot , 2^J)

  for(s in (J ):1){

    nD = 2^(J-s+1)
    nDo2 = nD/2
    tt <-1
    for(l in 0:(2^(s-1)-1)){
      ind = (l*nD+1):((l+1)*nD) # all "sub index for coef s,l (here s=D)

      print(ind)
      ind_l <-  ind[1:nDo2] #all "sub index in the left for coef s,l (here s=D)
      ind_r <-  ind[(nDo2+1):nD] # all "sub index in the right for coef s,l (here s=D)
      out[ind_l] <- out[ind_l]+ lp[indx_lst[[(s )]][tt]]
      out[ind_r] <- out[ind_r]+ lq[indx_lst[[(s )]][tt]]
      tt <- tt+1
    }
  }
  return(exp(out))

}

#alpha, lambda_tot




#' @export
logit = function(x){
  log(x/(1-x))
}
#' @export
sigmoid = function(x){
  1/(1+exp(-x))
}

get_post_log_int <- function(Mu_pm,
                             Mu_pv,
                             Y_min,
                             Y_tot,
                             sigma2_bin,
                             sigma2_pois,
                             b_pm,
                             gh_points,
                             tol=1e-5){





  Mu_pm[Mu_pm==-Inf] =  logit(0.1)
  Mu_pm[Mu_pm==Inf] =  logit(0.9)
  ### basic working exemple
  init_val_bin = c(c(Mu_pm[ ,-ncol(Y_min)]),log(c(Mu_pv[ ,-ncol(Y_min)])))
  init_val_pois =  Mu_pm[,ncol(Y_min)]


  beta_bin  <-  c(b_pm[ ,-ncol(Y_min)])
  beta_pois <-  c(b_pm[ , ncol(Y_min)])


  opt_binomial <- vga_binomial(init_val  = init_val_bin,
                               x         = c(Y_min [ ,-ncol(Y_min)]),
                               nb        = c(Y_tot[ ,-ncol(Y_min)]),
                               beta      = beta_bin,
                               sigma2    = sigma2_bin,
                               gh_points = gh_points)


  opt_Poisson  <- vga_pois_solver(init_val = init_val_pois ,
                                  x        = Y_min[,ncol(Y_min)],
                                  s        = rep( 1, nrow(Y)),
                                  beta     = beta_pois,
                                  sigma2   = sigma2_pois,
                                  maxiter  = 10,
                                  tol      = tol,
                                  method   = 'newton')


  A_pm <- cbind(matrix(opt_binomial$m, ncol = (ncol(Y_min)-1)), opt_Poisson$m) # we are missing C column

  A_pv <- cbind(matrix(opt_binomial$v, ncol = (ncol(Y_min)-1)), opt_Poisson$v) # we are missing C column

  #need Beta and beta posterior variance
  # Update sigma2
  ##sigma2 = mean(opt_binomial$m^2+opt_binomial$v+beta_bin^2+b_pv-2*b_pm*opt_binomial$m) #
  #Posterior variance fitted value??

  return(  list(A_pm = A_pm,
                A_pv = A_pv)
        )
}



# compute the individual fitted inhomogeneous Poisson process
#In put is a matrix of posterior mean of logit of fitted probability
# for each binomial regression, the last column corresponds to
#the log intensity of the overall process
get_ind_fitted_Poisproc <- function(post_mat,indx_lst ){

  fitted_log_prob <- log(sigmoid(post_mat [ ,-ncol(post_mat )]))
  fitted_log_int  <- post_mat [ ,  ncol(post_mat )]
  fitted_Pois     <- cbind(fitted_log_prob,fitted_log_int)

  fitted_Pois <- lapply(1:nrow( fitted_Pois),
                        function( i)
                          reverse_intensity_transform(vec_int  = fitted_Pois[i,],
                                                      indx_lst = indx_lst,
                                                      is.logprob=TRUE,
                                                      is.log_int =TRUE) )
  fitted_Pois <- do.call(rbind,fitted_Pois)
  return( fitted_Pois)
}

get_fitted_Poisproc <-function( EBmvFR.obj,c_mean,indx_lst){

  fitted_log_prob   <- log(sigmoid( EBmvFR.obj$fitted_wc[[1]] [ ,-ncol(EBmvFR.obj$fitted_wc[[1]]  )]))
  fitted_log_int    <-EBmvFR.obj$fitted_wc[[1]] [ ,  ncol(EBmvFR.obj$fitted_wc[[1]]  )]
  fitted_effect     <- cbind(fitted_log_prob,fitted_log_int)
  fitted_effect     <- fitted_effect+ matrix(c_mean,
                                             byrow = TRUE,
                                             nrow=nrow(fitted_effect),
                                             ncol=ncol(fitted_effect)
  )

  fitted_effect <- lapply(1:nrow( fitted_effect),
                          function( i)
                            reverse_intensity_transform(vec_int  = fitted_effect[i,],
                                                        indx_lst = indx_lst,
                                                        is.logprob=TRUE,
                                                        is.log_int =TRUE) )
  fitted_effect <- do.call(rbind,fitted_effect)
  return( fitted_Pois)
}



extend_vec  <- function (x)
{
  n = length(x)
  J = log2(n)
  if ((J%%1) == 0) {
    x = c(x, x[n:1])
    return(list(x = x, idx = 1:n))
  }  else {
    n.ext = 2^ceiling(J)
    lnum = round((n.ext - n)/2)
    rnum = n.ext - n - lnum
    if (lnum == 0) {
      x.lmir = NULL
    }    else {
      x.lmir = x[lnum:1]
    }
    if (rnum == 0) {
      x.rmir = NULL
    }    else {
      x.rmir = x[n:(n - rnum + 1)]
    }


    x =  c(x.lmir, x, x.rmir)
    return(list(x = x, idx = (lnum + 1):(lnum + n)))
  }
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

























#' Interleave two vectors.
#' @keywords internal
interleave=function(x,y){
  return(as.vector(rbind(x,y)))
}
