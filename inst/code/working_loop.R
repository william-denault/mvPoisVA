
#data simulation
N=100
Tp <- 63
Y <-matrix( rpois(N*Tp, lambda = 2), ncol=Tp)+1#offset to avoid problem due to 0 count


### parameters/arguments
#b_pm_init start posterior mean
reflect =FALSE
verbose=TRUE
n_gh = 10 #nb points for Gauss Hermite quadrature


### data formating ------

init=TRUE
gh_points = fastGHQuad::gaussHermiteData(n_gh)
J = log2(ncol(Y)); if((J%%1) != 0) reflect=TRUE
if(reflect){
  tl <- lapply(1:nrow(Y), function(i) reflect_vec(Y[i,]))
  Y <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$x))
  idx_out <- tl[[1]]$idx #### indx of interest at the end
}



indx_lst <-  susiF.alpha::gen_wavelet_indx(log2(ncol(Y)))

tl <-  lapply(1:nrow(Y), function(i)
  get_empirical_intensity(Y[i,],
                          indx_lst = indx_lst)
)

Y_min <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$Y_min))
Y_tot <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$Y_tot))
b_pm = 0*Y
# if(missing(b_pm_init)){
#   b_pm = 0*Y
# }else{
#   b_pm = b_pm_init #posterior mean
# }



indx_lst <-  susiF.alpha::gen_wavelet_indx(log2(ncol(Y)))

tl <-  lapply(1:nrow(Y), function(i)
                              get_empirical_intensity(Y[i,],
                                                      indx_lst = indx_lst)
                )

Y_min <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$Y_min))
Y_tot <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$Y_tot))
rm(tl)

if(verbose){
  print("done transforming data")
}
# get the matrix for Binomial reg (and Pois reg for top coefficient)



sum(is.na(Y_min))

#deal with case exactly one or exactly 0

Mu_pm = logit((Y_min/Y_tot) ) #remove last column contain C coeff
Mu_pm[Mu_pm==-Inf] =  logit(0.1)
Mu_pm[Mu_pm==Inf]  =  logit(0.9)

Mu_pm[,ncol(Y_min)] <- log(Y_min[,ncol(Y_min)])
Mu_pv = 1/Y_tot
  Mu_pm[Mu_pm==-Inf] =  logit(0.1)
  Mu_pm[Mu_pm==Inf] =  logit(0.9)



  ### basic working exemple
  init_val_bin = c(c(Mu_pm[ ,-ncol(Y_min)]),log(c(Mu_pv[ ,-ncol(Y_min)])))
  init_val_pois = log(Y_min[,ncol(Y_min)])

  sigma2_bin  =1
  sigma2_pois =1
  beta_bin  <-  c(b_pm[ ,-ncol(Y_min)])
  beta_pois <-  c(b_pm[ , ncol(Y_min)])


  init_val=init_val_bin
  x=c(Y_min [ ,-ncol(Y_min)])
  nb=c(Y_tot[ ,-ncol(Y_min)])
  beta= beta_bin
  sigma2=sigma2_bin

  opt_binomial<- vga_binomial(init_val=init_val_bin,
                              x=c(Y_min [ ,-ncol(Y_min)]),
                              nb=c(Y_tot[ ,-ncol(Y_min)]),
                              beta= beta_bin,
                              sigma2=sigma2_bin,
                              gh_points=gh_points)
  opt_Poisson <- vga_pois_solver (init_val= init_val_pois,
                                  x=Y_min[,ncol(Y_min)],
                                  s= rep( 1, nrow(Y)),
                                  beta= beta_pois,
                                  sigma2=sigma2_pois,
                                  maxiter=10,
                                  tol=tol,
                                  method = 'newton')
  A_pm <- cbind(matrix(opt_binomial$m, ncol = (ncol(Y_min)-1)), opt_Poisson$m) # we are missing C column

  A_pv <- cbind(matrix(opt_binomial$v, ncol = (ncol(Y_min)-1)), opt_Poisson$v) # we are missing C column
  plot( Mu_pm[1,], A_pm[1,])
  abline(a=0,b=1)


  #need Beta and beta posterior variance
  # Update sigma2
  ##sigma2 = mean(opt_binomial$m^2+opt_binomial$v+beta_bin^2+b_pv-2*b_pm*opt_binomial$m) #
  #Posterior variance fitted value??





  tt <- reverse_intensity_transform (vec_int = c(log( sigmoid( A_pm[1,-ncol(A_pm)])), A_pm[1,ncol(A_pm)]) ,
                                     indx_lst = indx_lst,
                                     is.logprob=TRUE,
                                     is.log_int =TRUE)
  plot(tt)
  lines(tt, col="green")
  lines(Y[1,])

