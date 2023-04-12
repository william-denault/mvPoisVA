#' @export
binomial_mean_splitting = function(x,nb=1,
                                   tol=1e-5,
                                   maxiter=1e3,
                                   ebnm_params=NULL,
                                   optim_method ='L-BFGS-B',
                                   n_gh = 20,
                                   printevery = 100,
                                   b_pm_init = NULL,
                                   sigma2_init = NULL,
                                   convergence_criteria = 'obj'){
  n = length(x)
  if(length(nb)==1){
    nb = rep(nb,n)
  }
  obj = rep(0,maxiter+1)
  obj[1] = -Inf
  if(is.null(ebnm_params)){
    ebnm_params = ebnm_params_default()
  }else{
    temp = ebnm_params_default()
    for(i in 1:length(ebnm_params)){
      temp[[names(ebnm_params)[i]]] = ebnm_params[[i]]
    }
    ebnm_params = temp
  }

  # const in objective function
  const = sum(lfactorial(nb)-lfactorial(x)-lfactorial(nb-x))

  gh_points = gaussHermiteData(n_gh)

  if(is.null(b_pm_init)){
    b_pm = rep(0,n)
  }else{
    b_pm = b_pm_init
  }

  #b_pv = rep(1/n,n)
  mu_pm = logit(x/nb)
  mu_pm[mu_pm==-Inf] = logit(0.1)
  mu_pm[mu_pm==Inf] = logit(0.9)
  mu_pv = rep(1/n,n)
  if(is.null(sigma2_init)){
    sigma2 = var(mu_pm-b_pm)
  }else{
    sigma2 = sigma2_init
  }
  t_start = Sys.time()
  for (iter in 1:maxiter) {

    opt = vga_binomial(c(mu_pm,log(mu_pv)),x,nb,b_pm,sigma2,gh_points=gh_points)
    mu_pm = opt$m
    mu_pv = opt$v
    # EBNM
    res = ebnm(mu_pm,sqrt(sigma2),
               mode=ebnm_params$mode,
               prior_family=ebnm_params$prior_family,
               scale = ebnm_params$scale,
               g_init = ebnm_params$g_init,
               fix_g = ebnm_params$fix_g,
               output = ebnm_params$output,
               optmethod = ebnm_params$optmethod)
    b_pm = res$posterior$mean
    b_pv = res$posterior$sd^2
    H = res$log_likelihood + n*(log(2*pi*sigma2)/2)+sum((mu_pm^2-2*mu_pm*b_pm+b_pm^2+b_pv)/sigma2/2)

    # Update sigma2
    sigma2 = mean(mu_pm^2+mu_pv+b_pm^2+b_pv-2*b_pm*mu_pm)

    # ELBO
    obj[iter+1] = sum(x*mu_pm-nb*Elog1pexp(mu_pm,mu_pv,gh_points)) +const - n/2*log(2*pi*sigma2)
    - sum(mu_pm^2 + mu_pv + b_pm^2 + b_pv - 2*mu_pm*b_pm)/2/sigma2 + H + sum(log(2*pi*mu_pv))/2 - n/2

    if(iter%%printevery==0){
      print(paste('At iter', iter, 'elbo=',round(obj[iter+1],3),'sigma2=',round(sigma2,3)))
    }
    if(convergence_criteria=='obj'){
      if((obj[iter+1]-obj[iter])<tol){
        obj = obj[1:(iter+1)]
        if((obj[iter+1]-obj[iter])<0){
          warning('An iteration decreases ELBO. This is likely due to numerical issues.')
        }
        break
      }
    }
    if(convergence_criteria=='objabs'){
      if(abs(obj[iter+1]-obj[iter])<tol){
        obj = obj[1:(iter+1)]
        break
      }
    }


  }
  t_end = Sys.time()
  return(list(posterior = list(mean_logit = mu_pm,
                               mean_b = b_pm,
                               mean = Esigmoid(mu_pm,mu_pv,gh_points)),
              fitted_g = list(sigma2=sigma2,g_b = res$fitted_g),
              elbo=obj[length(obj)],
              obj_trace = obj,
              fit = res,
              run_time = difftime(t_end,t_start,units='secs')))

}



#' @export

#data simulation
#N=100
#Tp <- 63
#Y <-matrix( rpois(N*Tp, lambda = 20), ncol=Tp)+1#offset to avoid problem due to 0 count


### parameters/arguments
#b_pm_init start posterior mean
#reflect =FALSE
#verbose=TRUE
#n_gh = 10 #nb points for Gauss Hermite quadrature
#

### data formating ------


#b_pm = 0*Y
# if(missing(b_pm_init)){
#   b_pm = 0*Y
# }else{
#   b_pm = b_pm_init #posterior mean
# }


#gh_points = fastGHQuad::gaussHermiteData(n_gh)
#J = log2(ncol(Y)); if((J%%1) != 0) reflect=TRUE
#if(reflect){
#  #  tl <- lapply(1:nrow(Y), function(i) reflect_vec(Y[i,]))
#  Y <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$x))
#  idx_out <- tl[[1]]$idx #### indx of interest at the end
#}
#
#indx_lst <-  susiF.alpha::gen_wavelet_indx(log2(ncol(Y)))
#
#tl <-  lapply(1:nrow(Y), function(i)
# get_empirical_intensity(Y[i,],
#                          indx_lst = indx_lst)
#)

#Y_min <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$Y_min))
#Y_tot <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$Y_tot))
#rm(tl)

#if(verbose){
#  print("done transforming data")
#}
# get the matrix for Binomial reg (and Pois reg for top coefficient)





#deal with case exactly one or exactly 0


#Mu_pm = logit((Y_min/Y_tot) ) #remove last column contain C coeff#

#Mu_pm[Mu_pm==-Inf] =  logit(0.1)
#Mu_pm[Mu_pm==Inf] =  logit(0.9)
#Mu_pv = 1/Y_tot




#### basic working exemple

#x <- Y_min [1,-ncol(Y_min)]
#nb <- Y_tot[1,-ncol(Y_min)]
#mu_pm <- Mu_pm[1,]#init value
#mu_pv <- rep(1/length(x),length(x))[-ncol(Y_min)]#init value
#b_pm <- 0*x
#sigma2 =2
# here everything is a vector (not sigma^2) mu_pm and mu_
#
#check line 52 VA_binomial
#opt<- vga_binomial(c(mu_pm,log(mu_pv)),x,nb,b_pm,sigma2,gh_points=gh_points)
#plot(opt$m,mu_pm)
#abline(a=0,b=1)
#mu_pm = opt$m
#mu_pv = opt$v


vga_binomial = function(init_val,x,nb,beta,sigma2,method='lbfgs',gh_points){
  if(method!="lbfgs"){
    stop('only lbfgs is implemented')
  }
  n = length(x)
  opt = try(optim(init_val,
                  fn = vga_binom_obj,
                  gr = vga_binom_obj_grad,
                  x=x,
                  nb=nb,
                  beta=beta,
                  sigma2=sigma2,
                  n=n,
                  gh_points=gh_points,
                  method = 'L-BFGS-B'),silent = TRUE)

  if(class(opt)=='try-error'){
    warning(paste(opt[1],'; returning the initialization values'))
    return(list(m=init_val[1:n],v=exp(init_val[(n+1):(2*n)])))
  }else{
    return(list(m=opt$par[1:n],v=exp(opt$par[(n+1):(2*n)]),opt=opt))
  }


}

#' @export
vga_binom_obj = function(theta,x,nb,beta,sigma2,n,gh_points){
  m = theta[1:n]
  lv = theta[(n+1):(2*n)]
  val = - sum(x*m - nb*Elog1pexp(m,exp(lv),gh_points) - (m^2+exp(lv)-2*m*beta)/2/sigma2 + lv/2)
  return(val)
}
#' @export
vga_binom_obj_grad = function(theta,x,nb,beta,sigma2,n,gh_points){
  m = theta[1:n]
  lv = theta[(n+1):(2*n)]
  dm = - (x - nb*Elog1pexp_dm(m,exp(lv),gh_points) - (m-beta)/sigma2)
  dlv = - (- nb*Elog1pexp_dlv(m,lv,gh_points) - exp(lv)/2/sigma2 + 1/2)
  return(c(dm,dlv))
}
#' @export
Elog1pexp = function(m,v,gh_points){
  mat = outer(sqrt(2*v),gh_points$x,FUN='*') + m
  return(c(log1pexp(mat)%*%gh_points$w)/sqrt(pi))
}
#' @export
Elog1pexp_dm = function(m,v,gh_points){
  mat = outer(sqrt(2*v),gh_points$x,FUN='*') + m
  return(c(sigmoid(mat)%*%gh_points$w)/sqrt(pi))
}
#' @export
Elog1pexp_dlv = function(m,lv,gh_points){
  temp = outer(sqrt(2*exp(lv)),gh_points$x,FUN='*')
  mat = temp + m
  mat2 = temp/2
  return(c((mat2*sigmoid(mat))%*%gh_points$w)/sqrt(pi))
}


