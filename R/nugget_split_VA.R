#'@title Smooth over-dispersed Poisson sequence
#'@param x data vector
#'@param maxiter,tol max iteration and tolerance for stopping it.
#'@param Eb_init,sigma2_init initial values of smooth mean and nugget effect.
#'@examples
#' set.seed(12345)
#' n=2^9
#' sigma=0.5
#' mu=c(rep(0.3,n/4), rep(3, n/4), rep(10, n/4), rep(0.3, n/4))
#' x = rpois(n,exp(log(mu)+rnorm(n,sd=sigma)))
#' fit = pois_smooth_split(x,maxiter=30)
#' plot(x,col='grey80')
#' lines(exp(fit$Eb))
#' fit$sigma2
#' plot(fit$obj)
#'@details The problem is
#'\deqn{x_i\sim Poisson(\exp(\mu_i)),}
#'\deqn{\mu_i\sim N(b_i,\sigma^2),}
#'\deqn{\b_i\sim g(.).}
#'@export

pois_smooth_split = function(x,
                             s = NULL,
                             Eb_init = NULL,
                             sigma2_init = NULL,
                             est_sigma2 = TRUE,
                             maxiter = 100,
                             tol=1e-5,
                             filter.number = 1,
                             family = 'DaubExPhase',
                             verbose=FALSE,
                             printevery = 10,
                             ebnm_params=list(mode=0),
                             optim_method='L-BFGS-B'){

  n = length(x)
  if(is.null(s)){
    s = 1
  }
  if(length(s)==1){
    s = rep(s,n)
  }
  if(is.null(Eb_init)){
    Eb = log(runmed(x/s,1 + 2 * min((n-1)%/% 2, ceiling(0.1*n)))+0.01)
  }else{
    Eb = Eb_init
  }
  if(is.null(sigma2_init)){
    sigma2 = var(log(x/s+0.01)-Eb)
  }else{
    sigma2 = sigma2_init
  }

  W = (t(GenW(n,filter.number,family)))[-1,]

  mu_pm = rep(0,n)
  mu_pv = rep(1/n,n)
  obj = -Inf

  for(iter in 1:maxiter){
    # get m, s^2
    opt = optim(c(mu_pm,log(mu_pv)),
                fn = pois_mean_GG_opt_obj,
                gr = pois_mean_GG_opt_obj_gradient,
                x=x,
                #    s=s,
                beta=Eb,
                sigma2=sigma2,
                n=n,
                method = optim_method)
    mu_pm = opt$par[1:n]
    mu_pv = exp(opt$par[(n+1):(2*n)])
    qb = smash_dwt(mu_pm,sqrt(sigma2),filter.number=filter.number,family=family,ebnm_params=ebnm_params,W=W)
    Eb = qb$mu.est
    Eb2 = qb$mu.est.var + Eb^2
    # get sigma2
    if(est_sigma2){
      sigma2 = mean(mu_pm^2+mu_pv+Eb2-2*mu_pm*Eb)
    }


    # calc obj
    obj[iter+1] = pois_smooth_split_obj(x,s,mu_pm,mu_pv,Eb,Eb2,sigma2,qb$dKL)

    if(verbose){
      if(iter%%printevery==0){
        print(paste("Done iter",iter,"obj =",obj[iter+1]))
      }

    }

    if((obj[iter+1]-obj[iter])<tol){
      break
    }

  }
  return(list(Emean = exp(mu_pm+mu_pv/2),
              Vmean = exp(mu_pv-1)*exp(2*mu_pm+mu_pv),
              Emu=mu_pm,
              Vmu=mu_pv,
              Eb=Eb,
              Eb2=Eb2,
              sigma2=sigma2,
              obj=obj,
              H = qb$dKL + sum(log(2*pi*mu_pv)/2-log(2*pi*sigma2)/2-(mu_pm^2+mu_pv-2*mu_pm*Eb+Eb2)/2/sigma2)))
}

pois_smooth_split_obj = function(x,s,m,s2,Eb,Eb2,sigma2,KLb){
  return(sum(x*m-s*exp(m+s2/2)+log(s2)/2-log(sigma2)/2-(m^2+s2-2*m*Eb+Eb2)/2/sigma2)+KLb)
}



#'@title Empirical Bayes wavelet smoothing via DWT
#'@description Smooth homogeneous Gaussian data.
#'@param x data
#'@param sigma known standard error
#'@param filter.number,family wavelet family and filter number as in wavethresh package
#'@param ebnm_params a list of `ebnm` parameters
#'@param W the dwt matrix for calc posterior variance. Remove the first row which is all 1/sqrt(n)
#'@return a list of
#'  \item{mu.est:}{posterior mean}
#'  \item{mu.est.var:}{posterior variance}
#'  \item{loglik:}{log likelihood}
#'  \item{dKL:}{KL divergence between g(the prior) and q(the posterior)}
#'@import wavethresh
#'@import ebnm
#'@export
smash_dwt = function(x,sigma,filter.number=1,
                     family="DaubExPhase",
                     ebnm_params=list(mode=0),W=NULL){

  n = length(x)
  J = log(n,2)
  if(ceiling(J)!=floor(J)){
    stop('Length of x must be power of 2')
  }
  if(is.null(ebnm_params)){
    ebnm_params = ebnm_params_default()
  }else{
    temp = ebnm_params_default()
    for(i in 1:length(ebnm_params)){
      temp[[names(ebnm_params)[i]]] = ebnm_params[[i]]
    }
    ebnm_params = temp
  }
  tsum = sum(x)/sqrt(n)
  x.w = wd(x, filter.number = filter.number,
           family = family, type = "wavelet")

  data.var = sigma^2
  if(length(data.var==1)){
    data.var = rep(data.var,n)
  }

  if(is.null(W)){
    W = (t(GenW(n,filter.number,family)))[-1,]
  }

  if(length(sigma)==1){
    x.w.v = rep(sigma^2,n-1)
    tsum.var = sigma^2
  }else{
    x.w.v =  data.var
    tsum.var = x.w.v[1]
    x.w.v = x.w.v[-1]
  }

  dKL = 0
  loglik.scale = c()
  x.w.v.s = rep(0, 2^J-1)
  for (j in 0:(J - 1)) {
    x.pm = rep(0, 2^j)
    #index = (((J - 1) - j) * n + 1):((J - j) * n)
    index = (n-2^(j+1)+1):(n-2^j)
    x.w.j = accessD(x.w, j)
    x.w.v.j = x.w.v[index]
    ind.nnull = (x.w.v.j != 0)

    a = ebnm(x.w.j[ind.nnull],sqrt(x.w.v.j[ind.nnull]),
             mode=ebnm_params$mode,
             prior_family=ebnm_params$prior_family,
             scale = ebnm_params$scale,
             g_init = ebnm_params$g_init,
             fix_g = ebnm_params$fix_g,
             output = ebnm_params$output,
             optmethod = ebnm_params$optmethod)

    dKL = dKL + a$log_likelihood - Eloglik(x.w.j[ind.nnull], sqrt(x.w.v.j[ind.nnull]),a$posterior$mean, a$posterior$mean^2+a$posterior$sd^2)
    x.pm[ind.nnull] = a$posterior$mean
    x.pm[!ind.nnull] = 0
    x.w = putD(x.w, j, x.pm)
    loglik.scale[j + 1] = a$log_likelihood
    x.w.v.s[index[ind.nnull]] = a$posterior$sd^2
    x.w.v.s[index[!ind.nnull]] = 0
  }
  mu.est = wr(x.w)
  loglik = sum(loglik.scale)
  #x.w.v.s = c(tsum.var,x.w.v.s)
  mu.est.var = colSums(W^2*x.w.v.s)
  return(list(mu.est=mu.est,mu.est.var=mu.est.var,loglik = loglik,dKL = dKL))
}



Eloglik = function(x, s, Et, Et2) {
  # Deal with infinite SEs:
  idx = is.finite(s)
  x = x[idx]
  s = s[idx]
  Et = Et[idx]
  Et2 = Et2[idx]
  return(-0.5 * sum(log(2*pi*s^2) + (1/s^2) * (Et2 - 2*x*Et + x^2)))
}




## previously, I used a for loop to solve for posterior mean, variance for each observation.
## This function is a vectorised version of the for loop + pois_mean_GG1

## This is not faster than the 1-by-1 version. Very likely due to the numerical calc of hessian matrix.
pois_mean_GG_opt = function(x,
                            beta,
                            sigma2,
                            optim_method = 'BFGS',
                            m_init  = NULL,
                            s2_init = NULL){
  n = length(x)
  # init m, s2
  if(is.null(m_init)){
    m = rep(0,n)
  }else{
    m = m_init
  }

  if(is.null(s2_init)){
    s2 = rep(1,n)
  }else{
    s2 = s2_init
  }


  opt = optim(c(m,log(s2)),
              fn = pois_mean_GG_opt_obj,
              gr = pois_mean_GG_opt_obj_gradient,
              x=x,
              beta=beta,
              sigma2=sigma2,
              n=n,
              method = optim_method)

  return(list(m=opt$par[1:n],s2=exp(opt$par[(n+1):(2*n)]),obj=-opt$value))

  # opt = nlm(f_obj_nlm,c(m,log(s2)),x=x,
  #           beta=beta,
  #           sigma2=sigma2,
  #           n=n)
  # return(list(m=opt$estimate[1:n],s2=exp(opt$estimate[(n+1):(2*n)]),obj=-opt$minimum))

}

#'calculate objective function
pois_mean_GG_opt_obj = function(theta,x,beta,sigma2,n){
  m = theta[1:n]
  v = theta[(n+1):(2*n)]
  return(-sum(x*m-exp(m+exp(v)/2)-(m^2+exp(v)-2*m*beta)/2/sigma2+log(exp(v))/2))
}

#'calculate gradient vector
pois_mean_GG_opt_obj_gradient = function(theta,x,beta,sigma2,n){
  m = theta[1:n]
  v = theta[(n+1):(2*n)]
  g1 = -(x-exp(m+exp(v)/2)-m/sigma2+beta/sigma2)
  g2 = -(-exp(v)/2*exp(m+exp(v)/2) - exp(v)/2/sigma2 + 1/2)
  return(c(g1,g2))
}

f_obj_nlm = function(theta,x,beta,sigma2,n){
  m = theta[1:n]
  v = theta[(n+1):(2*n)]
  out = -sum(x*m-exp(m+exp(v)/2)-(m^2+exp(v)-2*m*beta)/2/sigma2+log(exp(v))/2)
  g1 = -(x-exp(m+exp(v)/2)-m/sigma2+beta/sigma2)
  g2 = -(-exp(v)/2*exp(m+exp(v)/2) - exp(v)/2/sigma2 + 1/2)
  attr(out,'gradient') = c(g1,g2)
  return(out)
}
output_default <- ebnm:::ebnm_output_default

