#'@title Solve Gaussian approximation to Poisson mean problem
#'@description Gaussian prior, Gaussian posterior in Poisson mean problem.
#'@param x data vector
#'@param s scaling vector
#'@param w prior weights
#'@param prior_mean prior mean
#'@param prior_var prior variance
#'@param optim_method optimization method in `optim` function
#'@param maxiter max number of iterations
#'@param tol tolerance for stopping the updates
#'@return a list of
#'  \item{posteriorMean:}{posterior mean}
#'  \item{posteriorVar:}{posterior variance}
#'  \item{obj_value:}{objective function values}
#'  \item{prior_mean:}{prior mean}
#'  \item{prior_var:}{prior variance}
#'  @example
#'  n = 10000
#'  mu = rnorm(n)
#'  x = rpois(n,exp(mu))
#'  pois_mean_GG(x)
#'@details The problem is
#'\deqn{x_i\sim Poisson(\exp(\mu_i)),}
#'\deqn{\mu_i\sim N(\beta,\sigma^2).}
#'@export
pois_mean_GP = function(x,
                        s = NULL,
                        prior_mean = NULL,
                        prior_var=NULL,
                        optim_method = 'L-BFGS-B',
                        maxiter = 1000,
                        tol = 1e-5){

  # init the posterior mean and variance?
  n = length(x)
  m = log(x+0.1)
  v = rep(1/ sqrt(n),n)
  if(is.null(s)){
    s = 1
  }
  if(length(s)==1){
    s = rep(s,n)
  }
  #
  if(is.null(prior_mean) | is.null(prior_var)){

    if(is.null(prior_mean)){
      est_beta = TRUE
    }else{
      est_beta = FALSE
      beta = prior_mean
    }
    if(is.null(prior_var)){
      est_sigma2=TRUE
    }else{
      est_sigma2 = FALSE
      sigma2=prior_var
    }

    obj = rep(0,maxiter+1)
    obj[1] = -Inf
    for(iter in 1:maxiter){
      if(est_beta){
        beta = mean(m)
      }
      if(est_sigma2){
        sigma2 = mean(m^2+v-2*m*beta+beta^2)
      }
      # for(i in 1:n){
      #   temp = pois_mean_GG1(x[i],s[i],beta,sigma2,optim_method,m[i],v[i])
      #   m[i] = temp$m
      #   v[i] = temp$v
      # }
      opt = optim(c(m,log(v)),
                  fn = pois_mean_GP_opt_obj,
                  gr = pois_mean_GP_opt_obj_gradient,
                  x=x,
                  s=s,
                  beta=beta,
                  sigma2=sigma2,
                  n=n,
                  method = optim_method)
      m = opt$par[1:n]
      v = exp(opt$par[(n+1):(2*n)])
      obj[iter+1] = pois_mean_GG_obj(x,s,beta,sigma2,m,v)
      if((obj[iter+1] - obj[iter])<tol){
        obj = obj[1:(iter+1)]
        break
      }
    }

  }else{
    beta = prior_mean
    sigma2 = prior_var
    # for(i in 1:n){
    #   temp = pois_mean_GG1(x[i],s[i],prior_mean,prior_var,optim_method,m[i],v[i])
    #   m[i] = temp$m
    #   v[i] = temp$v
    # }
    opt = optim(c(m,log(v)),
                fn = pois_mean_GP_opt_obj,
                gr = pois_mean_GP_opt_obj_gradient,
                x=x,
                s=s,
                beta=beta,
                sigma2=sigma2,
                n=n,
                method = optim_method)
    m = opt$par[1:n]
    v = exp(opt$par[(n+1):(2*n)])
    obj = pois_mean_GG_obj(x,s,prior_mean,prior_var,m,v)

  }

  return(list(posterior = list(posteriorMean_latent = m,
                               posteriorVar_latent = v,
                               posteriorMean_mean = exp(m + v/2)),
              fitted_g = list(mean = beta, var=sigma2),
              obj_value=obj))

  #return(list(posteriorMean=m,priorMean=beta,priorVar=sigma2,posteriorVar=v,obj_value=obj))

}
#'calculate objective function
pois_mean_GP_opt_obj = function(theta,x,s,beta,sigma2,n){
  m = theta[1:n]
  v = theta[(n+1):(2*n)]
  return(-sum(x*m-s*exp(m+exp(v)/2)-(m^2+exp(v)-2*m*beta)/2/sigma2+v/2))
}
#'calculate gradient vector
pois_mean_GP_opt_obj_gradient = function(theta,x,s,beta,sigma2,n){
  m = theta[1:n]
  v = theta[(n+1):(2*n)]
  g1 = -(x-s*exp(m+exp(v)/2)-m/sigma2+beta/sigma2)
  g2 = -(-exp(v)/2*s*exp(m+exp(v)/2) - exp(v)/2/sigma2 + 1/2)
  return(c(g1,g2))
}


pois_mean_GG_obj = function(x,s,beta,sigma2,m,v){
  return(sum(x*m-s*exp(m+v/2)-log(sigma2)/2-(m^2+v-2*m*beta+beta^2)/2/sigma2+log(v)/2))
}
