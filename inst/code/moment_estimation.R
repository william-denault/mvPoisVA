#' Compute the third and fourth moments of a normal or a mixture of normals distribution.
#' @param mu mean
#' @param sigma standard deviation
#' @param moment: if moment=3 compute the third moment, if moment=4 compute the forth moment
#' @keywords internal
tfmoment=function(mu,sigma,moment,pi=NULL){
  if(moment==3){
    moment=mu^3+3*mu*sigma^2
  }else if(moment==4){
    moment=mu^4+6*mu^2*sigma^2+3*sigma^4
  }
  #if pi is not null compute the third and fourth moments of a mixture of normals
  if(!is.null(pi))
    moment=colSums(pi*moment)
  return(moment)
}



fwrapper <- function(mean, var, Mean, Sq, wSq, beta.tm, beta.fm, g){
  ffdash = ffdash.moments(mean,var)
  ffdd = ffdd.moments(mean,var)
  ffdashdd = ffdashdd.moments(mean,var)
  lpratio.mean = ffdash$mean*Mean + is.factor(g)*wSq*0.5*ffdd$mean*Sq
  lpratio.var = pmax(ffdash$meansq*Sq + is.factor(g)*wSq*(0.25*beta.fm*ffdd$meansq*wSq + beta.tm*ffdashdd$mean) - lpratio.mean^2,0)

  return(list(mean=lpratio.mean, var=lpratio.var))
}






##Multiseq wrap up
 #baseline.mean = reverse.pwave(res.rate$lp.mean, matrix(res$lp.mean, J, n, byrow=TRUE))
#baseline.var = reverse.pwave(res.rate$lp.var, matrix(res$lp.var, J, n, byrow=TRUE), matrix(res$lq.var, J, n, byrow=TRUE))
effect.mean = reverse.pwave(res.rate$lpratio.mean, matrix(res$lpratio.mean, J, n, byrow=TRUE), matrix(res$lqratio.mean, J, n, byrow=TRUE))
effect.var = reverse.pwave(res.rate$lpratio.var, matrix(res$lpratio.var, J, n, byrow=TRUE), matrix(res$lqratio.var, J, n, byrow=TRUE))









ff_exp = function(x){
  return(exp(-log(1+exp(-x))))
}

ffdash = function(x){
  return(1/(1+exp(x)))
}

ffdd = function(x){
  ex=exp(x)
  return(-ex/(1+ex)^2)
}

ffdashdd = function(x){
  ex=exp(x)
  return(-ex/(1+ex)^3)
}

#first and second  derivatives of ff
ff.deriv_exp = function(x){
  ex = exp(x)
  return(list(first = ex/(1+2*ex+ex^2), second= (ex-ex^3)/(1+4*ex+6*ex^2+4*ex^3+ex^4)))
}

#first and second  derivatives of ffdash
ffdash.deriv = function(x){
  ex=exp(x)
  return(list(first= -ex/(1+ex)^2, second= ex*(ex-1)*(1+ex)^(-3)))
}

ffdd.deriv = function(x){
  ex=exp(x)
  return(list(first= ex*(ex-1)*(1+ex)^(-3), second= ex*(4*ex-ex^2-1)*(1+ex)^(-4)))
}

ffdashdd.deriv = function(x){
  ex=exp(x)
  return(list(first=ex*(2*ex-1)*(1+ex)^(-4), second=ex*(7*ex-4*ex^2-1)*(1+ex)^(-5)))
}

#return estimated E(ff(X)) and var(ff(X)) by delta method
#INPUT: mx =E(X) and vx = Var(X)
ff.moments_exp = function(mx,vx){
  d = ff.deriv_exp(mx)
  return(list(mean=ff_exp(mx)+d$second*vx/2, var=(d$first^2*vx), meansq=ff_exp(mx)^2+(d$second*vx/2)^2+(d$first^2+ff_exp(mx)*d$second)*vx))
}

#return estimated E(ff'(X)) and var(ff'(X)) by delta method
#INPUT: mx =E(X) and vx = Var(X)
ffdash.moments = function(mx,vx){
  d = ffdash.deriv(mx)
  return(list(mean=ffdash(mx)+d$second*vx/2, var=(d$first^2*vx), meansq=ffdash(mx)^2+(d$second*vx/2)^2+(d$first^2+ffdash(mx)*d$second)*vx))
}

#return estimated E(ff''(X)) and var(ff''(X)) by delta method
#INPUT: mx =E(X) and vx = Var(X)
ffdd.moments = function(mx,vx){
  d = ffdd.deriv(mx)
  return(list(mean=ffdd(mx)+d$second*vx/2, var=(d$first^2*vx), meansq=ffdd(mx)^2+(d$second*vx/2)^2+(d$first^2+ffdd(mx)*d$second)*vx))
}

#return estimated E(ff'(X)*ff''(X)) and var(ff'(X)*ff''(X)) by delta method
#INPUT: mx =E(X) and vx = Var(X)
ffdashdd.moments = function(mx,vx){
  d = ffdashdd.deriv(mx)
  return(list(mean=ffdashdd(mx)+d$second*vx/2, var=(d$first^2*vx), meansq=ffdashdd(mx)^2+(d$second*vx/2)^2+(d$first^2+ffdashdd(mx)*d$second)*vx))
}






#### Addtionnal moment functions

fl = function(x,n){
  return(logit(x/n))
}


fl.deriv = function(x,n){
  return(list(first=n/(x*(n-x))))
}

fl.moments = function(mx,vx,n){
  d = fl.deriv(mx,n)
  return(d$first^2*vx)
}



ff = function(x){
  return(-log(1+exp(-x)))
}

ffdash = function(x){
  return(1/(1+exp(x)))
}

#first and second  derivatives of ff
ff.deriv = function(x){
  ex = exp(x)
  return(list(first = 1/(1+ex), second= -ex/(1+ex)^2))
}

#first and second  derivatives of ffdash
ffdash.deriv = function(x){
  ex=exp(x)
  return(list(first= -ex/(1+ex)^2, second= ex*(ex-1)*(1+ex)^(-3)))
}

#return estimated E(ff(X)) and var(ff(X)) by delta method
#INPUT: mx =E(X) and vx = Var(X)
ff.moments = function(mx,vx){
  d = ff.deriv(mx)
  return(list(mean=ff(mx)+d$second*vx/2, var=(d$first^2*vx), meansq=ff(mx)^2+(d$second*vx/2)^2+(d$first^2+ff(mx)*d$second)*vx))
}

#return estimated E(ff'(X)) and var(ff'(X)) by delta method
#INPUT: mx =E(X) and vx = Var(X)
ffdash.moments = function(mx,vx){
  d = ffdash.deriv(mx)
  return(list(mean=ffdash(mx)+d$second*vx/2, var=(d$first^2*vx), meansq=ffdash(mx)^2+(d$second*vx/2)^2+(d$first^2+ffdash(mx)*d$second)*vx))
}








#' Compute posterior mean and var for log(p), log(q), log(p0/p1) and log(q0/q1).
#'
#' This function returns posterior means and variances of log(p), log(q), log(p0/p1) and log(q0/q1) as lp, lq, lpratio and lqratio, respectively, where p
#' is the probability of going left and q=1-p.
#' @return a list with elements "lp.mean", "lp.var", "lq.mean", "lq.var" [and "lpratio.mean", "lpratio.var", "lqratio.mean", "lqratio.var"
#' if covariate is present, i.e. \code{g} is not NULL]
#' @keywords internal
compute.res <- function(zdat.ash.intercept, repara, baseline=NULL, w=NULL, g=NULL, zdat=NULL, zdat.ash=NULL){
  alpha=list(mean=zdat.ash.intercept$PosteriorMean, var=zdat.ash.intercept$PosteriorSD^2) #find mean and variance of alpha
  alpha.prior=list(mean=rep(0,length(zdat.ash.intercept$PosteriorMean)), var=sum(zdat.ash.intercept$fitted.g$pi*zdat.ash.intercept$fitted.g$sd^2))
  if(is.null(g)){#if covariate is absent
    lp = ff.moments(alpha$mean, alpha$var)
    lq = ff.moments(-alpha$mean, alpha$var)  #find mean and variance of log(q)
    lp.prior = ff.moments(alpha.prior$mean, alpha.prior$var)
    lq.prior = ff.moments(-alpha.prior$mean, alpha.prior$var)
    return(list(lp.mean=lp$mean, lp.var=lp$var, lq.mean=lq$mean, lq.var=lq$var, lp.prior.mean=lp.prior$mean, lp.prior.var=lp.prior$var, lq.prior.mean=lq.prior$mean, lq.prior.var=lq.prior$var))
  }else{#if covariate is present
    if(repara==TRUE){  #if reparametrization is used then we want gamma returned as well
      mbvar.ind=is.na(zdat[5,])
      mbvar=zdat[5,]
      mbvar[mbvar.ind]=0
    }else{
      mbvar=0
    }
    if(baseline=="grp")
      w1=w[1]
    else if (baseline=="inter")
      w1=0
    else
      w1=baseline
    #apply ash to vector of slope estimates and SEs
    zdat.ash_post=posterior_dist(zdat.ash$fitted.g,zdat[3,],zdat[4,])
    #compute the posterior third and fourth moments of beta
    gamma=list(mean=alpha$mean+(w1+mbvar)*zdat.ash$PosteriorMean, var=alpha$var+((w1+mbvar)*zdat.ash$PosteriorSD)^2)
    gamma.prior=list(mean=alpha.prior$mean, var=alpha.prior$var+(w1+mbvar)^2*sum(zdat.ash$fitted.g$pi*zdat.ash$fitted.g$sd^2))
    lp = ff.moments(gamma$mean, gamma$var)    #find mean and variance of p in baseline estimate
    lq = ff.moments(-gamma$mean, gamma$var)  #find mean and variance of q in baseline estimate
    lp.prior = ff.moments(gamma.prior$mean, gamma.prior$var)
    lq.prior = ff.moments(-gamma.prior$mean, gamma.prior$var)

    beta.tm=tfmoment(zdat.ash_post$mu, zdat.ash_post$sigma,3,zdat.ash_post$pi)
    beta.fm=tfmoment(zdat.ash_post$mu, zdat.ash_post$sigma,4,zdat.ash_post$pi)
    wSq=(w[2]+mbvar)^2-(w[1]+mbvar)^2
    PosteriorSq=zdat.ash$PosteriorSD^2+zdat.ash$PosteriorMean^2
    lpratio=ffwrapper(alpha$mean, alpha$var, zdat.ash$PosteriorMean, PosteriorSq, wSq, beta.tm, beta.fm, g)
    lqratio=ffwrapper(-alpha$mean, alpha$var, -zdat.ash$PosteriorMean, PosteriorSq, wSq, -beta.tm, beta.fm, g)
    lpratio.prior=ffwrapper(alpha.prior$mean, alpha.prior$var, zdat.ash$PosteriorMean, PosteriorSq, wSq, beta.tm, beta.fm, g)
    lqratio.prior=ffwrapper(-alpha.prior$mean, alpha.prior$var, -zdat.ash$PosteriorMean, PosteriorSq, wSq, -beta.tm, beta.fm, g)
    return(list(lp.mean=lp$mean, lp.var=lp$var, lq.mean=lq$mean, lq.var=lq$var, lpratio.mean=lpratio$mean, lpratio.var=lpratio$var, lqratio.mean=lqratio$mean, lqratio.var=lqratio$var, lp.prior.mean=lp.prior$mean, lp.prior.var=lp.prior$var, lq.prior.mean=lq.prior$mean, lq.prior.var=lq.prior$var, lpratio.prior.mean=lpratio.prior$mean, lpratio.prior.var=lpratio.prior$var, lqratio.prior.mean=lqratio.prior$mean, lqratio.prior.var=lqratio.prior$var))
  }
}
