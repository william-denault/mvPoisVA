 set.seed(12345)
 n=2^9
 sigma=0.2
 mu=c(rep(0.3,n/4), rep(3, n/4), rep(10, n/4), rep(0.3, n/4))

est_ebps <- list()
est_smooth <- list()
 for ( k in 1:100 ){

   x = rpois(n,exp(log(mu)+rnorm(n,sd=sigma)))
   fit = pois_smooth_split(x,maxiter=30)
   fit_ebps = ebps(x)

   est_ebps [[k]] <- sqrt(fit_ebps$fitted_g$sigma2)
   est_smooth  [[k]] <-  sqrt(fit$sigma2)

   print(k)
 }
est_ebps <- do.call(c, est_ebps)
est_smooth <-  do.call(c, est_smooth )
df =data.frame( est=  c(est_ebps, est_smooth),
            method=rep(c("EBPS", "Poisson smoothing"),
                       each=length(est_smooth))
)
boxplot(est~method, data=df, main="Estimated over-dispersion standard deviation")
abline(h=sigma, col="red")



set.seed(12345)
n=2^9
sigma=0.5
mu=c(rep(0.3,n/4), rep(3, n/4), rep(10, n/4), rep(0.3, n/4))

est_ebps <- list()
est_smooth <- list()
for ( k in 1:100 ){

  x = rpois(n,exp(log(mu)+rnorm(n,sd=sigma)))
  fit = pois_smooth_split(x,maxiter=30)
  fit_ebps = ebps(x)

  est_ebps [[k]] <- sqrt(fit_ebps$fitted_g$sigma2)
  est_smooth  [[k]] <-  sqrt(fit$sigma2)

  print(k)
}
est_ebps <- do.call(c, est_ebps)
est_smooth <-  do.call(c, est_smooth )
df =data.frame( est=  c(est_ebps, est_smooth),
                method=rep(c("EBPS", "Poisson smoothing"),
                           each=length(est_smooth))
)
boxplot(est~method, data=df, main="Estimated over-dispersion standard deviation")
abline(h=sigma, col="red")




set.seed(12345)
n=2^9
sigma=1
mu=c(rep(0.3,n/4), rep(3, n/4), rep(10, n/4), rep(0.3, n/4))

est_ebps <- list()
est_smooth <- list()
for ( k in 1:100 ){

  x = rpois(n,exp(log(mu)+rnorm(n,sd=sigma)))
  fit = pois_smooth_split(x,maxiter=30)
  fit_ebps = ebps(x)

  est_ebps [[k]] <- sqrt(fit_ebps$fitted_g$sigma2)
  est_smooth  [[k]] <-  sqrt(fit$sigma2)

  print(k)
}
est_ebps <- do.call(c, est_ebps)
est_smooth <-  do.call(c, est_smooth )
df =data.frame( est=  c(est_ebps, est_smooth),
                method=rep(c("EBPS", "Poisson smoothing"),
                           each=length(est_smooth))
)
boxplot(est~method, data=df, main="Estimated over-dispersion standard deviation")
abline(h=sigma, col="red")


