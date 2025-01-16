
library(mvPoisVA)

n=2^9
sigma=1.5
mu=c(rep(0.3,n/4), rep(30, n/4), rep(10, n/4), rep(0.3, n/4))
set.seed(1)
s=1.2
x = rpois(length(mu),s*exp(log(mu)+rnorm(n=length(mu),sd=sigma)))
x1=x
fit_scaled = pois_smooth_split(x,maxiter=30,s=s )
fit_ebps_scaled = ebps(x ,s=s )
fit_pois_scaled = pois_mean_split(x,s=s  )



set.seed(1)
s=1.2
x = rpois(length(mu),s*exp(log(mu)+rnorm(n=length(mu),sd=sigma)))
plot(mu, type="l", main="Underlying intensity")
fit_unscaled = pois_smooth_split(x/s,maxiter=30 )
fit_ebps_unscaled = ebps(x/s )
fit_pois_unscaled = pois_mean_split(x/s )
plot (log(fit_pois_scaled$posterior$mean),log(fit_pois_unscaled $posterior$mean), col= "magenta" )

points(log(fit_ebps_scaled$posterior$mean_smooth),
       log(fit_ebps_unscaled $posterior$mean_smooth),col= "orange2" )
points(( fit_scaled$Eb),  (fit_unscaled $Eb),col= "blue3" )
abline(a=0,b=1)



plot (log(fit_pois_scaled$posterior$mean),log(x/s+1), col= "magenta" )

points(log(fit_ebps_scaled$posterior$mean_smooth)
       ,log(x/s+1),col= "orange2" )
points(( fit_scaled$Eb),  log(x/s+1),col= "blue3" )
