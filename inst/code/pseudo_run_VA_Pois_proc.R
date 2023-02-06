
#data simulation
N=1000
Tp <- 63
Y <-matrix( rpois(N*Tp, lambda = 20), ncol=Tp)


### parameters/arguments
#b_pm_init start posterior mean
reflect =FALSE
verbose=TRUE
n_gh = 10 #nb points for Gauss Hermite quadrature


### data formating ------


b_pm = 0*Y
# if(missing(b_pm_init)){
#   b_pm = 0*Y
# }else{
#   b_pm = b_pm_init #posterior mean
# }


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
rm(tl)

if(verbose){
  print("done transforming data")
}
# get the matrix for Binomial reg (and Pois reg for top coefficient)





#deal with case exactly one or exactly 0

Mu_pm = logit((Y_min/Y_tot)[,-ncol(Y_min)]) #Y_min/Y_tot
#PB which 0 count bin
Mu_pm[which(is.nat(Mu_pm))] =logit(0.5)

Mu_pm[Mu_pm==-Inf] =  logit(0.1)
Mu_pm[Mu_pm==Inf] =  logit(0.9)
Mu_pv = 1/Y_tot





Mu_pm = logit(x/nb) #Y_min/Y_tot

Mu_pm[mu_pm==-Inf] = logit(0.1)
Mu_pm[mu_pm==Inf] = logit(0.9)
Mu_pv = rep(1/n,n)


# here everything is a vector (not sigma^2) mu_pm and mu_

#check line 52 VA_binomial
opt<- vga_binomial(c(mu_pm,log(mu_pv)),x,nb,b_pm,sigma2,gh_points=gh_points)
mu_pm = opt$m
mu_pv = opt$v
opt = vga_binomial(c(mu_pm,log(mu_pv)),x,nb,b_pm,sigma2,gh_points=gh_points)
mu_pm = opt$m
mu_pv = opt$v


b_pv = res$posterior$sd^2
H = res$log_likelihood + n*(log(2*pi*sigma2)/2)+sum((mu_pm^2-2*mu_pm*b_pm+b_pm^2+b_pv)/sigma2/2)

# Update sigma2
sigma2 = mean(mu_pm^2+mu_pv+b_pm^2+b_pv-2*b_pm*mu_pm)

