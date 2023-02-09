library(ashr)
library(wavethresh)
library(susiF.alpha)
library(mvPoisVA)

set.seed(1)
#Example using curves simulated under the Mixture normal per scale prior

N <- 100    #Number of individuals
P <- 10     #Number of covariates/SNP
pos1 <- 1   #Position of the causal covariate for effect 1
pos2 <- 5   #Position of the causal covariate for effect 2
lev_res <- 7#length of the molecular phenotype (2^lev_res)
f1 <- sim_intenisty(lev_res )$sim_intens[-1]#first effect
f2 <- sim_intenisty(lev_res )$sim_intens[-1]#second effect

plot( f1, type ="l", ylab="effect", col="blue")
abline(a=0,b=0)
lines(f2, type="l", col="green")

legend(x=100,
       y=30,
       lty = rep(1,2),
       legend= c("effect 1", "effect 2" ),
       col=c( "blue","yellow"))
G = matrix(sample(c(0, 1,2), size=N*P, replace=TRUE), nrow=N, ncol=P) #Genotype
beta0       <- 0
beta1       <- 1
beta2       <- 1
count.data  <- list()

for ( i in 1:N)
{
  predictor <-beta1*G[i,pos1]*f1 + beta2*G[i,pos2]*f2
  count.data [[i]] <-   rpois(n= length(f1) ,
                              lambda =predictor  )

}
count.data <- do.call(rbind, count.data)




plot( count.data[1,], type = "l", col=(G[1, pos1]*3+1),
      main="Observed curves \n colored by the causal effect", ylim= c(-1,100), xlab="")
for ( i in 2:N)
{
  lines( count.data[i,], type = "l", col=(G[i, pos1]*3+1))

}
legend(x=0.3,
       y=-20,
       lty = rep(1,3),
       legend= c("0", "1","2"),
       col=c("black","blue","yellow"))



Y <-count.data
X <- G




reflect =FALSE
verbose=TRUE
n_gh = 10

tol= 1e-3
tol_vga_pois=1e-5
maxit=100
maxit.mrash=100
tol.mrash=1e-3
cal_obj.mrash=FALSE
verbose.mrash=FALSE
control_mixsqp=  list(verbose=FALSE,
                      eps = 1e-6,
                      numiter.em = 4
)
thresh_lowcount=0
prior_mv=  "mixture_normal_per_scale"
gridmult=sqrt(2)



## static param  ----

nullweight.mrash=10

init_pi0_w.mrash=10

X <- susiF.alpha:::colScale(X)
gh_points = fastGHQuad::gaussHermiteData(n_gh)

##initiatilzation for count data -----


Y_org <- Y
### dealing with non 2^S data ----
J = log2(ncol(Y)); if((J%%1) != 0) reflect=TRUE
if(reflect){
  tl <- lapply(1:nrow(Y), function(i) reflect_vec(Y[i,]))
  Y <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$x))
  idx_out <- tl[[1]]$idx #### indx of interest at the end
}

#### to avoid 0 in Y_min to correct at the end
Y <- Y+1


indx_lst <-  susiF.alpha::gen_wavelet_indx(log2(ncol(Y)))
### Wavelet like transform -----
tl <-  lapply(1:nrow(Y), function(i)
  get_empirical_intensity(Y[i,],
                          indx_lst = indx_lst)
)
### Cal Ymin Ytot -----
Y_min <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$Y_min))
Y_tot <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$Y_tot))
rm(tl)

emp_prop <- Y_min[1,]/Y_tot[1,]
emp_log_prop <- log(emp_prop[- length(emp_prop)])
plot(emp_prop)
emp_lambda_tot <- Y_tot[1,ncol(Y_tot)]
emp_log_lambda_tot <- log(emp_lambda_tot )




emp_log_prop <- log(sigmoid(EBmvFR.obj$fitted_wc[[1]][5,-256]))
emp_log_lambda_tot <-8# EBmvFR.obj$fitted_wc[[1]][5, 256]

lp <- emp_log_prop
lq = log(1-pmin(exp(lp),1-1e-10))# complementary prop

est <- emp_log_lambda_tot

J=log2(length(Y_min[1,]))
out <- rep( emp_log_lambda_tot, 2^J)
for(s in (J ):1){
  #print(exp(est))
  #readline("press a key")
  nD = 2^(J-s+1)
  nDo2 = nD/2
  tt <-1

  for(l in 0:(2^(s-1)-1)){
    ind = (l*nD+1):((l+1)*nD) # all "sub index for coef s,l (here s=D)
   # print(ind)
    ind_l <-  ind[1:nDo2] #all "sub index in the left for coef s,l (here s=D)
    ind_r <-  ind[(nDo2+1):nD] # all "sub index in the right for coef s,l (here s=D)

    out[ind_l] <- out[ind_l]+ lp[indx_lst[[(s )]][tt]]
    print(lp[indx_lst[[(s+1)]][tt]])
    out[ind_r] <- out[ind_r]+ lq[indx_lst[[(s )]][tt]]
    tt <- tt+1
  }
}

plot( exp(out)[idx_out])
plot(f1[idx_out])
lines(f1)
plot( exp(out)[idx_out])
lines(Y_org[1,])


tt <- reverse_intensity_transform(vec_int = c(emp_prop[-length(emp_prop)] ,
                                              emp_lambda_tot),
                                  indx_lst = indx_lst,
                                  is.prob =TRUE,
                                  is.log_int = FALSE)

plot( tt)
lines(Y[1,])
