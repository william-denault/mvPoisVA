library(ashr)
library(wavethresh)
library(susiF.alpha)
library(mvPoisVA)
init=TRUE
set.seed(1)
#Example using curves simulated under the Mixture normal per scale prior
check=1
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

if(verbose){
  print("done transforming data")
}

###initialization for functional reg object -----
lowc_wc <-   susiF.alpha::which_lowcount(Y_f=Y_min ,thresh_lowcount)

###Initialization dynamic param
Mu_pm = logit((Y_min/Y_tot) ) #remove last column contain C coeff
Mu_pm[Mu_pm==-Inf] =  logit(0.1)
Mu_pm[Mu_pm==Inf]  =  logit(0.9)

Mu_pm[,ncol(Y_min)] <- log(Y_min[,ncol(Y_min)])
Mu_pv = 1/Y_tot

  b_pm <- 0*Y

sigma2_bin  = 1
sigma2_pois = 1

iter =1
l_bpm <- list()
while( check >tol & iter <5){

  #### Check potential pb due to centering
  post_mat <- get_post_log_int(Mu_pm       = Mu_pm,
                               Mu_pv       = Mu_pv,
                               Y_min       = Y_min,
                               Y_tot       = Y_tot,
                               sigma2_bin  = sigma2_bin,
                               sigma2_pois = sigma2_pois,
                               b_pm        = b_pm,
                               gh_points   = gh_points,
                               tol         = tol_vga_pois)
  if(verbose){
    print( paste('Posterior log intensity computed for iter ',iter))
  }
  Mu_pm <- post_mat$A_pm
  Mu_pv <- post_mat$A_pv


  tmp_Mu_pm <- Mu_pm - colMeans(Mu_pm, na.rm = TRUE) #potentially run smash on colmean
  W <- list( D = tmp_Mu_pm [, -ncol(tmp_Mu_pm )],
             C = tmp_Mu_pm [,  ncol(tmp_Mu_pm )])



  if(init){

    temp <- susiF.alpha:: init_prior(Y              = tmp_Mu_pm,
                                     X              = X,
                                     prior          = prior_mv ,
                                     v1             = v1,
                                     indx_lst       = indx_lst,
                                     lowc_wc        = lowc_wc,
                                     control_mixsqp = control_mixsqp,
                                     nullweight     = nullweight.mrash,
                                     gridmult       = gridmult )
    G_prior     <- temp$G_prior


    #Recycled for the first step of the while loop
    EBmvFR.obj   <-  susiF.alpha::init_EBmvFR_obj(G_prior = G_prior,
                                                  Y       = Y,
                                                  X       = X
                                                  )

    init=FALSE
  }
  ### TODO: Maybe use better restarting point for EBmvFR.obj
  EBmvFR.obj   <- susiF.alpha::EBmvFR.workhorse(EBmvFR.obj     = EBmvFR.obj,
                                                W              = W,
                                                X              = X,
                                                tol            = tol.mrash,
                                                lowc_wc        = lowc_wc  ,
                                                init_pi0_w     = init_pi0_w.mrash ,
                                                control_mixsqp = control_mixsqp ,
                                                indx_lst       = indx_lst,
                                                nullweight     = nullweight.mrash,
                                                cal_obj        = cal_obj.mrash,
                                                verbose        = verbose.mrash,
                                                maxit          = maxit.mrash
                                                )
  if(verbose){
    print( paste('Posterior of regression coefficient computed for iter ',iter))
  }
  b_pm <-   X%*%EBmvFR.obj$fitted_wc[[1]]+ colMeans(Mu_pm, na.rm = TRUE)
  l_bpm[[iter]] =b_pm
  iter=iter+1
  ##include mr.ash
}


fitted_proc <- lapply(1:nrow( b_pm ),
                      function( i) reverse_intensity_transform(vec_int = b_pm[i,],
                                                             indx_lst = indx_lst) )
fitted_proc <- do.call(rbind,fitted_proc)

fitted_effect <- lapply(1:nrow( EBmvFR.obj$fitted_wc[[1]] ),
                       function( i) reverse_intensity_transform(vec_int =EBmvFR.obj$fitted_wc[[1]][i,],
                                                                indx_lst = indx_lst) )
 fitted_effect <- do.call(rbind,fitted_effect)
plot( fitted_effect[5,])
