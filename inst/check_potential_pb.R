load("C:/Document/Serieux/Travail/Data_analysis_and_papers/enhancer_fine_mapping/results_fsusie/MRPL18_locus_4_ATACseq.RData")

source("C:/Document/Serieux/Travail/Data_analysis_and_papers/enhancer_fine_mapping/analysis_res/fsusie_res_tools.R" )

plot(apply(out$Y,2,mean)~out$pos)

which.max(apply(out$Y,2,mean))
boxplot(out$Y[,913] ~
          out$X[, out$res$cs[[3]][1]])

library(mvPoisVA)

devtools::load_all(".")


tol=0.0000001
Y <- out$Y
Mu_pm<- Y
iter=1
beta_pois <- 0* c(log(Mu_pm +1))
check <- 3*tol
sigma2_pois=0.1
Mu_pm_init <- log(Mu_pm+1)
tt <- pois_mean_split(c(Y),
                      mu_pm_init= c(Mu_pm_init))
Mu_pm <- matrix( tt$posterior$mean_log,byrow = FALSE, ncol=ncol(Y))

par(mfrow=c(2,1))
plot(apply(Mu_pm,2,mean)~out$pos)
plot(apply(log(out$Y+1),2,mean)~out$pos)

plot(apply(out$Y,2,mean)~out$pos)



out$res$cs
X <- out$X[, which(colnames(out$X) %in% do.call(c,
                                                lapply(1:length(out$res$cs), function(l) names(out$res$cs[[l]])))
                   )]
X

library(fsusieR)
res<- acc_Pois_fSuSiE2 (X=X, Y=out$Y, L=10, tol=1e-6 , verbose=TRUE)
res_fsusie<- fsusieR:::susiF(X=X, Y=Mu_pm ,
                             L=10,
                             tol=1e-6 ,
                             verbose=TRUE,
                             cal_obj   = TRUE,
                             maxit = 20,
                             cor_small = TRUE,
                             quantile_trans = TRUE,
                             post_processing = "HMM")

diff(res_fsusie$ELBO)
plot(res_fsusie$ELBO)
plot_susiF(res_fsusie)
plot_susiF(out$res)

plot_fsusie(res_fsusie)

tt_out<-out
tt_out$res <- res_fsusie

plot_fsusie(out)
plot(res_fsusie$ fitted_func[[1]])
plot(out$res$ fitted_func[[2]])

out$res$cs
res_fsusie$cs
