library(susieR)

set.seed(12345)
n=2^9
sigma=1
data <- readRDS("C:/Document/Serieux/Travail/Data_analysis_and_papers/anjing_data/xiong_atact_seq_multi_trait.12regions.51.2kb_expanded.binned_coverage_count.rds")
X=data[[1]]$filtered_geno
sum(is.na(X))

X
for ( k in 1:ncol(X)){
  if ( sum(is.na(X[,k]))>0){
    X[which(is.na(X[,k])),k]= median(X[-which(is.na(X[,k])),k])
  }
}

X= X[,-which (apply(X,2,sum)<10) ]


true_pos=1

dl= list()


mu0= mu=c(rep(0.3,n/4), rep(3, n/4), rep(10, n/4), rep(0.3, n/4))
  mu1= mu=c(rep(0.3,n/4), rep(3, n/4), rep(5, n/4), rep(0.3, n/4))
  mu2=  mu=c(rep(0.3,n/4), rep(3, n/4), rep(2.5, n/4), rep(0.3, n/4))
for ( i in 1:nrow(X)){

  if(X[i,true_pos]==0){
    dl[[i]]=  rpois(length(mu0),exp(log(mu0)+rnorm(n=length(mu),sd=sigma)))
  }
  if(X[i,true_pos]==1){
    dl[[i]]=  rpois(length(mu1),exp(log(mu1)+rnorm(n=length(mu),sd=sigma)))
  }
  if(X[i,true_pos]==2){
    dl[[i]]=  rpois(length(mu2),exp(log(mu2)+rnorm(n=length(mu),sd=sigma)))
  }

}

  Y= do.call(rbind, dl)

  library(fsusieR)
  res= acc_Pois_fSuSiE2(Y=log1p(Y),X=X, post_processing = "TI" )



  plot(res$fitted_func[[1]])
lines(res$cred_band[[1]][1,])
lines(res$cred_band[[1]][2,])




res=  acc_Pois_fSuSiE2 (Y= (Y),X=X , post_processing = "smash" )



plot(res$fitted_func[[1]])
lines(res$cred_band[[1]][1,])
lines(res$cred_band[[1]][2,])

plot(res$susiF.obj$fitted_func[[1]])
lines(res$susiF.obj$cred_band[[1]][1,])
lines(res$susiF.obj$cred_band[[1]][2,])

