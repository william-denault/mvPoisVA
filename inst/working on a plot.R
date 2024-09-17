 load("C:/Document/Serieux/Travail/Data_analysis_and_papers/enhancer_fine_mapping/results_fsusie/HLA-DRB5_locus_6_ATACseq.RData")

plot( (apply(out$Y,2,mean)))


X <- out$X
Y= out$Y
library(mvPoisVA)
fit  = acc_Pois_fSuSiE2(Y=Y, X=X[,100:150])


plot( (apply(out$Y,2,mean)))

for ( l in 1:length(fit$susiF.obj$fitted_func)){
  lines( fit$susiF.obj$fitted_func[[l]]  ,col=l+1)
}
lines(exp(fit$susiF.obj$ind_fitted_func[l,])-1,col=l+1)
plot(exp(fit$susiF.obj$fitted_func[[1]])-1,col=l+1)




fit  = acc_Pois_fSuSiE2(Y=Y, X=X[,1:50])

plot( (apply( Y,2,mean)))

for (  i in 1:nrow(fit$susiF.obj$ind_fitted_func) ){
  lines( exp(fit$susiF.obj$ind_fitted_func[i,] ) ,col=i+1)
}



