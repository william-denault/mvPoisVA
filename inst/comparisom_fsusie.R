
path_mv_pois="C:/Document/Serieux/Travail/Data_analysis_and_papers/GTEX_analysis_Fsusie_res/new_pois_split_pois"
lf_mv_pois= list.files(path_mv_pois)
lf_mv_pois[1:10]

extract_gene_id <- function(x) {
  sub("\\..*", "", x) # Removes everything after the first dot
}

lt =list()
n_cs_fsusie=0
overlapp= 0
for ( o in 1: length (lf_mv_pois)){


  load( paste0(path_mv_pois,"/",lf_mv_pois[o] ))


    cs_fsusie =out$res_fsusie$cs
    n_cs_fsusie= length(out$res_fsusie$cs)

    cs_size_fsusie= mean (lengths(cs_fsusie))





    n_cs_mv_pois= length(out$res_pois $cs)

    cs_mv_pois=out$res_pois $cs
    cs_size_pois= mean (lengths(cs_mv_pois))

  overlapp= 0
  if( n_cs_fsusie>0 & n_cs_mv_pois>0 ){



    for ( i in 1: length(  cs_fsusie)){
      for ( j in 1:length(cs_mv_pois)){

        if(  length( intersect(cs_fsusie[[i]] , cs_mv_pois[[j]]))>0){
          overlapp=overlapp+1
        }
      }



    }
  }




  lt[[o]]= c( n_cs_fsusie,cs_size_fsusie,n_cs_mv_pois,  cs_size_pois, overlapp)


}





sum_res= do.call(rbind, lt)



colnames(sum_res)= c( "n_cs_fsusie","cs_size_fsusie","n_cs_mv_pois",  "cs_size_pois", "overlapp")
sum_res= data.frame(sum_res)
sum_res$n_cs_mv_pois[which(sum_res$cs_size_pois>100)]=0
sum_res$n_cs_fsusie[which(sum_res$cs_size_fsusie>100)]=0
idx= which (sum_res$n_cs_fsusie==0 & sum_res$n_cs_mv_pois==0)
length(idx)


#sum_res_res=sum_res[-idx,]
table(sum_res$overlapp)

table(sum_res$overlapp[-idx])[1]


1- (table(sum_res$overlapp[-idx])[1]/(sum(sum_res$overlapp[-idx])+table(sum_res$overlapp[-idx])[1]))


idx= which (sum_res$n_cs_fsusie==0 | sum_res$n_cs_mv_pois==0)
apply(sum_res[-idx,], 2, sum)
apply(sum_res[-idx,], 2, mean)

length(which(sum_res[,1]==0))
length(which(sum_res[,3]==0))

which (sum_res$n_cs_fsusie==0 & !(sum_res$n_cs_mv_pois==0))

source("C:/Document/Serieux/Travail/Data_analysis_and_papers/GTEX_analysis_Fsusie/code/plot_log.R")

# example where mv POIS found effect -----


which (sum_res$n_cs_fsusie==0 & !(sum_res$n_cs_mv_pois==0))
load( paste0(path_mv_pois,"/",lf_mv_pois[32] ))



fsusie_log_plot(obj=out$res_pois ,chr = paste0("chr",out$chr),
                pos0 = out$locus[1],pos1 = out$locus[2],
                X=out$X,Y=out$Y,snp_info = out$info_SNP,cs = 1,
                effect_log = TRUE,
                log1p_count=TRUE
)


fsusie_log_plot(obj=out$res_fsusie ,chr = paste0("chr",out$chr),
                pos0 = out$locus[1],pos1 = out$locus[2],
                X=out$X,Y=out$Y,snp_info = out$info_SNP,cs = 1,
                effect_log = TRUE,
                log1p_count=TRUE
)

out$res_pois$cs

out$res_fsusie$cs

hist(c(cor(out$X[,out$res_fsusie$cs[[1]] ])))

hist(c(cor(out$X[,out$res_pois$cs[[1]] ])))
plot(out$res_fsusie$pip )
plot(out$res_pois$pip )




size_factor=rowSums(out$Y)/ mean(rowSums(out$Y))

plot (log1p( out$Y /as.vector(size_factor)) ,out$Mu_pm)
for (i in 1:nrow(out$Y)){

  points (log1p( out$Y /as.vector(size_factor))[i,] ,out$Mu_pm[i,],
        col=i)
}
abline(a=0,b=1)


# 5 -----


load( paste0(path_mv_pois,"/",lf_mv_pois[5] ))



fsusie_log_plot(obj=out$res_pois ,chr = paste0("chr",out$chr),
                pos0 = out$locus[1],pos1 = out$locus[2],
                X=out$X,Y=out$Y,snp_info = out$info_SNP,cs = 1,
                effect_log = TRUE,
                log1p_count=TRUE
)



size_factor=rowSums(out$Y)/ mean(rowSums(out$Y))

plot (log1p( out$Y /as.vector(size_factor)) ,out$Mu_pm)
for (i in 1:nrow(out$Y)){

  points (log1p( out$Y /as.vector(size_factor))[i,] ,out$Mu_pm[i,],
          col=i)
}
abline(a=0,b=1)







#2  5  8 15 24 25 33
# 15 -----


load( paste0(path_mv_pois,"/",lf_mv_pois[33] ))


out$res_pois$cs
plot(out$res_pois$pip )
out$res_fsusie$cs
fsusie_log_plot(obj=out$res_pois ,chr = paste0("chr",out$chr),
                pos0 = out$locus[1],pos1 = out$locus[2],
                X=out$X,Y=out$Y,snp_info = out$info_SNP,cs = 1,
                effect_log = TRUE,
                log1p_count=TRUE
)



size_factor=rowSums(out$Y)/ mean(rowSums(out$Y))

plot (log1p( out$Y /as.vector(size_factor)) ,out$Mu_pm)
for (i in 1:nrow(out$Y)){

  points (log1p( out$Y /as.vector(size_factor))[i,] ,out$Mu_pm[i,],
          col=i)
}
abline(a=0,b=1)














load( paste0(path_mv_pois,"/",shared_genes[o], ".csv.gz.RData"))
out$res$susiF$cs

fsusie_log_plot(obj=out$res$susiF.obj,chr = paste0("chr",out$chr),
                pos0 = out$locus[1],pos1 = out$locus[2],
                X=out$X,Y=out$Y,snp_info = out$info_SNP,cs = 1,
                effect_log = TRUE,
                log1p_count=TRUE
)





rm(list=ls())
path_mv_pois="C:/Document/Serieux/Travail/Data_analysis_and_papers/GTEX_analysis_Fsusie_res/new_pois_split_pois"
lf_mv_pois= list.files(path_mv_pois)
lf_mv_pois[1:10]
load( paste0(path_mv_pois,"/",lf_mv_pois[24] ))
library(fsusieR)

size_factor=rowSums(out$Y)/ mean(rowSums(out$Y))
res1= susiF(X= out$X, Y=log1p( out$Y /as.vector(size_factor)),
            L=3)


tt <- mvPoisVA::pois_mean_split(c(out$Y),s= rep(  size_factor, ncol(out$Y)))
Mu_pm <- matrix( tt$posterior$mean_log,byrow = FALSE, ncol=ncol(out$Y))


plot(log1p( out$Y /as.vector(size_factor)) , Mu_pm  )
for (i in 1:nrow(out$Y)){

  points (log1p( out$Y /as.vector(size_factor))[i,] , Mu_pm[i,],
          col=i)
}
abline(a=0,b=1)
res2= susiF(X= out$X, Y= Mu_pm,
            L=3)
res1$cs
res2$cs
plot(res1$fitted_func[[1]])
plot(res1$fitted_func[[2]])


plot(res2$fitted_func[[1]])
plot(res2$fitted_func[[3]])



