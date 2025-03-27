rm(list=ls())

load("C:/Document/Serieux/Travail/Package/mvPoisVA/simredoing_mv_sd_0.RData")
res= do.call(rbind, res_list)
res=as.data.frame(res)
apply(res, 2,mean)
library(dplyr)
 est0 =res%>%
  group_by(L) %>%
  summarise(cov_pois= mean(false_cs_pois),
            cov = mean(false_cs ),
            n = n())

est0$overdisp =0


load("C:/Document/Serieux/Travail/Package/mvPoisVA/simredoing_mv_sd_0.5.RData")

res= do.call(rbind, res_list)
load("C:/Document/Serieux/Travail/Package/mvPoisVA/simredoing_mv.RData")
res= rbind(res,
           do.call(rbind, res_list))

res=as.data.frame(res)
apply(res, 2,mean)
library(dplyr)
est05 =res%>%
  group_by(L) %>%
  summarise(cov_pois= mean(false_cs_pois),
            cov = mean(false_cs ),
             n = n())
est05$overdisp =0.5




load("C:/Document/Serieux/Travail/Package/mvPoisVA/simredoing_mv_sd_1.RData")

res= do.call(rbind, res_list)
res=as.data.frame(res)
apply(res, 2,mean)
library(dplyr)
est1 =res%>%
  group_by(L) %>%
  summarise(cov_pois= mean(false_cs_pois),
            cov = mean(false_cs ),
            n = n())

est1$overdisp =1




df=data.frame(L= rep(est0$L,6),
              cov= 1- c(est0$cov_pois,est0$cov,
                     est05$cov_pois,est05$cov,
                     est1$cov_pois,est1$cov
                     ),
              overdisp= c(rep(est0$overdisp,2),
                          rep(est05$overdisp,2),
                          rep(est1$overdisp,2)
                          ),
              method= c( rep("PoisFsusie", length(est0$L)),rep("Fsusie", length(est0$L)),
                         rep("PoisFsusie", length(est0$L)),rep("Fsusie", length(est0$L)),
                         rep("PoisFsusie", length(est0$L)),rep("Fsusie", length(est0$L))
                         )
              )
library(ggplot2)
ggplot(df, aes(x=L, y=cov, colour = as.factor(method)))+
  geom_point(size=2)+
  geom_hline(yintercept = 0.95)+
  facet_wrap(overdisp~.)
