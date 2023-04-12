library(ashr)
library(wavethresh)
library(susiF.alpha)
library(mvPoisVA)
init=TRUE
set.seed(1)
#Example using curves simulated under the Mixture normal per scale prior
check=1
N <- 20    #Number of individuals
P <- 10     #Number of covariates/SNP
pos1 <- 1   #Position of the causal covariate for effect 1
pos2 <- 7   #Position of the causal covariate for effect 2
lev_res <- 7#length of the molecular phenotype (2^lev_res)
f1 <- 0.1*sim_intenisty(lev_res )$sim_intens[-1]#first effect
f2 <- 0.1*sim_intenisty(lev_res )$sim_intens[-1]#second effect

f1 <- ifelse(f1>8, 8,f1)
f2 <- ifelse(f2>8, 8,f2)

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
X <-   G ###### Include intercept----


Z <-   matrix(sample(c(0, 1,2), size=N*P, replace=TRUE), nrow=N, ncol=P) #Genotype



out_susiF <- susiF(Y=Y,X=X, L=10)
out_susiF$cs

out_mvpois <- mv_Poisproc_reg(Y=Y,X=X,  L=10)
out_mvpois$cs

out_susiF <- HF_susiF(Y=Y,X=X, L=10)
out_susiF$cs
