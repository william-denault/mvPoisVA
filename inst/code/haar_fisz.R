library(haarfisz)
library(ashr)
library(wavethresh)
library(susiF.alpha)
library(mvPoisVA)
init=TRUE
set.seed(1)
#Example using curves simulated under the Mixture normal per scale prior
check=1
N <- 200    #Number of individuals
P <- 10     #Number of covariates/SNP
pos1 <- 1   #Position of the causal covariate for effect 1
pos2 <- 7   #Position of the causal covariate for effect 2
lev_res <- 7#length of the molecular phenotype (2^lev_res)
f1 <- sim_intenisty(lev_res )$sim_intens #first effect
f2 <- sim_intenisty(lev_res )$sim_intens #second effect

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



lst <- list()
for ( i in 1:nrow(count.data)){
  lst[[i]] <-hft(count.data[i,])
}



hft_dat <- do.call( rbind,lst)
inv_lst <- list()
for ( i in 1:nrow(count.data)){
  inv_lst[[i]] <-hft.inv(hft_dat[i,])
}
plot( inv_lst[[1]], count.data[1,])


HF_susiF (Y=count.data, X=G)
