rm(list = ls())
library(mvPoisVA)
library(fsusieR)
library(susieR)
library(ebnm)
data(N3finemapping)
X <- N3finemapping$X
mysd=0.2
N =50
genotype <-X[1:N,1:100]

idx <- which( apply( genotype,2, var ) <1e-15)
genotype <- genotype [, -idx]
lev_res =6
count.data  <- list()
L <-2# sample(1:2, size =1)#actual number of effect

lf <-  list()
lf[[1]]<- rep(0.1, 2^6)
lf[[1]][10:20] <-2

lf[[2]]<- rep(0.1, 2^6)
lf[[2]][50:60] <-2



mu= lf[[2]]
x = rpois(length(mu),exp(log(mu)+rnorm(n=length(mu),sd=mysd)))

fit = pois_smooth_split(x,maxiter=30)
plot(x,col='grey80')
lines(exp(fit$Eb))
lines( smashr::smash(haarfisz::hft(x))  , col="green")
lines(lf[[2]])

data(N3finemapping)
X <- N3finemapping$X
genotype <-X[sample(1:nrow(X), size=N),1:100]

idx <- which( apply( genotype,2, var ) <1e-15)
if( length(idx)==0){
  X <-genotype

  Rtrue <- cor (genotype )
}else{
  genotype <- genotype [, -idx]
  X <-genotype

}
G<- genotype
X <- (X -0.99*min(X))/(0.5*max(X ))

G <-  (G -0.99*min(G ))/(0.5*max(G ))

tpos <- sample(1:ncol(genotype), replace = FALSE,size=2)
true_pos <- tpos
pos1 <- tpos[1]
pos2 <- tpos[2]
if( length(which(apply(G,2,var)==0))>0){
  G <- G[,-which(apply(G,2,var)==0)]
}
# G <- matrix( rnorm(nrow(genotype)*300), nrow = nrow(genotype))


predictor <-rep (0, length(lf[[1]] ))
count.data  <- list()

G[ , true_pos[1]] <-G[ , true_pos[1]] -min(G[ , true_pos[1]] )
G[ , true_pos[2]] <-G[ , true_pos[2]] -min(G[ , true_pos[2]] )
for ( i in 1:N)
{

  predictor <-rep (0, length(lf[[1]] ))

  for ( l in 1:L){
    predictor <-predictor + G[i, true_pos[l]]*lf[[l]]+0.3
  }
  predictor <- exp( predictor + rnorm(  length(lf[[1]]), sd=mysd))

  count.data [[i]] <-   rpois(n= length(lf[[1]]) ,
                              lambda =predictor  )

}
count.data <- do.call(rbind, count.data)


Y <- count.data
image(Y)
hist(Y)


tt <- pois_mean_split(c(Y),
                      mu_pm_init= c(log(Y+1)))

Y_t <- matrix( tt$posterior$mean_log,byrow = FALSE, ncol=ncol(Y))


Y_t1  <-fit_latent_nugget(Y )

tl =list()

for ( i in 1:nrow(Y)){

  tl[[i]] =pois_mean_split(Y[i,],
                           mu_pm_init= c(log(Y[i,]+1)))$posterior$mean_log
}
Y2 = do.call(rbind, tl)


Y3 = matrix( pois_mean_split(c(Y),
                             mu_pm_init= c(log(Y+1)))$posterior$mean_log, byrow = FALSE, nrow = nrow(Y) )


plot( (Y2 +1), (Y3 +1) )
plot( (Y2 +1), log(Y+1) )
points( (Y3 +1), log(Y+1), col="magenta")

points( (Y_t1 $Y), log(Y+1), col="green")
points(Y_t  , log(Y+1), col="red")

abline(a=0,b=1)

