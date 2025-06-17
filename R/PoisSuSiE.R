



#' @export



library(susieR)

X_all= N3finemapping$X
X= X_all[sample(size=200, 1:nrow(X_all)),]
L= sample(size=1, 1:3)

true_pos= sample(size=L, 1:ncol(X))
pos_causal= rep(0, ncol(X))
pos_causal[true_pos]=1

beta <- rep(0,ncol(X))
beta[true_pos] <-1
log_mu <- X %*% beta+ rnorm( nrow(X), sd=0.5)
y <- rpois(n=nrow(X), lambda = exp(log_mu))
Z= matrix(rnorm(200*10), ncol=10)
library(mr.ash.alpha)
Pois_SuSiE <- function(y,
                        Z,
                        X,
                        L=3,
                        scaling= NULL,
                        L_start=3,
                        verbose=TRUE,
                        tol= 1e-3,
                        tol_vga_pois=1e-5,
                        max.iter=100,
                        control_mixsqp=  list(verbose=FALSE,
                                              eps = 1e-6,
                                              numiter.em = 4
                        )


)
{
  ####Changer les calcul d'objective -----
  if(missing(X)&missing(Z)){
    stop("Please provide a Z or a X matrix")
  }

  fit_approach <- "both"
  if(missing(X)){
    print("No correlated covariate provided, the algorithm will perform penalized regression only")
    fit_approach <- "penalized"
  }
  if(missing(Z)){
    print("No Z matrix provided,  the algorithm will perform fine-mapping only")
    fit_approach <- "fine_mapping"

  }




  if(fit_approach %in% c('both',"fine_mapping")){
    tidx <- which(apply(X,2,var)==0)
    if( length(tidx)>0){
      warning(paste("Some of the columns of X are constants, we removed" ,length(tidx), "columns"))
      X <- X[,-tidx]
    }

  }

  if(fit_approach %in% c('both',"penalized")){
    tidx <- which(apply(Z,2,var)==0)
    if( length(tidx)>0){
      warning(paste("Some of the columns of Z are constants, we removed" ,length(tidx), "columns"))
      Z <- Z[,-tidx]
    }
  }

  if( is.null( scaling)){
    scaling = rep(1, length(y))
  }else{
    if( ! (length(scaling)== nrow(y))){
      warning(paste("scaling shoudl have a length equal to number of row  in y"  ))
    }
  }



  init_val_pois<- c(log(y+1))
  beta_pois <- 0* c(log(y+1))
  sigma2_pois=1

  ##initiatilzation for count data -----
  Mu_pm<- y
  iter=1
  beta_pois <- 0* c(log(Mu_pm +1))
  check <- 3*tol

  b_pm <- 0* Mu_pm
  fm_pm <- 0* Mu_pm
  fitted_mrash = 0*y
  fitted_susie = 0*y

  while( check >tol & iter <=  max.iter ){


    if ( iter ==1 ){
      tt= ebpm_normal(c(y),s=  scaling   )
      Mu_pm <-  tt$posterior$mean_log


    }else{



      tt <-    pois_mean_GG(c(y), prior_mean = c(Mu_pm_init),
                            prior_var = sigma2_pois )
      Mu_pm <-  tt$posterior$posteriorMean_latent
      Mu_pv <-  tt$posterior$posteriorVar_latent
    }



    if(verbose){
      print( paste('Posterior log intensity computed for iter ',iter))
    }

    if( fit_approach %in% c('both',"penalized")){

      partial_res  =  Mu_pm-fitted_susie
      fitmrash     =  mr.ash(X=Z,y=partial_res)
      fitted_mrash =  predict(fitmrash,Z)
    }else{
      fitted_mrash= 0*y
    }


    if( fit_approach %in% c('both',"fine_mapping")){
      partial_res= Mu_pm-fitted_mrash
      fitsusie= susie(X=X,y=partial_res,L=L)
      fitted_susie=  predict(fitsusie,X)
    }

    resid <- Mu_pm -fitted_susie - fitted_mrash
    #not correct to work on later
    sigma2_pois <- var( resid )
    #print(sigma2_pois)
    Mu_pm <-fitted_susie +fitted_mrash#update
    Mu_pm_init <-Mu_pm
    iter=iter+1
    ##include mr.ash

    #  par (mfrow=c(1,2))

    # plot ( y[1,], col="blue")
    #points ( exp(Mu_pm  [1,]))
    # lines(exp(Mu_pm  [1,]), col="green")


     plot( y ,exp(Mu_pm   ))

    #abline(a=0,b=1)
    #par (mfrow=c(1,1))
  }







  if( fit_approach ==   "both" )
  {
    out <- list( Mu_pm=Mu_pm,
                 susie.obj= fitsusie,
                 fitted_susie=fitted_susie,
                 mr.ash.obj= fitmrash,

                 fitted_mrash=fitted_mrash )
  }

  if( fit_approach ==   "fine_mapping" )
  {
    out<- list( Mu_pm=Mu_pm,
             susie.obj= fitsusie,
             fitted_susie=fitted_susie  )
  }
  if( fit_approach ==   "penalized")
  {
    list( Mu_pm=Mu_pm,
          mr.ash.obj= fitmrash,

          fitted_mrash=fitted_mrash )
  }
  return(out)

}





X_all= N3finemapping$X
X= X_all[sample(size=200, 1:nrow(X_all)),]
L= sample(size=1, 1:3)

true_pos= sample(size=L, 1:ncol(X))
pos_causal= rep(0, ncol(X))
pos_causal[true_pos]=1

beta <- rep(0,ncol(X))
beta[true_pos] <-1
log_mu <- X %*% beta+ rnorm( nrow(X), sd=0.5)
y <- rpois(n=nrow(X), lambda = exp(log_mu))
Z= matrix(rnorm(200*10), ncol=10)


res= Pois_SuSiE(X=X,y=y )
res$susie.obj$sets
res0= susie(X=X,y=log1p(y ))
res0$sets
true_pos
