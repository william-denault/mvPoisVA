



#' @export


acc_Pois_fSuSiE2 <- function(Y,
                            Z,
                            X,
                            L=3,
                            scaling= NULL,
                            ebps_method=c('pois_mean_split',
                                          'ind_pois_mean_split',
                                          'ind_ebps',
                                          'ind_poisson_smoothing',
                                           'nugget'),
                            L_start=3,
                            reflect =FALSE,
                            verbose=TRUE,
                            n_gh = 10,
                            init_b_pm,
                            tol= 1e-3,
                            tol_vga_pois=1e-5,
                            maxit=10,
                            control_mixsqp=  list(verbose=FALSE,
                                                  eps = 1e-6,
                                                  numiter.em = 4
                            ),
                            thresh_lowcount= 1e-2,
                            prior_mv=  "mixture_normal",
                            gridmult=sqrt(2),
                            nullweight.mrash=10,
                            init_pi0_w.mrash=10,
                            cov_lev=0.95,
                            min_purity     =0.5,
                            greedy=TRUE,
                            backfit=TRUE,
                            tol.mrash=1e-3,
                            verbose.mrash=TRUE,
                            maxit.mrash=10,
                            cal_obj.mrash=FALSE,
                            maxit.fsusie=100,
                            cal_obj.fsusie=FALSE,
                            max_SNP_EM     = 100,
                            max_step_EM    = 1,
                            cor_small=TRUE,
                            post_processing="HMM"

)
{
  ####Changer les calcul d'objective -----
  if(missing(X)&missing(Z)){
    stop("Please provide a Z or a X matrix")
  }

  #static arguements
  ebps_method  <- match.arg(ebps_method)
  fit_approach <- "both"


  if (missing(X) & missing(Z)){
    stop("Please provide a Z or a X matrix")
  }
  if(missing(X)){
    print("No covariate to fine-map provided, the algorithm will perform penalized regression only")
    fit_approach <- "penalized"
  }
  if(missing(Z)){
    print("No Z matrix provided,  the algorithm will perform fine-mapping only")
    fit_approach <- "fine_mapping"

  }
  ##initiatilzation -----
  init=TRUE
  J = log2(ncol(Y)); if((J%%1) != 0) reflect=TRUE
  if(reflect){
    tl <- lapply(1:nrow(Y), function(i) reflect_vec(Y[i,]))
    Y <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$x))
    idx_out <- tl[[1]]$idx #### indx of interest at the end
  }else{
    idx_out <- 1: ncol(Y)
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
    scaling = rep(1, nrow(Y))
  }else{
    if( ! (length(scaling)== nrow(Y))){
      warning(paste("scaling shoudl have a length equal to number of row  in Y"  ))
    }
  }

  indx_lst <-  fsusieR::gen_wavelet_indx(log2(ncol(Y)))






  init_val_pois<- c(log(Y+1))
  beta_pois <- 0* c(log(Y+1))
  sigma2_pois=1

  ##initiatilzation for count data -----
  Mu_pm<- Y
  iter=1
  beta_pois <- 0* c(log(Mu_pm +1))
  check <- 3*tol
  sigma2_pois=0.1
  Mu_pm_init <- log(Mu_pm+1)
  ##### Poisson Part ----

  if (ebps_method =="pois_mean_split" ){
    tt <- pois_mean_split(c(Y),s= rep( scaling, ncol(Y)),
                          mu_pm_init= c(Mu_pm_init))
    Mu_pm <- matrix( tt$posterior$mean_log,byrow = FALSE, ncol=ncol(Y))
  }
  if(ebps_method =="ind_pois_mean_split" ){

   Mu_pm <-0*Y
    for ( i in 1:nrow(Y)){
     Mu_pm[i,]     =   pois_mean_split(Y[i,],s=scaling[i],
                                       mu_pm_init= c(log(Y[i,]+1)) )$posterior$mean_log
    }
  }

  if(ebps_method =="ind_ebps" ){


    Mu_pm <-0*Y
    for ( i in 1:nrow(Y)){
      Mu_pm[i,]     =  ebps(Y[i,],s=scaling[i])$posterior$mean_log
    }
  }
 if(ebps_method =="ind_poisson_smoothing" ){

   Mu_pm <-0*Y
   for ( i in 1:nrow(Y)){
     Mu_pm[i,]     =  log(pois_smooth_split(Y[i,],s=scaling[i],
                                        Eb_init= c(log(Y[i,]+1)) ) $Emean)
   }

  }
  if (ebps_method =="nugget" ){
    Mu_pm <- (fit_latent_nugget(Y)$Y)
  }



  #### SuSiE part ----


    tmp_Mu_pm_pen <- 0* Mu_pm
    tmp_Mu_pm_fm  <- 0* Mu_pm
    fm_pm          <- 0* Mu_pm
    init=FALSE

    #### fit EBmvFR ----
    if(fit_approach%in% c("both", "penalized")){
      tmp_Mu_pm_pen <- Mu_pm  -  fm_pm#potentially run smash on colmean


      ### TODO: Maybe use better restarting point for EBmvFR.obj
      EBmvFR.obj   <-  EBmvFR ( Y=tmp_Mu_pm_pen,
                                X              = Z,
                                tol            = tol.mrash,
                                control_mixsqp = control_mixsqp  ,
                                nullweight     = nullweight.mrash,
                                cal_obj        = cal_obj.mrash,
                                verbose        = FALSE,
                                maxit          = maxit.mrash,
                                max_step_EM =1
      )
      if(verbose){
        print( paste('Posterior of EB regression coefficient computed for iter ',iter))
      }
      b_pm <-   Z%*%  EBmvFR.obj$fitted_func

      if( fit_approach== "penalized")
        mat_mean <-   matrix(b_pm  , byrow = TRUE,
                              nrow=nrow(Z), ncol=ncol(Y))

    }else{
      b_pm <- 0* tmp_Mu_pm_pen

    }

    if(fit_approach%in% c("both", "fine_mapping")){
      tmp_Mu_pm_fm <- Mu_pm -  b_pm#potentially run smash on colmean


      susiF.obj     <- susiF (
        Y              =  tmp_Mu_pm_fm ,
        X               = X ,
        L               = L,
        tol             = tol,
        control_mixsqp  = control_mixsqp ,
        nullweight      = nullweight.mrash,
        cal_obj         = cal_obj.fsusie,
        verbose         = verbose,
        cov_lev         = cov_lev,
        min_purity      = min_purity,
        maxit           = maxit.fsusie ,
        cor_small       = cor_small,
        post_processing = post_processing)



      fm_pm <- X%*%Reduce("+",lapply(1:length(susiF.obj$cs),
                                     function(l)
                                       t(susiF.obj$fitted_func[[l]]%*% t(susiF.obj$alpha[[l]]))
      )
      )
      mat_mean <-   matrix( fm_pm , byrow = TRUE,
                            nrow=nrow(X), ncol=ncol(Y))
    }else{
      fm_pm <-0* tmp_Mu_pm_fm
      susiF.obj   <- NULL
    }




#print(    susiF.obj$cs)
  iter=iter+1




  tt_all <-exp(Mu_pm   )





  if( fit_approach ==   "both" )
  {
    susiF.obj <- fsusieR::update_cal_pip(susiF.obj)
    out <- list( Mu_pm=Mu_pm,
                 susiF.obj=susiF.obj,
                 EBmvFR.obj=EBmvFR.obj,
                 fitted = tt_all[,idx_out] )
  }

  if( fit_approach ==   "fine_mapping" )
  {
    susiF.obj <- fsusieR::update_cal_pip(susiF.obj)
    out <- list( Mu_pm=Mu_pm,
                 susiF.obj=susiF.obj,
                 fitted = tt_all[,idx_out]  )
  }
  if( fit_approach ==   "penalized")
  {
    out <- list( Mu_pm=Mu_pm,
                 EBmvFR.obj=EBmvFR.obj,
                 fitted = tt_all[,idx_out]  )
  }
  return(out)

}
