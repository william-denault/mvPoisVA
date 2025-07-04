


#' Implementation of the Poisson fSuSiE
#' # using the Poisson multiscale approach use in multiseq by Shim et al 2024


mv_Poisproc_reg <- function(Y,
                            Z,
                            X,
                            L=3,
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
                            thresh_lowcount=1e-2,
                            prior_mv=  "mixture_normal_per_scale",
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
                            maxit.fsusie=50,
                            cal_obj.fsusie=FALSE,
                            max_SNP_EM     = 100,
                            max_step_EM    = 1,
                            cor_small=TRUE,
                            max.iter=3
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

  ##initiatilzation -----
  init=TRUE
  gh_points = fastGHQuad::gaussHermiteData(n_gh)
  J = log2(ncol(Y)); if((J%%1) != 0) reflect=TRUE
  if(reflect){
    tl <- lapply(1:nrow(Y), function(i) reflect_vec(Y[i,]))
    Y <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$x))
    idx_out <- tl[[1]]$idx #### indx of interest at the end
  }else{
    idx_out <- 1: ncol(Y)
  }
  #### to avoid 0 in Y_min to correct at the end
  Y <- Y+1

  if(fit_approach %in% c('both',"fine_mapping")){
    tidx <- which(apply(X,2,var)==0)
    if( length(tidx)>0){
      warning(paste("Some of the columns of X are constants, we removed" ,length(tidx), "columns"))
      X <- X[,-tidx]
    }
    X <- fsusieR:::colScale(X)
  }
  if(fit_approach %in% c('both',"penalized")){
    tidx <- which(apply(Z,2,var)==0)
    if( length(tidx)>0){
      warning(paste("Some of the columns of Z are constants, we removed" ,length(tidx), "columns"))
      Z <- Z[,-tidx]
    }
    Z <- fsusieR:::colScale(Z)
  }

  indx_lst <-  fsusieR::gen_wavelet_indx(log2(ncol(Y)))

  tl <-  lapply(1:nrow(Y), function(i)
    get_empirical_intensity(Y[i,],
                            indx_lst = indx_lst)
  )

  Y_min <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$Y_min))
  Y_tot <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$Y_tot))
  rm(tl)

  if( length(which(apply(Y_min,2, var )<= thresh_lowcount))==0)
  {
    lowc_wc <-NULL
  }else{
    lowc_wc <- which(apply(Y_min,2, var )<= thresh_lowcount)
  }





  if(verbose){
    print("done transforming data")
  }


  if(init){
    Mu_pm = logit((Y_min/Y_tot) ) #remove last column contain C coeff
    #Mu_pm[which(is.na(Mu_pm))]<-Inf
    Mu_pm[Mu_pm==-Inf] =  logit(0.1)
    Mu_pm[Mu_pm==Inf]  =  logit(0.9)

    Mu_pm[,ncol(Y_min)] <- log(Y_min[,ncol(Y_min)])
    Mu_pv = 1/Y_tot
    if(missing(init_b_pm)){
      b_pm <- 0*Y
    }
    sigma2_bin  = 1
    sigma2_pois = 1

  }





  iter <- 1
  check <- 3*tol
  Mu_0 <- Mu_pm
  while( check >tol & iter <max.iter ){

    #### Check potential pb due to centering
    post_mat <- get_post_log_int(Mu_pm       = Mu_pm,
                                 Mu_pv       = Mu_pv,
                                 Y_min       = Y_min,
                                 Y_tot       = Y_tot,
                                 sigma2_bin  = 1/sigma2_bin,
                                 sigma2_pois =   1/sigma2_pois,
                                 b_pm        = Mu_pm,
                                 gh_points   = gh_points,
                                 tol         = tol_vga_pois)
    if(verbose){
      print( paste('Posterior log intensity computed for iter ',iter))
    }

    Mu_pm <- post_mat$A_pm

    Mu_pv <- post_mat$A_pv


    b_pm <- 0* Mu_pm
    fm_pm <- 0* Mu_pm

    if(init){

      tmp_Mu_pm <- fsusieR::colScale(Mu_pm, scale = FALSE)#potentially run smash on colmean

      lowc_wc <-  which_lowcount(tmp_Mu_pm,
                                 thresh_lowcount=thresh_lowcount)
      W <- list( D = tmp_Mu_pm [, -ncol(tmp_Mu_pm )],
                 C = tmp_Mu_pm [,  ncol(tmp_Mu_pm )])
      if (fit_approach %in% c("both", "penalized")){
        temp <- fsusieR:: init_prior(Y              = tmp_Mu_pm,
                                         X              = Z ,
                                         prior          = prior_mv ,
                                         v1             = v1,
                                         indx_lst       = indx_lst,
                                         lowc_wc        = lowc_wc,
                                         control_mixsqp = control_mixsqp,
                                         nullweight     = nullweight.mrash,
                                         gridmult       = gridmult )
        G_prior     <- temp$G_prior


        #Recycled for the first step of the while loop
        EBmvFR.obj   <-  fsusieR::init_EBmvFR_obj(G_prior = G_prior,
                                                      Y       = Y,
                                                      X       = Z
        )
        print('Done initializing EBmvFR.obj')
      }
      if(fit_approach %in%c("both","fine_mapping")){
        temp <- fsusieR:: init_prior(Y              = tmp_Mu_pm,
                                         X              = X ,
                                         prior          = prior_mv ,
                                         v1             = v1,
                                         indx_lst       = indx_lst,
                                         lowc_wc        = lowc_wc,
                                         control_mixsqp = control_mixsqp,
                                         nullweight     = nullweight.mrash,
                                         gridmult       = gridmult )
        G_prior     <- temp$G_prior


        #Recycled for the first step of the while loop
        susiF.obj   <-  fsusieR::init_susiF_obj(L_max   = L,
                                                    G_prior = G_prior,
                                                    Y       = tmp_Mu_pm,
                                                    X       = X,
                                                    L_start = L_start,
                                                    greedy  = greedy,
                                                    backfit = backfit
        )
        print('Done initializing susiF.obj')

      }
      tmp_Mu_pm_pen <- 0*tmp_Mu_pm
      tmp_Mu_pm_fm  <- 0*tmp_Mu_pm
      init=FALSE
    }
    #### fit EBmvFR ----
    if(fit_approach%in% c("both", "penalized")){
      tmp_Mu_pm_pen <- Mu_pm  -  fm_pm#potentially run smash on colmean

      t_mean_EBmvFR <-  apply(tmp_Mu_pm_pen,2, mean )
      tmp_Mu_pm_pen <- fsusieR::colScale(tmp_Mu_pm_pen, scale=FALSE)
      W <- list( D = tmp_Mu_pm [, -ncol(tmp_Mu_pm_pen )],
                 C = tmp_Mu_pm [,  ncol(tmp_Mu_pm_pen )])


      ### TODO: Maybe use better restarting point for EBmvFR.obj
      EBmvFR.obj   <- fsusieR::EBmvFR.workhorse( obj     = EBmvFR.obj,
                                                    W              = W,
                                                    X              = Z,
                                                    tol            = tol.mrash,
                                                    lowc_wc        = lowc_wc  ,
                                                    init_pi0_w     = init_pi0_w.mrash ,
                                                    control_mixsqp = control_mixsqp ,
                                                    indx_lst       = indx_lst,
                                                    nullweight     = nullweight.mrash,
                                                    cal_obj        = cal_obj.mrash,
                                                    verbose        = FALSE,
                                                    maxit          = maxit.mrash,
                                                    max_step_EM =1
      )
      if(verbose){
        print( paste('Posterior of EB regression coefficient computed for iter ',iter))
      }
      b_pm <-   Z%*%  EBmvFR.obj$fitted_wc[[1]]

      if( fit_approach== "penalized")
        mat_mean <-   matrix( t_mean_EBmvFR , byrow = TRUE,
                              nrow=nrow(X), ncol=ncol(Y))

    }else{
      b_pm <- 0* tmp_Mu_pm_pen

    }

    if(fit_approach%in% c("both", "fine_mapping")){
      tmp_Mu_pm_fm <- Mu_pm -  b_pm#potentially run smash on colmean

      t_mean_susiF <-  apply(tmp_Mu_pm_fm,2, mean )

      tmp_Mu_pm_fm <- fsusieR::colScale(tmp_Mu_pm_fm, scale=FALSE)
      W <- list( D = tmp_Mu_pm [, -ncol(tmp_Mu_pm_fm )],
                 C = tmp_Mu_pm [,  ncol(tmp_Mu_pm_fm )])

      susiF.obj     <- fsusieR:::susiF.workhorse( obj      = susiF.obj,
                                       W              = W,
                                       X              = X,
                                       tol            = tol,
                                       init_pi0_w     = init_pi0_w.mrash ,
                                       control_mixsqp = control_mixsqp ,
                                       indx_lst       = indx_lst,
                                       lowc_wc        = lowc_wc,
                                       nullweight     = nullweight.mrash,
                                       cal_obj        = cal_obj.fsusie,
                                       verbose        = verbose,
                                       cov_lev        = cov_lev,
                                       min_purity     = min_purity,
                                       maxit          = maxit.fsusie,
                                       tt             = temp$tt,
                                       cor_small      =cor_small)





      fm_pm <- X%*%Reduce("+",lapply(1:length(susiF.obj$cs),
                                     function(l)
                                       sweep(  susiF.obj$fitted_wc[[l]]  ,
                                               1,
                                               susiF.obj$alpha[[l]],
                                               "*")
      )
      )
      mat_mean <-   matrix( t_mean_susiF , byrow = TRUE,
                            nrow=nrow(X), ncol=ncol(Y))
    }else{
      fm_pm <-0* tmp_Mu_pm_fm
      susiF.obj   <- NULL
    }



    resid <- Mu_pm -mat_mean -fm_pm-b_pm
    #not correct to work on later
    sigma2_pois <- var(c(resid[,ncol(resid)]))
    #print(sigma2_pois)
    sigma2_bin  <- var(c(resid[,-ncol(resid)]))
    #print(sigma2_bin)
    Mu_pm <- mat_mean +fm_pm+b_pm#update

    # print(    susiF.obj$cs)
    if(verbose){print( paste(   iter, " done"))}
    iter=iter+1
    ##include mr.ash

    # par (mfrow=c(1,2))
    tt <- reverse_intensity_transform (vec_int = c(log( sigmoid(  Mu_pm[1,-ncol( Mu_pm)])),
                                                   log(Y_min[1,ncol(Y_min)])  ),
                                       indx_lst = indx_lst,
                                       is.logprob = TRUE,
                                       is.log_int = TRUE)
    #plot ( Y[1,], col="blue")
    # points (tt)
    # lines(tt, col="green")


    #plot( Y[1,],tt)

    # abline(a=0,b=1)
    # par (mfrow=c(1,1))


  }



  tt_all <- do.call(rbind,lapply( 1: nrow(Y), function(i)  reverse_intensity_transform (vec_int = c(log( sigmoid(  Mu_pm[i,-ncol( Mu_pm)])),
                                                                                                    log(Y_min[i,ncol(Y_min)])  ),
                                                                                        indx_lst = indx_lst,
                                                                                        is.logprob = TRUE,
                                                                                        is.log_int = TRUE))
  )

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




out_prep_mvpoisva <- function( Mu_pm=Mu_pm,
                               susiF.obj=susiF.obj,
                               EBmvFR.obj=EBmvFR.obj)
{
  out <- list(cs= susiF.obj$cs,
              pip= fsusieR:::update_cal_pip( susiF.obj)$pip)

  return(out)
}
