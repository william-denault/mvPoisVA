





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
                            maxit=100,
                            control_mixsqp=  list(verbose=FALSE,
                                                  eps = 1e-6,
                                                  numiter.em = 4
                            ),
                            thresh_lowcount=0,
                            prior_mv=  "mixture_normal_per_scale",
                            gridmult=sqrt(2),
                            nullweight.mrash=10,
                            init_pi0_w.mrash=10,
                            cov_lev=0.95,
                            min.purity     =0.5,
                            greedy=TRUE,
                            backfit=TRUE,
                            tol.mrash=1e-3,
                            verbose.mrash=TRUE,
                            maxit.mrash=10,
                            cal_obj.mrash=FALSE,
                            maxit.fsusie=50,
                            cal_obj.fsusie=FALSE
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
  }



  #### to avoid 0 in Y_min to correct at the end
  Y <- Y+1


  indx_lst <-  susiF.alpha::gen_wavelet_indx(log2(ncol(Y)))

  tl <-  lapply(1:nrow(Y), function(i)
    get_empirical_intensity(Y[i,],
                            indx_lst = indx_lst)
  )

  Y_min <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$Y_min))
  Y_tot <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$Y_tot))
  rm(tl)

  ###initialization for functional reg object -----
  lowc_wc <-   susiF.alpha::which_lowcount(Y_f=Y_min ,thresh_lowcount)

  if(verbose){
    print("done transforming data")
  }


  if(init){
    Mu_pm = logit((Y_min/Y_tot) ) #remove last column contain C coeff
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


  X <- susiF.alpha:::colScale(X)


  iter <- 1
  check <- 3*tol
  while( check >tol & iter <5){

    #### Check potential pb due to centering
    post_mat <- get_post_log_int(Mu_pm       = Mu_pm,
                                 Mu_pv       = Mu_pv,
                                 Y_min       = Y_min,
                                 Y_tot       = Y_tot,
                                 sigma2_bin  = sigma2_bin,
                                 sigma2_pois = sigma2_pois,
                                 b_pm        = b_pm,
                                 gh_points   = gh_points,
                                 tol         = tol_vga_pois)
    if(verbose){
      print( paste('Posterior log intensity computed for iter ',iter))
    }
    Mu_pm <- post_mat$A_pm
    Mu_pv <- post_mat$A_pv


    tmp_Mu_pm <- Mu_pm -  apply(Mu_pm, 2, median)#potentially run smash on colmean
    W <- list( D = tmp_Mu_pm [, -ncol(tmp_Mu_pm )],
               C = tmp_Mu_pm [,  ncol(tmp_Mu_pm )])



    if(init){


      if (fit_approach %in% c("both", "penalized")){
        temp <- susiF.alpha:: init_prior(Y              = tmp_Mu_pm,
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
        EBmvFR.obj   <-  susiF.alpha::init_EBmvFR_obj(G_prior = G_prior,
                                                      Y       = Y,
                                                      X       = X
        )
        print('Done initializing EBmvFR.obj')
      }
      if(fit_approach %in%c("both","fine_mapping")){
        temp <- susiF.alpha:: init_prior(Y              = tmp_Mu_pm,
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
        susiF.obj   <-  susiF.alpha::init_susiF_obj(L_max   = L,
                                                    G_prior = G_prior,
                                                    Y       = tmp_Mu_pm,
                                                    X       = X,
                                                    L_start = L_start,
                                                    greedy  = greedy,
                                                    backfit = backfit
        )
        print('Done initializing susiF.obj')

      }

      init=FALSE
    }
    #### fit EBmvFR ----
    if(fit_approach%in% c("both", "penalized")){
      ### TODO: Maybe use better restarting point for EBmvFR.obj
      EBmvFR.obj   <- susiF.alpha::EBmvFR.workhorse(EBmvFR.obj     = EBmvFR.obj,
                                                    W              = W,
                                                    X              = Z,
                                                    tol            = tol.mrash,
                                                    lowc_wc        = lowc_wc  ,
                                                    init_pi0_w     = init_pi0_w.mrash ,
                                                    control_mixsqp = control_mixsqp ,
                                                    indx_lst       = indx_lst,
                                                    nullweight     = nullweight.mrash,
                                                    cal_obj        = cal_obj.mrash,
                                                    verbose        = verbose.mrash,
                                                    maxit          = maxit.mrash
      )
      if(verbose){
        print( paste('Posterior of EB regression coefficient computed for iter ',iter))
      }
      b_pm <-   Z%*%  EBmvFR.obj$fitted_wc[[1]]    +  matrix( apply(Mu_pm, 2, median), byrow = TRUE,
                                                              nrow=nrow(X), ncol=ncol(Y))



    }else{
      b_pm <- matrix( apply(Mu_pm, 2, median), byrow = TRUE,
                      nrow=nrow(X), ncol=ncol(Y))
      EBmvFR.obj   <- NULL
    }

    ### fit fsusie -----
    if(fit_approach%in% c("both", "fine_mapping")){
      tmp_Mu_pm <- Mu_pm -  b_pm#potentially run smash on colmean
      W <- list( D = tmp_Mu_pm [, -ncol(tmp_Mu_pm )],
                 C = tmp_Mu_pm [,  ncol(tmp_Mu_pm )])



      susiF.obj     <- susiF.workhorse(susiF.obj      = susiF.obj,
                                       W              = W,
                                       X              = X,
                                       tol            = tol,
                                       low_wc         = low_wc,
                                       init_pi0_w     = init_pi0_w.mrash ,
                                       control_mixsqp = control_mixsqp ,
                                       indx_lst       = indx_lst,
                                       lowc_wc        = lowc_wc,
                                       nullweight     = nullweight.mrash,
                                       cal_obj        = cal_obj.fsusie,
                                       verbose        = verbose,
                                       cov_lev        = cov_lev,
                                       min.purity     = min.purity,
                                       maxit          = maxit.fsusie,
                                       tt             = temp$tt)




      fm_pm <- X%*%Reduce("+",lapply(1:length(susiF.obj$cs),
                                     function(l)
                                       sweep( sweep( susiF.obj$fitted_wc[[l]] ,
                                                     1,
                                                     1/(susiF.obj$csd_X ), "*"),1,susiF.obj$alpha[[l]])))

      susiF.obj$greedy  <- TRUE
      susiF.obj$backfit <- TRUE
    }else{
      fm_pm <-0* b_pm
      susiF.obj   <- NULL
    }


    resid <- Mu_pm-fm_pm+b_pm
    #not correct to work on later
    sigma2_pois <- var(c(resid[,ncol(resid)]))
    print(sigma2_pois)
    sigma2_bin  <- var(c(resid[,-ncol(resid)]))
    print(sigma2_bin)
    Mu_pm <- fm_pm+b_pm#update



    iter=iter+1
    ##include mr.ash
  }




  out <- out_prep_mvpoisva( Mu_pm=Mu_pm,
                            susiF.obj=susiF.obj,
                            EBmvFR.obj=EBmvFR.obj)


}
