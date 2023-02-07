mv_Poisproc_reg <- function(Y,
                            X,
                            reflect =FALSE,
                            verbose=TRUE,
                            n_gh = 10,
                            init_b_pm,
                            tol= 1e-3,
                            tol_vga_pois=1e-5,
                            maxit=100,
                            maxit.mrash=100,
                            tol.mrash=1e-3,
                            cal_obj.mrash=FALSE,
                            verbose.mrash=FALSE,
                            nullweight.mrash,
                            init_pi0_w.mrash,
                            control_mixsqp=  list(verbose=FALSE,
                                                  eps = 1e-6,
                                                  numiter.em = 4
                                                  ),
                            thresh_lowcount=0
                            )
{

  ## static param  ----

  if(missing(nullweight.mrash)){
    nullweight.mrash=10
  }
  if(missing(init_pi0_w.mrash)){
    init_pi0_w.mrash=10
  }
  init=TRUE
  gh_points = fastGHQuad::gaussHermiteData(n_gh)

##initiatilzation for count data -----



  ### dealing with non 2^S data ----
  J = log2(ncol(Y)); if((J%%1) != 0) reflect=TRUE
  if(reflect){
    tl <- lapply(1:nrow(Y), function(i) reflect_vec(Y[i,]))
    Y <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$x))
    idx_out <- tl[[1]]$idx #### indx of interest at the end
  }
 indx_lst <-  susiF.alpha::gen_wavelet_indx(log2(ncol(Y)))
 ### Wavelet like transform -----
  tl <-  lapply(1:nrow(Y), function(i)
    get_empirical_intensity(Y[i,],
                            indx_lst = indx_lst)
  )
 ### Cal Ymin Ytot -----
  Y_min <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$Y_min))
  Y_tot <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$Y_tot))
  rm(tl)

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
  ##initialization for functional reg object -----
  lowc_wc <-   susiF.alpha::which_lowcount(Y_f,thresh_lowcount)


 while( check >tol & iter <maxit){

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
   Mu_pm <- post_mat$A_pm
   Mu_pv <- post_mat$A_pv


   tmp_Mu_pm <- Mu_pm - colMeans(Mu_pm, na.rm = TRUE) #potentially run smash on colmean


   EBmvFR.obj   <- EBmvFR.workhorse(EBmvFR.obj     = EBmvFR.obj,
                                    W              = tmp_Mu_pm ,
                                    X              = X,
                                    tol            = tol.mrash,
                                    low_wc         = low_wc,
                                    init_pi0_w     = init_pi0_w.mrash ,
                                    control_mixsqp = control_mixsqp ,
                                    indx_lst       = indx_lst,
                                    nullweight     = nullweight.mrash,
                                    cal_obj        = cal_obj.mrash,
                                    verbose        = verbose.mrash,
                                    maxit          = maxit.mrash)

   b_pm <-   X%*%EBmvFR.obj$fitted_func+ colMeans(Mu_pm, na.rm = TRUE)
   ##include mr.ash
 }

### !!!! idx_out



}
