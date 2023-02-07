mv_Poisproc_reg <- function(Y,
                            X,
                            reflect =FALSE,
                            verbose=TRUE,
                            n_gh = 10,
                            init_b_pm,
                            tol= 1e-3,
                            tol_vga_pois=1e-5,
                            maxit=100)
{

  ##initiatilzation -----
  init=TRUE
  gh_points = fastGHQuad::gaussHermiteData(n_gh)
  J = log2(ncol(Y)); if((J%%1) != 0) reflect=TRUE
  if(reflect){
    tl <- lapply(1:nrow(Y), function(i) reflect_vec(Y[i,]))
    Y <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$x))
    idx_out <- tl[[1]]$idx #### indx of interest at the end
  }



  indx_lst <-  susiF.alpha::gen_wavelet_indx(log2(ncol(Y)))

  tl <-  lapply(1:nrow(Y), function(i)
    get_empirical_intensity(Y[i,],
                            indx_lst = indx_lst)
  )

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

  while( check >tol & iter <maxit){


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



    ##include fsusie
  }





}
