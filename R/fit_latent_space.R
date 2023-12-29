fit_latent_space <- function(Y,tol=1e-4,verbose=TRUE,reflect =FALSE){
  ##initiatilzation -----

  J = log2(ncol(Y)); if((J%%1) != 0) reflect=TRUE
  if(reflect){
    tl <- lapply(1:nrow(Y), function(i) reflect_vec(Y[i,]))
    Y <- do.call(rbind, lapply(1:length(tl), function(i) tl[[i]]$x))
    idx_out <- tl[[1]]$idx #### indx of interest at the end
  }else{
    idx_out <- 1: ncol(Y)
  }
  #### to avoid 0 in Y_min to correct at the end


  indx_lst <-  susiF.alpha::gen_wavelet_indx(log2(ncol(Y)))




  init_val_pois<- c(log(Y+1))
  beta_pois <- 0* c(log(Y+1))
  sigma2_pois=1

  ##initiatilzation for count data -----
  Y_c <- Y[complete.cases(Y),]
  Mu_pm<-Y_c
  iter=1
  beta_pois <- 0* c(log(Mu_pm +1))
  check <- 3*tol

  ##### Poisson Part ----


  while(  iter <15  ){#check >tol &

    init_val_pois<- c(log(Y_c+1))
    beta_pois <- c(Mu_pm)
    sigma2_pois=1
    opt_Poisson  <- vga_pois_solver(init_val = init_val_pois ,
                                    x        = c(Y_c),
                                    s        = rep( 1, prod (dim(Y_c))),
                                    beta     = beta_pois,
                                    sigma2   = sigma2_pois,
                                    maxiter  = 10,
                                    tol      = tol,
                                    method   = 'newton')



    Mu_pm <- matrix(opt_Poisson$m,byrow = FALSE, ncol=ncol(Y_c))



    if(verbose){
      print( paste('Posterior log intensity computed for iter ',iter))
    }


    tt <-  ashr::ash(opt_Poisson$m,opt_Poisson$v)

    resid <- Mu_pm -matrix( tt$result$PosteriorMean,byrow = FALSE, ncol=ncol(Y_c))
    #not correct to work on later
    sigma2_pois <- var(c(resid ))
    #print(sigma2_pois)
    Mu_pm <- matrix( tt$result$PosteriorMean,byrow = FALSE, ncol=ncol(Y_c))
    iter=iter+1

  }


  out <- matrix(NA, ncol=ncol(Y), nrow = nrow(Y))
  out [complete.cases(Y),] <- Mu_pm
  return( list(Y=out,
          reflect=reflect )
  )
}
