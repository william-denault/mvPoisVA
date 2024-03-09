
#' @export
#'
#' @examples
#'library(mvPoisVA)
#'library(fsusieR)
#'library(susieR)
#'data(N3finemapping)
#'X <- N3finemapping$X
#'mysd=0.1
#'N =30

#'lev_res =8

#'genotype <-X[1:N,1:500]

#'idx <- which( apply( genotype,2, var ) <1e-15)
#'genotype <- genotype [, -idx]
#'count.data  <- list()
#'L <-2# sample(1:2, size =1)#actual number of effect

#'lf <-  list()
#'for(l in 1:L){
#'  lf[[l]] <-log(abs(0.2*sim_intenisty(lev_res )$sim_intens) )#functional effect for effect l
#'}


#'data(N3finemapping)
#'X <- N3finemapping$X
#'genotype <-X[sample(1:nrow(X), size=N),]

#'idx <- which( apply( genotype,2, var ) <1e-15)
#'if( length(idx)==0){
#'  X <-genotype
#'
#'  Rtrue <- cor (genotype )
#'}else{
#'  genotype <- genotype [, -idx]
#'  X <-genotype
#'
#'}
#'G<- genotype
#'X <- (X -0.99*min(X))/(0.5*max(X ))
#'
#'G <-  (G -0.99*min(G ))/(0.5*max(G ))
#'
#'tpos <- sample(1:ncol(genotype), replace = FALSE,size=2)
#'true_pos <- tpos
#'pos1 <- tpos[1]
#'pos2 <- tpos[2]
#'if( length(which(apply(G,2,var)==0))>0){
#'  G <- G[,-which(apply(G,2,var)==0)]
#'}
#'# G <- matrix( rnorm(nrow(genotype)*300), nrow = nrow(genotype))
#'
#'
#'predictor <-rep (0, length(lf[[1]] ))
#'count.data  <- list()
#'for ( i in 1:N)
#'{
#'
#'  predictor <-rep (0, length(lf[[1]] ))
#'
#'  for ( l in 1:L){
#'    predictor <-predictor + G[i, true_pos[l]]*lf[[l]]
#'  }
#'  predictor <- exp(predictor+ rnorm(  length(lf[[1]]), sd=mysd))
#'
#'  count.data [[i]] <-   rpois(n= length(lf[[1]]) ,
#'                              lambda =predictor  )
#'
#'}
#'count.data <- do.call(rbind, count.data)
#'
#'
#'Y <- count.data
#'
#'Y[1,] <- NA
#'Y[10, ] <- NA
#'out <- fit_latent_space(Y)
#'
#'image (out$Y)
#'plot(out$Y, Y)


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


  indx_lst <-  fsusieR::gen_wavelet_indx(log2(ncol(Y)))




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
  sigma2_pois=0.1

  while(  iter <20  ){#check >tol &

    init_val_pois<- c(log(Y_c+1))
    beta_pois <- c(Mu_pm)

    opt_Poisson  <- vga_pois_solver(init_val = init_val_pois ,
                                    x        = c(Y_c),
                                    s        = rep( 1, prod (dim(Y_c))),
                                    beta     = beta_pois,
                                    sigma2   = sigma2_pois,
                                    maxiter  =  50,
                                    tol      = tol,
                                    method   = 'newton')



    Mu_pm <- matrix(opt_Poisson$m,byrow = FALSE, ncol=ncol(Y_c))



    if(verbose){
      print( paste('Posterior log intensity computed for iter ',iter))
    }

    opt_Poisson$v <- opt_Poisson$v +1e-10# prevent some underflow problem
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

  colnames(out) <- colnames(Y)
  rownames(out) <- rownames(Y)

  return( list(Y=out,
               reflect=reflect )
  )
}






fit_latent_nugget<- function(Y,tol=1e-4,verbose=TRUE,reflect =FALSE,
                             est_sigma2=NULL
                             ){
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


  indx_lst <-  fsusieR::gen_wavelet_indx(log2(ncol(Y)))




  init_val_pois<- c(log(Y+1))
  beta_pois <- 0* c(log(Y+1))
  sigma2_pois=1

  ##initiatilzation for count data -----
  Y_c <- Y[complete.cases(Y),]


 # Mu_pm <-do.call(rbind, lapply(1:nrow(Y_c), function (i)pois_smooth_split(Y[i,])$Emu ))
  Mu_pm <-do.call(rbind, lapply(1:nrow(Y_c), function (i) ebps(Y[i,]
                                                                          )$posterior$mean_log))

  out <- matrix(NA, ncol=ncol(Y), nrow = nrow(Y))
  out [complete.cases(Y),] <- Mu_pm

  colnames(out) <- colnames(Y)
  rownames(out) <- rownames(Y)

  return( list(Y=out,
               reflect=reflect )
  )
}
