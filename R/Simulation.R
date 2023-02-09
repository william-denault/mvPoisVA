sim_intenisty <- function( lev_res=7,
                           scaling=10,
                          length_grid= 10,
                          pi0,
                          alpha=0.8,
                          prop_decay=0.1){

  if(missing(pi0))
  {
    pi0 <- 1-exp(- (  prop_decay*(1:lev_res)))
  }
  temp <- susiF.alpha::simu_IBSS_per_level ( lev_res=lev_res,
                                             length_grid= length_grid,
                                             pi0=pi0,
                                             alpha=alpha,
                                             prop_decay=prop_decay)

  out <- list( sim_intens    = 10*abs( temp$sim_func),
               #true_coef     = abs(temp$tue_coef),
               mix_per_scale =  temp$G_level ,
               emp_pi0       =  temp$emp_pi0
  )
  return(out)
}
