

#out_prep_HF_fsusie <- function(susiF.obj,Y, X, indx_lst, filter.cs, lfsr_curve, outing_grid,...)
  #{

  # susiF.obj <-  fsusieR::update_cal_pip(susiF.obj)
  # susiF.obj <-  update_cal_fit_func_HF_fsusie(susiF.obj, indx_lst)
  # susiF.obj <-  update_cal_credible_band_HF_fsusie(susiF.obj, indx_lst)
  # susiF.obj <-  fsusieR::name_cs(susiF.obj,X)
  # if(filter.cs)
    # {
    #   susiF.obj <- fsusieR::check_cs(susiF.obj,min.purity=0.5,X=X)
    # }
  # susiF.obj <-  update_cal_indf_HF_fsusie(susiF.obj, Y, X, indx_lst)

  #  susiF.obj$outing_grid <- outing_grid
  #  susiF.obj$purity      <-   fsusieR::cal_purity(l_cs= susiF.obj$cs, X=X)
  #  return(susiF.obj)
#}



#To define
#update_cal_credible_band_HF_fsusie
#update_cal_fit_func_HF_fsusie
#update_cal_indf_HF_fsusie

out_prep_HF_fsusie <- function(susiF.obj  )
{

  out <- list(cs= susiF.obj$cs,
              pip= fsusieR:::update_cal_pip ( susiF.obj)$pip,
              fitted_func = susiF.obj$fitted_func )

  return(out)
}
