
out_prep_HF_fsusie <- function(susiF.obj  )
{

  out <- list(cs= susiF.obj$cs,
              pip= fsusieR:::update_cal_pip ( susiF.obj)$pip,
              fitted_func = susiF.obj$fitted_func )

  return(out)
}
