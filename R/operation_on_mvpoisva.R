
out_prep_mvpoisva <- function( Mu_pm=Mu_pm,
                          susiF.obj=susiF.obj,
                          EBmvFR.obj=EBmvFR.obj)
{
 out <- list(cs= susiF.obj$cs,
             pip= fsusieR:::update_cal_pip( susiF.obj)$pip)

 return(out)
}
