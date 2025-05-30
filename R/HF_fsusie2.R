#' @title Sum of Single Function (SuSiF) with Haar Fizs transform
#'
#' @description This function implements the SuSiF method for analyzing functional phenotypes using wavelet-based techniques using Haar Fisz transform for Poisson data.
#'
#' @details The SuSiF method is designed to work with functional phenotypes, where the underlying algorithm employs wavelets. This function accommodates various parameters for fine-tuning the analysis, such as specifying priors, convergence tolerance, and effect-fitting bounds.
#'
#' @param Y A matrix of size N by J representing the functional phenotype. The algorithm assumes that J is a power of 2. If J is not a power of 2, SuSiF internally remaps the data to a grid of length 2^J.
#' @param Z A matrix of covariates. Details were not provided in the original documentation, please add relevant details.
#' @param X A matrix of size N by P containing the covariates.
#' @param L An upper bound on the number of effects to fit. Defaults to 2 if not specified.
#' @param nugget Logical, if `TRUE`, the nugget effect is included in the model. Defaults to `TRUE`.
#' @param L_start The number of effects initialized at the start of the algorithm. Defaults to 3.
#' @param reflect Logical, if `TRUE`, reflection of the data is performed. Defaults to `FALSE`.
#' @param verbose Logical, if `TRUE`, the function prints progress information. Defaults to `TRUE`.
#' @param n_gh Number of Gauss-Hermite quadrature points. Defaults to 10.
#' @param init_b_pm Initial values for the posterior mean of the coefficients. This parameter needs to be documented.
#' @param tol A small, non-negative number specifying the convergence tolerance for the fitting procedure. Defaults to 1e-3.
#' @param tol_vga_pois Convergence tolerance for the variational Gaussian approximation (VGA) for Poisson distribution. Defaults to 1e-5.
#' @param maxit Maximum number of iterations for the fitting procedure. Defaults to 10.
#' @param control_mixsqp A list of control parameters for the `mixsqp` function. Defaults are provided.
#' @param thresh_lowcount A numeric threshold for handling wavelet coefficients with problematic distributions. Defaults to 0.
#' @param gridmult A numeric value used to control the number of components in the mixture prior. Defaults to `sqrt(2)`.
#' @param nullweight.mrash A numeric value for penalizing likelihood at point mass 0 in the `mrash` method. Defaults to 10.
#' @param init_pi0_w.mrash Initial weight on the null component in the `mrash` method. Defaults to 10.
#' @param cov_lev A numeric value between 0 and 1 for the expected level of coverage of credible sets. Defaults to 0.95.
#' @param min_purity Minimum purity for estimated credible sets. Defaults to 0.5.
#' @param greedy Logical, if `TRUE`, allows for a greedy search for extra effects. Defaults to `TRUE`.
#' @param backfit Logical, if `TRUE`, allows discarding effects via backfitting. Defaults to `TRUE`.
#' @param tol.mrash Convergence tolerance for the `mrash` method. Defaults to 1e-3.
#' @param verbose.mrash Logical, if `TRUE`, prints progress information for the `mrash` method. Defaults to `TRUE`.
#' @param maxit.mrash Maximum number of iterations for the `mrash` method. Defaults to 10.
#' @param cal_obj.mrash Logical, if `TRUE`, computes ELBO for `mrash`. Defaults to `FALSE`.
#' @param maxit.fsusie Maximum number of iterations for the `fsusie` method. Defaults to 50.
#' @param cal_obj.fsusie Logical, if `TRUE`, computes ELBO for `fsusie`. Defaults to `FALSE`.
#' @param max_SNP_EM Maximum number of SNPs considered in the EM algorithm. Defaults to 100.
#' @param max_step_EM Maximum number of steps in the EM algorithm. Defaults to 1.
#' @param cor_small Logical, if `TRUE`, corrects for small values in the correlation matrix. Defaults to `TRUE`.
#' @param post_processing Specifies the post-processing method. Defaults to "HMM".
#'
#' @return The function returns a list containing the SuSiF results and related components.
#' @importFrom stats var
#' @export
HF_susiF2 <- function(Y,
                     Z,
                     X,
                     L=3,
                     nugget=TRUE,
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
                     thresh_lowcount=0,

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
                     post_processing="HMM"
)
{


  Y_0 <-  Y
  Y    <- HFT(Y)

  susiF.obj     <- fsusieR::susiF (
    Y               =  Y ,
    X               = X,
    L               = L,
    tol             = tol,
    control_mixsqp  = control_mixsqp ,
    nullweight      = nullweight.mrash,
    cal_obj         = cal_obj.fsusie,
    verbose         = verbose,
    cov_lev         = cov_lev,
    min_purity      = min_purity,
    maxit           = maxit.fsusie ,
    cor_small       = cor_small,
    post_processing =post_processing)




  #preparing output
  # susiF.obj <- fsusieR::update_cal_pip(susiF.obj)
  out <-  susiF.obj
  return(  out)
}
