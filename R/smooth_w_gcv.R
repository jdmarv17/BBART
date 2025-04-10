#' Function to get weighted least squares fit using B-splines
#'
#' @param covar_df A functional covariate data frame consisting of rows for time points and columns for subjects.
#' @param basis A basis object of the type bspline.
#' @param time_points Vector of observed time points. Method assumes all subjects and functional covariates are observed at a common set of time points.
#' @param lambda_vals A vector of positive numbers to cross-validate over to get roughness penalty.
#' @param min_method A value to minimize during cross-validation, either "mean" or "quantile" (over the subjects GCV values).
#' @param quantile if min_method = "quantile", a value for quantile to use
#'
#' @returns A list of two objects. Object (1) is an fdSmooth object using object (2) lambda_opt to smooth.
#' @export
#'
#'
#'
smooth_w_gcv = function(covar_df, basis, time_points,
                        lambda_vals, min_method = "mean", quantile = NULL) {



  if (min_method == "mean") {
    gcv = sapply(lambda_vals, function(lambda) {
      fd = fda::fdPar(basis, Lfdobj = int2Lfd(2), lambda = lambda)
      smooth_covar = fda::smooth.basis(argvals = time_points, y = as.matrix(covar_df), fdParobj = fd)
      mean(smooth_covar$gcv)
    })
  } else if (min_method == "quantile") {
    if (is.null(quantile)) {
      stop("You must provide a non-null quantile value when min_method = 'quantile'")
    }

    gcv = sapply(lambda_vals, function(lambda) {
      fd = fda::fdPar(basis, Lfdobj = int2Lfd(2), lambda = lambda)
      smooth_covar = fda::smooth.basis(argvals = time_points, y = as.matrix(covar_df), fdParobj = fd)
      stats::quantile(smooth_covar$gcv, probs = quantile)
    })
  }

  # get lambda with min gcv
  lambda_opt = lambda_vals[which.min(gcv)]
  # get final smoothed covariate
  fd_opt = fda::fdPar(basis, Lfdobj = int2Lfd(2), lambda = lambda_opt)
  smoothed_fd = fda::smooth.basis(argvals = time_points, y = as.matrix(covar_df), fdParobj = fd_opt)$fd

  to_return = list(smoothed_fd, lambda_opt)

  return(to_return)
}
