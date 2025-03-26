#' Function to get fitted B-spline coefficients for smoothed functional covariates
#'
#' @param func_data_list List of functional covariates. Each list element should be a dataframe with rows indicating time points and columns for each subject. Importantly, each element of the list of functional covariates must be named.
#' @param time_points Vector of observed time points. Method assumes all subjects and functional covariates are observed at a common set of time points.
#' @param break_spacing A value in c("quantiles", "equal", "manual") to set break spacing for B-spline basis.
#' @param n_breaks A number indicating number of breaks to use for B-spline basis (includes endpoints, results in n_breaks + 2 B-spline polynomials).
#' @param lambda_vals A vector of positive numbers to cross-validate over to get roughness penalty.
#' @param min_method A value to minimize during cross-validation, either "mean" or "quantile" (over the subjects GCV values)
#' @param min_quantile if min_method = "quantile", a value for quantile to use
#' @param num_threads Number of threads to smooth functional covariates similtaneously
#' @param manual_breaks An optional vector of manual break values.
#'
#' @returns A list of four objects. Object (1) is a data frame of fitted B-spline coefficients for all subjects. Object (2) is a list of the fdSmooth objects. Object (3) is a vector of the optimal lambdas. Object (4) is a basis object of the type bspline.
#' @export
#'
#' @examples
#'
#' library(fds)
#' library(fda)
#' library(dplyr)
#'
#' data("Fatspectrum")
#' data_list = list(fat = Fatspectrum$y) # list elements must be named
#' observed_wavelengths = Fatspectrum$x
#'
#' fit_bspline = get_bspline_coefs(func_data_list = data_list, time_points = observed_wavelengths,
#'                                 n_breaks = 10, lambda_vals = 10^seq(-6, 6, length.out = 250),
#'                                 min_method = "mean", num_threads = 1)
#'
#'
#'
get_bspline_coefs = function(func_data_list, time_points, break_spacing = "quantiles",
                             n_breaks = 10, lambda_vals = 10^seq(-6, 6, length.out = 100),
                             min_method = "mean", min_quantile,
                             num_threads = parallel::detectCores(),
                             manual_breaks) {
  # This function fits a smooth function to noisy functional data using WLS with a roughness
  # penalty. We assume that subjects are observed at the same time points
  ##### Return 1: a data frame of b-spline coefficients for all functional covariates
  ##### Return 2: a list of the smoothed covariate fd objects
  ##### Return 3: the optimal lambda for penalizing roughness
  ##### Return 4: the basis object used
  # ARGS:
  # `func_data_list` = list of functional data frame (1 for each covariate), each data frame should have a column
  # for each subject and rows representing the values over time
  # NOTE: DATA LIST MUST HAVE NAMED ENTRIES
  # `time_points` = a vector of time points i.e. 1:nrow(data[[1]])
  # `break_spacing` = a value in c("quantiles", "equal", "manual") to set break spacing
  # `n_breaks` = a number indicating number of breaks to use
  # `lambda_vals` = a vector of positive numbers to cross-validate over to get roughness penalty
  # `min_method` = value to minimize during cross-validation, either "mean" or "quantile
  # `min_quantile` = if min_method = "quantile", a value for quantile to use
  # `num_threads` = number of threads to use for processing
  # `manual_breaks` = if break_spacing = "manual", a vector of breaks to use, must include first and last time point


  # specify evenly spaced breaks or breaks at evenly spaced quantiles
  if (break_spacing == "quantiles") {
    quantiles = seq(from = 0, to = 1, length.out = n_breaks) # quantiles by number breaks
    quantile_ind = floor(quantiles * length(time_points)) # indices of breaks ()

    breaks = time_points[quantile_ind]
    breaks = sort(unique(c(min(time_points), breaks, max(time_points)))) # adds first and last point
  } else if (break_spacing == "equal") {
    breaks = seq(from = min(time_points), to = max(time_points), length.out = n_breaks)
  } else if (break_spacing == "manual") {
    breaks = manual_breaks
  }

  # create basis
  basis = fda::create.bspline.basis(rangeval = range(time_points),
                                    breaks = breaks, norder = 4)

  # lapply() over covariates using the GCV function
  smooth_covars_list = parallel::mclapply(names(func_data_list), function(name) {
    # get one covariate
    covar_df = func_data_list[[name]]

    # get GCV fd object and return with optimal lambda
    if (min_method == "mean") {
      smooth_opt = smooth_w_gcv(covar_df = covar_df, basis = basis,
                                time_points = time_points, lambda_vals = lambda_vals,
                                min_method = min_method)
    } else if (min_method == "quantile") {
      smooth_opt = smooth_w_gcv(covar_df = covar_df, basis = basis,
                                time_points = time_points, lambda_vals = lambda_vals,
                                min_method = min_method, quantile = min_quantile)
    }

    fd_opt = smooth_opt[[1]]
    lambda_opt = smooth_opt[[2]]
    # fix names and return
    names(fd_opt$fdnames) = c("time", "reps", name)
    list(smooth_fd = fd_opt, lambda = lambda_opt)
  }, mc.cores = num_threads)

  # fix names
  names(smooth_covars_list) = names(func_data_list)

  # lapply() to get coefficients for splines from smooth_covars_list
  coef_df_list = parallel::mclapply(names(smooth_covars_list), function(name) {
    fd_obj = smooth_covars_list[[name]][["smooth_fd"]]
    coefs = t(fd_obj$coefs)
    colnames(coefs) = paste0(name, "_bspl_", seq_len(ncol(coefs)))
    as.data.frame(coefs)
  }, mc.cores = num_threads)

  # get vector of optimal lambdas
  lambdas = sapply(smooth_covars_list, function(x) x$lambda)

  # combine all covariates into one data frame
  to_return = list(coef_df = do.call(cbind, coef_df_list),
                   smooth_covars_list = smooth_covars_list,
                   lambda_opt = lambdas,
                   basis = basis)

  return(to_return)
}
