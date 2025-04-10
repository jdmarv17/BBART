#' Functional Bayesian Additive Regression Trees (fBART)
#'
#' @param y A vector of outcome values.
#' @param x A dataframe of scalar covariates. If 'nxp' then 'y' must have length 'n' and each dataframe in 'func_data_list' must have 'n' columns.
#' @param func_data_list List of functional covariates. Each list element should be a dataframe with rows indicating time points and columns for each subject. Importantly, each element of the list of functional covariates must be named.
#' @param num_trees Number of BART trees per sample.
#' @param num_samps Number of posterior samples to take after burn-in.
#' @param num_burn Number of burn-in samples.
#' @param num_thin Keep every 'num_thin' draw.
#' @param num_chains Number of independent MCMC chains to run.
#' @param num_threads_bart Number of threads used to run BART - not recommended to be greater than 'num_chains'.
#' @param num_threads_wrangle Number of threads for post processing - can be larger than 'num_chains'.
#' @param prior_power Power parameter for tree prior.
#' @param prior_base Base parameter for tree prior.
#' @param time_points Vector of observed time points. Method assumes all subjects and functional covariates are observed at a common set of time points.
#' @param break_spacing A value in c("quantiles", "equal", "manual") to set break spacing for B-spline basis.
#' @param n_breaks A number indicating number of breaks to use for B-spline basis (includes endpoints, results in n_breaks + 2 B-spline polynomials).
#' @param lambda_vals A vector of positive numbers to cross-validate over to get roughness penalty.
#' @param min_method A value to minimize during cross-validation, either "mean" or "quantile" (over the subjects GCV values)
#' @param min_quantile if min_method = "quantile", a value for quantile to use
#' @param manual_breaks An optional vector of manual break values. Must include domain endpoints.
#'
#' @returns A list of three objects. Object (1) is an indicator data frame of variable tree inclusion, object (2) is a data frame of posterior trees without terminal nodes, and object (3) is the dbarts::bart2 object.
#' @export
#'
#' @examples
#' library(fds)
#' library(fda)
#' library(dplyr)

#' data("Fatspectrum")
#' data_list = list(fat = Fatspectrum$y) # list elements must be named
#' observed_wavelengths = Fatspectrum$x
#' y = unname(Fatvalues)

#' fbart_run = fbart(y = y, func_data_list = data_list, time_points = observed_wavelengths,
#'                   num_trees = 100, num_samps = 5000, num_burn = 5000, num_thin = 5, num_chains = 1,
#'                   num_threads_bart = 1, num_threads_wrangle = 1, break_spacing = "quantiles",
#'                   n_breaks = 10, lambda_vals = 10^seq(-6, 6, length.out = 250), min_method = "mean")
#'
#'
fbart = function(y, x = NULL, func_data_list,
                 num_trees = 100, num_samps = 5000, num_burn = 5000,
                 num_thin = 5, num_chains = min(4, dbarts::guessNumCores()),
                 num_threads_bart = min(num_chains, dbarts::guessNumCores()),
                 num_threads_wrangle = dbarts::guessNumCores(),
                 prior_power = 4, prior_base = 0.95,
                 time_points, break_spacing = "quantiles",
                 n_breaks = 8, lambda_vals = 10^seq(-6, 6, length.out = 100),
                 min_method = "mean", manual_breaks = NULL, min_quantile = NULL) {

  # get b spline coefs
  bspline = get_bspline_coefs(func_data_list = func_data_list, time_points = time_points,
                              break_spacing = break_spacing, n_breaks = n_breaks,
                              manual_breaks = manual_breaks,
                              lambda_vals = lambda_vals, min_method = min_method,
                              num_threads = num_threads_wrangle,
                              min_quantile = min_quantile)

  # if available bind scalar covariates with outcome and bspline coefficients
  if (!is.null(x)) {
    data_df = cbind(y = y, x, bspline$coef_df)
  } else {
    data_df = cbind(y = y, bspline$coef_df)
  }

  # bart formula
  formula = as.formula(
    paste("y ~", paste(paste0("`", colnames(data_df[, !names(data_df) %in% "y"]), "`"), collapse = " + ")))

  # run bart
  bart_run = run_bart(formula, data = data_df,
                      num_trees = num_trees, num_samps = num_samps,
                      num_burn = num_burn, num_thin = num_thin,
                      num_chains = num_chains, num_threads_bart = num_threads_bart,
                      num_threads_wrangle = num_threads_wrangle,
                      prior_power = prior_power, prior_base = prior_base)

  to_return = list(indicators = bart_run[[1]], splits = bart_run[[2]], bart_model = bart_run[[3]],
                   basis = bspline[[4]], lambdas = bspline[[3]])

  return(to_return)
}
