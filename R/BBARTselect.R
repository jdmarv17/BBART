#' BBART variable selection
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
#' @param alpha_g Determines the 1 - alpha_g quantile for a global threshold
#' @param alpha_l Determines the 1 - alpha_l quantile for separate variable thresholds
#' @param alpha_c Determines the 1 - alpha_c quantile for the global multiplier method
#'
#' @returns A list of returns including any the appended data matrix, the TIPs from the real BBART run,
#' the basis used for the functional covariate basis expansions, any calculated thresholds, and any selected
#' scalar and functional covariates. Important domain regions can be identified by examining specific basis component
#' selections for functional covariates.
#' @export
#'
#' @examples
fBARTselect = function(y, x = NULL, func_data_list,
                       num_trees = 100, num_samps = 5000, num_burn = 5000,
                       num_thin = 5, num_chains = min(4, dbarts::guessNumCores()),
                       num_threads_bart = min(num_chains, dbarts::guessNumCores()),
                       num_threads_wrangle = dbarts::guessNumCores(),
                       num_null_run = 10,
                       prior_power = 4, prior_base = 0.95,
                       time_points, break_spacing = "quantiles",
                       n_breaks = 8, lambda_vals = 10^seq(-6, 6, length.out = 1000),
                       min_method = "mean", min_quantile = NULL,
                       manual_breaks = NULL, alpha_g = NULL, alpha_l = NULL,
                       alpha_c = NULL) {

  # quick check
  if (missing(time_points)) {
    stop("Observed time points required.")
  }

  # quick check
  if (num_threads_bart > num_threads_wrangle) {
    stop("BART cores cannot be greater than available cores.")
  }

  # find max number of cores per worker
  max_cores = floor(num_threads_wrangle / num_threads_bart)

  # print out what we are running
  message(sprintf(
    "Completing %d null BART runs with up to %d null runs in parallel, each using %d cores.",
    num_null_run, max_cores, num_threads_bart
  ))

  # list of null outcomes
  y_permuted = replicate(num_null_run, sample(y), simplify = FALSE)

  # list of arg lists
  arg_lists = lapply(y_permuted, function(y_perm) {
    arg_list = list(y = y_perm, func_data_list = func_data_list,
                    num_trees = num_trees, num_samps = num_samps,
                    num_chains = num_chains, num_threads_bart = num_threads_bart,
                    prior_power = prior_power, prior_base = prior_base,
                    time_points = time_points, break_spacing = break_spacing,
                    n_breaks = n_breaks, lambda_vals = lambda_vals,
                    min_method = min_method)

    # conditionally add optional arguments
    if (!is.null(x)) arg_list$x = x
    if (!is.null(min_quantile)) arg_list$min_quantile = min_quantile
    if (!is.null(manual_breaks)) arg_list$manual_breaks = manual_breaks

    return(arg_list)
  })

  # timer start
  start = Sys.time()

  # parallelize BART runs - only use 1 thread if on windows
  if (.Platform$OS.type == "windows" || num_threads_wrangle == 1) {

    null_fbart_list = lapply(arg_lists, function(args) {
      tryCatch({
        do.call(fbart_for_select, args)
        cat("Done \n")
      }, error = function(e) {
        cat("Worker error:", conditionMessage(e), "\n")
        return(NULL)
      })
    })
  } else {
    null_fbart_list = parallel::mclapply(arg_lists, function(args) {
      tryCatch({
        do.call(fbart_for_select, args)
      }, error = function(e) {
        cat("Worker error:", conditionMessage(e), "\n")
        return(NULL)
      })
    }, mc.cores = max_cores)
  }

  # memory issues
  null_tip_list = lapply(null_fbart_list, function(x) x[["means"]])
  rm(arg_lists, y_permuted, null_fbart_list)
  gc()

  # timer stop
  stop = Sys.time()
  elapsed = stop - start

  # print time info
  message(sprintf("Amount of time for %d null fBART runs: %s.", num_null_run, format(elapsed)))
  message(sprintf("Starting real fBART run."))

  # convert to dataframe
  tip_df = as.data.frame(do.call(rbind, null_tip_list))
  to_return = list(tip_df = tip_df, elapsed_time_null = elapsed)

  # do real fbart run and get TIPs
  real_fbart = fbart_for_select(y = y, x = x, func_data_list = func_data_list,
                                num_trees = num_trees, num_samps = num_samps, num_burn = num_burn,
                                num_thin = num_thin, num_chains = num_chains,
                                num_threads_bart = num_threads_bart,
                                prior_power = prior_power, prior_base = prior_base,
                                time_points = time_points, break_spacing = break_spacing,
                                n_breaks = n_breaks, lambda_vals = lambda_vals,
                                min_method = min_method, manual_breaks = manual_breaks,
                                min_quantile = min_quantile)
  real_tip = colMeans(real_fbart[["indicators"]])
  to_return$real_tip = real_tip
  to_return$indicators = real_fbart[["indicators"]]
  to_return$basis = real_fbart[["basis"]]
  to_return$data_df = real_fbart[["data_df"]]

  # # get co-inclusions
  ind_sparse = Matrix(as.matrix(real_fbart[["indicators"]]), sparse = TRUE)
  co_occur = summary(t(ind_sparse) %*% ind_sparse)
  co_occur_unique = co_occur[co_occur$i != co_occur$j, ]
  co_occur_unique$var1 = colnames(real_fbart[["indicators"]])[co_occur_unique$i]
  co_occur_unique$var2 = colnames(real_fbart[["indicators"]])[co_occur_unique$j]

  co_occur_unique = co_occur_unique[co_occur_unique$var1 < co_occur_unique$var2, ]
  co_occur_unique$x = co_occur_unique$x / nrow(ind_sparse)
  to_return$tree_co_occur = co_occur_unique[order(-co_occur_unique$x), c("var1", "var2", "x")]
  colnames(to_return$tree_co_occur) = c("var1", "var2", "co_occurences")

  # memory issues
  rm(real_fbart, ind_sparse, co_occur, co_occur_unique)
  gc()

  # set to NULL and write over if non-NULL
  selected_global = selected_local = selected_global_se = NULL

  # global threshold
  if (!is.null(alpha_g)) {
    global_thresh = unname(quantile(as.matrix(tip_df), probs = 1 - alpha_g))
    to_return$global_thresh = global_thresh
    to_return$selected_global = names(real_tip[real_tip > global_thresh])

    grouped_terms = to_return$selected_global[grepl("_bspl_", to_return$selected_global)]
    to_return$functional_select_global = unique(sub("_.*", "", grouped_terms))
  }

  # local threshold
  if (!is.null(alpha_l)) {
    local_thresh = apply(tip_df, MARGIN = 2,
                         FUN = quantile, probs = (1 - alpha_l), na.rm = TRUE)
    to_return$local_thresh = unname(unlist(local_thresh))
    to_return$selected_local = names(real_tip[real_tip > local_thresh])

    grouped_terms = to_return$selected_local[grepl("_bspl_", to_return$selected_local)]
    to_return$functional_select_local = unique(sub("_.*", "", grouped_terms))
  }

  # global SE threshold
  if (!is.null(alpha_c)) {
    global_se_thresh = get_multiplier_thresh(record = tip_df, alpha = alpha_c)
    to_return$global_se_thresh = global_se_thresh
    to_return$selected_global_se = names(real_tip[real_tip > global_se_thresh])

    grouped_terms = to_return$selected_global_se[grepl("_bspl_", to_return$selected_global_se)]
    to_return$functional_select_global_se = unique(sub("_.*", "", grouped_terms))
  }

  # if no alphas provided calculate the global threshold with alpha = 0.01
  if (is.null(alpha_g) && is.null(alpha_l) && is.null(alpha_c)) {
    global_thresh = unname(quantile(as.matrix(tip_df), probs = 0.99))
    to_return$global_thresh = global_thresh
    to_return$selected_global = names(real_tip[real_tip > global_thresh])
  }

  return(to_return)
}
