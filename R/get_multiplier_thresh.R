#' Generate global multiplier thresholds
#'
#' @param A A data frame of inclusion proportions ('p' columns and 'num_run' rows)
#' @param alpha Determines 1-alpha quantile in multiplier method
#'
#' @returns A vector of variable specific global multiplier thresholds
#'
#'
get_multiplier_thresh = function(record, alpha = 0.01) {
  ##### This function will take a probability record and an alpha and return a global
  ##### multiplier C* according to Bleich et al. (2014) global SE strategy
  ##### Return 1: df with 1 - alpha quantiles, mean, sd, and multiplier thresh for each variable
  # ARGS:
  # `record` = record of null inclusion proportions
  #          = contains 'p' columns and 'num_run' rows
  # `alpha` = determines quantile in multiplier method


  # get means and sd's for each
  means = colMeans(record, na.rm = TRUE)
  std_dev = apply(record, MARGIN = 2, FUN = sd, na.rm = TRUE)

  # get list of variable specific quantiles
  quant_record = apply(record, MARGIN = 2, FUN = quantile, probs = (1 - alpha), na.rm = TRUE)

  # get multipliers
  stats =
    data.frame(quant_record, means, std_dev) %>%
    mutate(C = (quant_record - means) / std_dev,
           C_thresh = (means + max(C) * std_dev)) %>%
    mutate(C_thresh = case_when(
      C_thresh < 0 ~ 0, # if less than 0 set to 0
      C_thresh > 1 ~ 1, # if greater than 1 set to 1
      TRUE ~ C_thresh # otherwise leave the same
    ))

  # return C's and max{C_i}
  # to_return = list(stats)
  c = stats$C_thresh
  return(c)
}
