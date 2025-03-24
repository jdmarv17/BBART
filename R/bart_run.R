#' Function to run dbarts::bart2
#'
#' @param formula Object of class formula.
#' @param data A data frame of original data.
#' @param num_trees Number of BART trees per sample.
#' @param num_samps Number of posterior samples to take after burn-in.
#' @param num_burn Number of burn-in samples.
#' @param num_thin Keep every 'num_thin' draw.
#' @param num_chains Number of independent MCMC chains to run.
#' @param num_threads_bart Number of threads used to run BART - not recommended to be greater than 'num_chains'.
#' @param num_threads_wrangle Number of threads for post processing - can be larger than 'num_chains'.
#' @param prior_power Power parameter for tree prior.
#' @param prior_base Base parameter for tree prior.
#'
#' @return A list of three objects. Object (1) is an indicator data frame of variable tree inclusion, object (2) is a data frame of posterior trees without terminal nodes, and object (3) is the dbarts::bart2 object
#'
#' @importFrom rlang .data
#' @export
#'
#'
#' @examples
#' library(MASS)
#' library(dbarts)
#' library(dplyr)
#' library(janitor)
#' data(birthwt)
#' data = birthwt %>% dplyr::select(., -low)
#' formula = bwt ~ .
#' dbarts_data = dbartsData(formula, data)
#' dbarts_x = clean_names(dbarts_data@x)
#' dbarts_y = dbarts_data@y
#' vars = colnames(dbarts_x)
#' bart_run = run_bart(formula, data, vars = vars)
#'
run_bart = function(formula, data,
                    num_trees = 50, num_samps = 5000, num_burn = 5000,
                    num_thin = 5, num_chains = min(4, dbarts::guessNumCores()),
                    num_threads_bart = min(num_chains, dbarts::guessNumCores()),
                    num_threads_wrangle = dbarts::guessNumCores(),
                    prior_power = 4, prior_base = 0.95) {


  # grab response in case we need sd
  response = as.character(formula[2])

  # if n <= p get a sigma estimate from std dev of the response, otherwise let bart2() estimate it
  sigma = ifelse(nrow(data) <= (ncol(data) - 1), sd(data[[response]], na.rm = TRUE), NA_real_)

  # run bart2()
  bart_run = dbarts::bart2(formula, data = data,
                           n.trees = num_trees, n.samples = num_samps,
                           n.burn = num_burn, n.thin = num_thin, n.chains = num_chains,
                           n.threads = num_threads_bart,
                           power = prior_power, base = prior_base,
                           keepTrees = TRUE, printCutoffs = 0, verbose = FALSE,
                           sigest = sigma)

  # extract the trees
  bart_trees = dbarts::extract(bart_run, type = c("trees"))

  # get rid of terminal nodes (THIS IS FIRST RETURN OBJECT)
  filtered_trees = bart_trees[bart_trees$var != -1,]

  # split up by chain, tree, sample to take advantage of mclapply
  var_list =
    filtered_trees %>%
    dplyr::group_by(., chain, tree, sample) %>%
    dplyr::group_split() %>% # list of filtered_tree tibbles by chain, tree, sample
    parallel::mclapply(., function(x) x[,"var"], mc.cores = num_threads_wrangle) %>% # select only var column
    lapply(., unlist) %>% # get as vector
    lapply(., unname) # get rid of names

  # run through make_indicators() with mclapply
  indicator_list = parallel::mclapply(X = var_list, FUN = make_indicators,
                                      p = (ncol(data) - 1), mc.cores = num_threads_wrangle)

  # get indicators and set column names appropriately
  indicator_df =
    as.data.frame(data.table::transpose(indicator_list))
  colnames(indicator_df) = colnames(bart_run[["fit"]][["data"]]@x)


  to_return = list(indicator_df, filtered_trees, bart_run)

  return(to_return)
}

