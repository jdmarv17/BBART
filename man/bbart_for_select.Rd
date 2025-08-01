% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bbart_for_select.R
\name{bbart_for_select}
\alias{bbart_for_select}
\title{BBART for selection}
\usage{
bbart_for_select(
  y,
  x = NULL,
  func_data_list,
  num_trees = 20,
  num_samps = 5000,
  num_burn = 5000,
  num_thin = 5,
  num_chains = min(4, dbarts::guessNumCores()),
  num_threads_bart = min(num_chains, dbarts::guessNumCores()),
  prior_power = 4,
  prior_base = 0.95,
  time_points,
  break_spacing = "quantiles",
  n_breaks = 8,
  lambda_vals = 10^seq(-6, 6, length.out = 100),
  min_method = "mean",
  manual_breaks = NULL,
  min_quantile = NULL
)
}
\arguments{
\item{y}{A vector of outcome values.}

\item{x}{A dataframe of scalar covariates. If 'nxp' then 'y' must have length 'n' and each dataframe in 'func_data_list' must have 'n' columns.}

\item{func_data_list}{List of functional covariates. Each list element should be a dataframe with rows indicating time points and columns for each subject. Importantly, each element of the list of functional covariates must be named.}

\item{num_trees}{Number of BART trees per sample.}

\item{num_samps}{Number of posterior samples to take after burn-in.}

\item{num_burn}{Number of burn-in samples.}

\item{num_thin}{Keep every 'num_thin' draw.}

\item{num_chains}{Number of independent MCMC chains to run.}

\item{num_threads_bart}{Number of threads used to run BART - not recommended to be greater than 'num_chains'.}

\item{prior_power}{Power parameter for tree prior.}

\item{prior_base}{Base parameter for tree prior.}

\item{time_points}{Vector of observed time points. Method assumes all subjects and functional covariates are observed at a common set of time points.}

\item{break_spacing}{A value in c("quantiles", "equal", "manual") to set break spacing for B-spline basis.}

\item{n_breaks}{A number indicating number of breaks to use for B-spline basis (includes endpoints, results in n_breaks + 2 B-spline polynomials).}

\item{lambda_vals}{A vector of positive numbers to cross-validate over to get roughness penalty.}

\item{min_method}{A value to minimize during cross-validation, either "mean" or "quantile" (over the subjects GCV values)}

\item{manual_breaks}{An optional vector of manual break values. Must include domain endpoints.}

\item{min_quantile}{if min_method = "quantile", a value for quantile to use}

\item{num_threads_wrangle}{Number of threads for post processing - can be larger than 'num_chains'.}
}
\value{
A list of seven objects. Object (1) is the marginal TIP for variables in the real BBART run. Object (2)
an indicators data frame of variable inclusion in each posterior tree. Object (3) is the dataframe of BBART splits and
terminal nodes. Object (4) is the dbarts::bart2() object from the real BBART run. Object (5) is the B-spline basis
used to represent functional covariates. Object (6) are the lambdas from the weighted fit process for the basis
expansion. Object (7) is the appended data matrix with scalar outcome, scalar covariates, and fit B-spline coefficients.
}
\description{
BBART for selection
}
