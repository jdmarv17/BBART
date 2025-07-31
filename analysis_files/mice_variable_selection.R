library(fda)
library(fda.usc)
library(dbarts)
library(tidyverse)
library(parallel)
library(doParallel)
library(Bolstad2)
library(sparsegl)
library(refund)
library(zoo)
library(Matrix)
library(BBART)

load("analysis_files/data/mice_irradiation_5_12.rda")
data_list = lapply(combined_list, function(x) x[, !(names(x) %in% "week")])
names(data_list) = c("x1", "x2", "x3", "x4", "x5")
radiation = substring(colnames(data_list[[1]]), 2, 2)
outcome = ifelse(radiation == "0", 0, 1)
time_points = 1:nrow(data_list[[1]])

# impute missing values with interpolation
data_list_impute = lapply(data_list, function(df) {
  as.data.frame(lapply(df, function(col) {
    na.approx(col, x = time_points, na.rm = FALSE, rule = 2)
  }))
})


# set seed and get CV folds
set.seed(61225)
fold_vec = sample(c(rep(1:9, each = 45), rep(10, 41)), replace = FALSE)

results = tibble(sim_num = integer(0), method = character(0),
                 type = character(0), selected = list())
pred_results = tibble(sim_num = integer(0), method = character(0), pfc = double(0))

for (i in 1:length(unique(fold_vec))) {

  # get test and train data
  keep_cols = fold_vec != i
  train_list = purrr::map(data_list_impute, ~ select(.x, which(keep_cols)))
  outcome_train = ifelse(substring(colnames(train_list[[1]]), 2, 2) == "0", 0, 1)

  keep_cols = fold_vec == i
  test_list = purrr::map(data_list_impute, ~ select(.x, which(keep_cols)))
  outcome_test = ifelse(substring(colnames(test_list[[1]]), 2, 2) == "0", 0, 1)


  # run BART selection
  select =
    BBARTselect(y = outcome_train, x = NULL, func_data_list = train_list,
                num_trees = 10, num_samps = 10000, num_burn = 5000,
                num_thin = 5, num_chains = 4,
                num_threads_bart = 4,
                num_threads_wrangle = 16,
                num_null_run = 10,
                time_points = time_points,
                prior_power = 4, prior_base = 0.99,
                break_spacing = "quantiles",
                n_breaks = 10, lambda_vals = 10^seq(-5, 5, length.out = 1000),
                min_method = "mean",
                alpha_l = 0.0025, alpha_g = 0.0025, alpha_c = 0.0025)

  grouped_terms = select$selected_global[grepl("_bspl_", select$selected_global)]
  functional_select = unique(sub("_.*", "", grouped_terms))

  row1 = tibble(sim_num = i, method = "bBART", type = "functional", selected = list(functional_select))
  row2 = tibble(sim_num = i, method = "bBART", type = "component", selected = list(grouped_terms))


  # run sparse group lasso selection
  x_lasso = as.matrix(select$data_df)
  x_lasso = x_lasso[, colnames(x_lasso) != "y"]
  y_lasso = select$data_df$y
  groups = rep(1:5, each = 12)
  lasso_cv = cv.sparsegl(x = x_lasso, y = y_lasso, group = groups, asparse = 0.25)
  lasso_coef = coef(lasso_cv, s = "lambda.1se")
  lasso_selected = lasso_coef@Dimnames[[1]][lasso_coef@i + 1]
  lasso_selected = lasso_selected[lasso_selected != "(Intercept)"]

  grouped_terms = lasso_selected[grepl("_bspl_", lasso_selected)]
  functional_select = unique(sub("_.*", "", grouped_terms))

  row3 = tibble(sim_num = i, method = "SG_LASSO", type = "functional", selected = list(functional_select))
  row4 = tibble(sim_num = i, method = "SG_LASSO", type = "component", selected = list(grouped_terms))


  # run group lasso selection
  lasso_cv = cv.sparsegl(x = x_lasso, y = y_lasso, group = groups, asparse = 0)
  lasso_coef = coef(lasso_cv, s = "lambda.1se")
  lasso_selected = lasso_coef@Dimnames[[1]][lasso_coef@i + 1]
  lasso_selected = lasso_selected[lasso_selected != "(Intercept)"]

  grouped_terms = lasso_selected[grepl("_bspl_", lasso_selected)]
  functional_select = unique(sub("_.*", "", grouped_terms))

  row5 = tibble(sim_num = i, method = "G_LASSO", type = "functional", selected = list(functional_select))

  results = bind_rows(results, row1, row2, row3, row4, row5)




  # run FLM models for each methods selections
  breaks = c(1, select$basis$params, 144)
  data_mats_pfr_train = lapply(train_list, function(x) as.matrix(data.table::transpose(x)))
  data_pfr_train = c(list(y = outcome_train), data_mats_pfr_train)
  data_mats_pfr_test = lapply(test_list, function(x) as.matrix(data.table::transpose(x)))
  data_pfr_test = c(list(y = outcome_test), data_mats_pfr_test)

  # bart selected model
  formula = paste0("y ~ ", paste("lf(", unlist(row1$selected), ", argvals = time_points, bs = 'bs', xt = list(knots = breaks), m = 4)",
                                 sep = "", collapse = " + "))

  pfr_model = refund::pfr(as.formula(formula), data = data_pfr_train, family = binomial(link = "logit"))
  pred_bart_test = tibble(prob = predict(pfr_model, data_pfr_test),
                         observed = data_pfr_test$y) %>%
    mutate(pred = ifelse(prob > 0.5, 1, 0))
  pfc_bart_test = mean(pred_bart_test$pred != pred_bart_test$observed)

  # SG_LASSO selected model
  formula = paste0("y ~ ", paste("lf(", unlist(row3$selected), ", argvals = time_points, bs = 'bs', xt = list(knots = breaks), m = 4)",
                                 sep = "", collapse = " + "))

  pfr_model = refund::pfr(as.formula(formula), data = data_pfr_train, family = binomial(link = "logit"))
  pred_sglasso_test = tibble(prob = predict(pfr_model, data_pfr_test),
                          observed = data_pfr_test$y) %>%
    mutate(pred = ifelse(prob > 0.5, 1, 0))
  pfc_sglasso_test = mean(pred_sglasso_test$pred != pred_sglasso_test$observed)

  # G_LASSO selected model
  formula = paste0("y ~ ", paste("lf(", unlist(row5$selected), ", argvals = time_points, bs = 'bs', xt = list(knots = breaks), m = 4)",
                                 sep = "", collapse = " + "))

  pfr_model = refund::pfr(as.formula(formula), data = data_pfr_train, family = binomial(link = "logit"))
  pred_glasso_test = tibble(prob = predict(pfr_model, data_pfr_test),
                             observed = data_pfr_test$y) %>%
    mutate(pred = ifelse(prob > 0.5, 1, 0))
  pfc_glasso_test = mean(pred_glasso_test$pred != pred_glasso_test$observed)

  new_rows = tibble(sim_num = rep(i, 3),
                    method = c("bBART", "SG_LASSO", "G_LASSO"),
                    pfc = c(pfc_bart_test, pfc_sglasso_test, pfc_glasso_test))
  pred_results = bind_rows(pred_results, new_rows)


  cat("Fold", i, "complete. \n")
}

pred_results %>% group_by(method) %>% summarise(mean = mean(pfc))


var_key = tibble(var = colnames(x_lasso))

tables =
  results %>%
  tidyr::unnest(selected) %>%
  count(method, type, selected) %>%
  group_by(method, type) %>%
  arrange(desc(n), .by_group = TRUE) %>%
  group_split()

(tables_comp = tables[c(2, 4)])
(tables_func = tables[c(1, 3, 5)])




