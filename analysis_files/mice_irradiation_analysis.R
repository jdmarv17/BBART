library(fda)
library(dbarts)
library(dplyr)
library(refund)
library(mgcv)
library(doParallel)
library(parallel)
library(caret)
library(randomForest)
library(fBART)


load("analysis_files/data/mice_irradiation.rda")
data_list = mice_irradiation$mice_covars
outcome = mice_irradiation$y
time_points = 1:nrow(data_list[[1]])

set.seed(30625)
cv_folds = 5
folds = rep(1:cv_folds, each = 35)
folds = sample(folds)  # Shuffle the indices

multivar_record = as.data.frame(matrix(nrow = cv_folds, ncol = 3))
colnames(multivar_record) = c("fbart", "rf", "pfr")

for (i in 1:cv_folds) {

  train_ind = which(folds != i)

  # b spline representation
  bspline = get_bspline_coefs(func_data_list = data_list[c("x1", "x2", "x3", "x4", "x5")],
                              time_points = time_points,
                              break_spacing = "quantiles", n_breaks = 8,
                              lambda_vals = 10^seq(-6, 6, length.out = 250), min_method = "mean",
                              num_threads = 8)
  breaks = bspline$basis$params
  data_df = cbind(y = outcome, bspline$coef_df)
  data_df_train = data_df[train_ind, ]
  data_df_test = data_df[-train_ind, ]

  # pfr data
  data_mats_train = lapply(data_list[c("x1", "x2", "x3", "x4", "x5")], function(x) x[, train_ind])
  data_mats_test = lapply(data_list[c("x1", "x2", "x3", "x4", "x5")], function(x) x[, -train_ind])
  data_mats_pfr_train = lapply(data_mats_train, function(x) as.matrix(data.table::transpose(x)))
  data_mats_pfr_test = lapply(data_mats_test, function(x) as.matrix(data.table::transpose(x)))
  data_pfr_train = c(list(y = outcome[train_ind]), data_mats_pfr_train)
  data_pfr_test = c(list(y = outcome[-train_ind]), data_mats_pfr_test)

  formula = paste0("y ~ ",
                   paste("lf(", names(data_mats_train),
                         ", argvals = time_points, bs = 'bs', xt = list(knots = breaks), m = 4)",
                         sep = "", collapse = " + "))


  # run bart
  bart_run = run_bart(y ~ ., data = data_df_train,
                      num_trees = 150, num_samps = 10000, num_burn = 10000,
                      num_thin = 5, num_chains = 4, num_threads_bart = 4,
                      num_threads_wrangle = 8, prior_power = 4, prior_base = 0.95)
  bart_mod = bart_run[[3]]
  pred_bart = tibble(prob = apply(predict(bart_mod, data_df_test, type = "ev"), 2, mean),
                     observed = data_df_test$y) %>%
    mutate(pred = ifelse(prob > 0.5, 1, 0))
  pfc_bart = mean(pred_bart$pred != pred_bart$observed)


  # run random forest, get predictions
  rf_model = randomForest(factor(y) ~ ., data_df_train, ntree = 150)

  pred_rf = tibble(pred = predict(rf_model, data_df_test),
                   observed = data_df_test$y)
  pfc_rf = mean(pred_rf$pred != pred_rf$observed)

  # run refund::pfr(), get predictions
  pfr_model = refund::pfr(as.formula(formula), family = binomial(link = "logit"), data = data_pfr_train)
  pred_pfr = tibble(prob = predict(pfr_model, newdata = data_pfr_test, type = "response"),
                    observed = data_pfr_test$y) %>%
    mutate(pred = ifelse(prob > 0.5, 1, 0))
  pfc_pfr = mean(pred_pfr$pred != pred_pfr$observed)



  multivar_record[i,] = c(pfc_bart, pfc_rf, pfc_pfr)
  cat("CV fold", i, "completed \n")

}

colMeans(multivar_record)






