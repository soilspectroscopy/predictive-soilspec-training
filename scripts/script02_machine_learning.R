# ## ----libraries, message=FALSE, warning=FALSE--------------------------------------
# library("tidyverse")
# library("mlr3verse")
# library("Cubist") # Cubist ML algorithm
# library("ranger") # Random Forest ML algorithm
# library("glmnet") # Elastic net ML Algorithm
# library("future") # Parallel processing
# library("yardstick") # Additional metrics
# library("mlr3extralearners") # Additional models and features
# library("iml") # Feature importance
# 
# 
# ## ----eval=FALSE-------------------------------------------------------------------
# ## # latest GitHub release
# ## remotes::install_github("mlr-org/mlr3extralearners@*release")
# 
# 
# ## ----working_directory------------------------------------------------------------
# my.wd <- "~/projects/local-files/soilspec_training"
# # Or within an RStudio project
# # my.wd <- getwd()
# 
# 
# ## ----setup, message=FALSE, warning=FALSE------------------------------------------
# #install.packages("mlr3verse")
# library("mlr3verse")
# library("tidyverse")
# 
# # Loading train and test data
# train_matrix <- read_csv(file.path(my.wd, "train.csv"))
# test_matrix <- read_csv(file.path(my.wd, "test.csv"))
# 
# 
# ## ----R6_object_example, message=FALSE, warning=FALSE------------------------------
# # Creating an example object
# example_object <- TaskRegr$new(id = "NIR_spectral",
#                               backend = train_matrix, 
#                               target = "oc_usda.c729_w.pct")
# 
# # Access to the $id field
# example_object$id
# 
# # Calling methods of the object
# example_object$set_col_roles("id.sample_local_c", role = "name")
# example_object$set_col_roles("location.country_iso.3166_txt", role = "group")
# example_object$col_roles
# 
# 
# ## ----histograms, message=FALSE, warning=FALSE-------------------------------------
# hist(train_matrix$oc_usda.c729_w.pct, breaks = 100)
# 
# 
# ## ----log1p, message=FALSE, warning=FALSE------------------------------------------
# train_matrix$oc_usda.c729_w.pct <- log1p(train_matrix$oc_usda.c729_w.pct)
# test_matrix$oc_usda.c729_w.pct <- log1p(test_matrix$oc_usda.c729_w.pct)
# hist(train_matrix$oc_usda.c729_w.pct, breaks = 100)
# 
# 
# ## ----task, message=FALSE, warning=FALSE-------------------------------------------
# # Setting a seed for reproducibility
# set.seed(349)
# 
# # Create a task
# train_matrix <- train_matrix %>%
#   select(-location.country_iso.3166_txt)
# 
# task_neospectra <- as_task_regr(train_matrix, 
#                               target = "oc_usda.c729_w.pct",
#                               id = "NIR_spectral")
# 
# # Set row names role to the id.sample_local_c column 
# task_neospectra$set_col_roles("id.sample_local_c", roles = "name")
# task_neospectra
# 
# 
# ## ----featureless_baseline, message=FALSE, warning=FALSE---------------------------
# # Load featureless learner
# # featureless = always predicts new values as the mean of the dataset
# lrn_baseline <- lrn("regr.featureless")
# 
# # Only one hyperparameter can be set: robust = calculate mean if false / calculate median if true
# # lrn_baseline$param_set
# 
# # Train learners using $train() method
# lrn_baseline$train(task_neospectra)
# 
# 
# ## ----featureless_baseline_evaluation, message=FALSE, warning=FALSE----------------
# # Setting a seed for reproducibility
# set.seed(349)
# 
# # Define goodness-of-fit metrics
# measures <- msrs(c("regr.rmse", "regr.bias", "regr.rsq"))
# 
# # Predict on a new data using $predict_newdata() method
# pred_baseline <- lrn_baseline$predict_newdata(test_matrix)
# 
# # Evaluation on the new data applying $score() method
# pred_baseline$score(measures)
# 
# 
# ## ----learners_cv_benchmarking, message=FALSE, warning=FALSE-----------------------
# lrn_cv <- lrns(c("regr.ranger", "regr.cv_glmnet", "regr.cubist"))
# cv10 <- rsmp("cv", folds = 10)
# 
# 
# ## ----simple_benchmarking, message=FALSE, warning=FALSE, results=FALSE-------------
# bench_cv_grid <- benchmark_grid(task_neospectra, lrn_cv, cv10)
# bench_cv <- benchmark(bench_cv_grid)
# 
# 
# ## ----simple benchmarking results, message=FALSE, warning=FALSE--------------------
# # Check the results of each fold using $score()
# bench_cv$score(measures)
# 
# # Check the aggregated results of the models using $aggregate()
# bench_cv$aggregate(measures)
# 
# 
# ## ----tune_cubist_model, message=FALSE, warning=FALSE, results=FALSE---------------
# # Define search space for committees and neighbors inside the Learner object using to_tune function
# lrn_cubist <- lrn("regr.cubist",
#                  committees = to_tune(1, 20),
#                  neighbors = to_tune(0, 5))
# 
# # Setting a seed for reproducibility
# set.seed(349)
# 
# # Define tuning instance using tune() helper function
# # tune() function will automatically calls $optimize() method of the instance
# simple_tune <- tune(
#   method = tnr("random_search"), # Picking random combinations of hyperparameters from the search space
#   task = task_neospectra, # The task with training data
#   learner = lrn_cubist, # Learner algorithm
#   resampling = cv10, # Apply 10-fold CV to calculate the evaluation metrics
#   measures = msr("regr.rmse"), # Use RMSE to select the best set of hyperparameters
#   term_evals = 10 # Terminate tuning after testing 10 combinations
# )
# 
# 
# ## ----tune_cubist_model_results, message=FALSE, warning=FALSE----------------------
# # Exploring tuning results
# simple_tune$result
# 
# 
# ## ----showing_tunig_spaces, message=FALSE, warning=FALSE---------------------------
# # Showing dictionary of available pre-trained tuning spaces
# as.data.table(mlr_tuning_spaces)[,.(label)]
# 
# 
# ## ----inner_cv_tuning, message=FALSE, warning=FALSE, results=FALSE-----------------
# # Selecting the spaces using sugar function ltss()
# # Not available for cubist
# spaces <- ltss(c("regr.ranger.default", "regr.glmnet.default"))
# 
# # Creating a list of learners and an empty list to write tuned learners
# lrn_cv <- lrns(c("regr.ranger", "regr.glmnet"))
# lrn_cv_tuned <- list()
# 
# # To avoid most of the logs
# lgr::get_logger("mlr3")$set_threshold("warn")
# 
# # For loop to tune hyperparameters of every learner using random search
# i=1
# for(i in 1:length(lrn_cv)){
# 
#   # Using multiple cores for tuning each model
#   future::plan("multisession")
# 
#   set.seed(349)
#   instance <- tune(task = task_neospectra,
#                    method = tnr("random_search"),
#                    learner = lrn_cv[[i]],
#                    resampling = cv10,
#                    measures = msr("regr.rmse"),
#                    term_evals = 10,
#                    search_space = spaces[[i]])
# 
#   # Writing tuned hyperparameters
#   lrn_cv_tuned[[i]] <- lrn_cv[[i]]
#   lrn_cv_tuned[[i]]$param_set$values <- instance$result_learner_param_vals
# 
# }
# 
# # Close parallel connection
# future:::ClusterRegistry("stop")
# 
# # Adding previously tuned cubist model to the list for benchmarking
# # Create a new instance of the cubist model
# lrn_cubist_tuned <- lrn("regr.cubist")
# 
# # Set the best performing hyperparameters to the new instance
# lrn_cubist_tuned$param_set$values <- simple_tune$result_learner_param_vals
# 
# # Bind it to the tuned models lists
# lrn_cv_tuned <- c(lrn_cv_tuned, lrn_cubist_tuned)
# 
# 
# ## ----outer_cv_benchmarking, message=FALSE, warning=FALSE, results=FALSE-----------
# # Creating a benchmark grid
# bench_cv_grid <- benchmark_grid(task_neospectra, lrn_cv_tuned, cv10)
# 
# # Implementing bechmarking using multiple cores
# future::plan("multisession")
# bench_cv <- benchmark(bench_cv_grid)
# future:::ClusterRegistry("stop")
# 
# 
# ## ----benchmarking_results, message=FALSE, warning=FALSE---------------------------
# # Evaluating the results
# bench_cv$aggregate(measures)
# 
# 
# ## ----final_cubist_model, message=FALSE, warning=FALSE-----------------------------
# # New Cubist learner
# lrn_cubist_tuned <- lrn("regr.cubist")
# 
# # Setting the best HPs
# lrn_cubist_tuned$param_set$values <- lrn_cv_tuned[[3]]$param_set$values # Cubist is the third learner in the list
# 
# # Fitting with whole train set
# lrn_cubist_tuned$train(task_neospectra)
# 
# # Predicting on test set
# pred_cubist_response <- lrn_cubist_tuned$predict_newdata(test_matrix)
# 
# # Final goodness-of-fit metrics
# pred_cubist_response$score(measures)
# 
# 
# ## ----importance, message=FALSE, warning=FALSE-------------------------------------
# # Grabbing the features in train data
# model_predictors <- task_neospectra$data(cols = task_neospectra$feature_names)
# 
# # Grabbing the outcome in train data
# model_outcome <- task_neospectra$data(cols = task_neospectra$target_names)
# 
# # Building a predictor object for importance evaluation
# predictor = Predictor$new(lrn_cubist_tuned,
#                           data = model_predictors,
#                           y = model_outcome)
# 
# importance = FeatureImp$new(predictor, loss = "mse", n.repetitions = 10)
# 
# importance$plot()
# 
# 
# ## ----accuracy_plot, message=FALSE, warning=FALSE----------------------------------
# accuracy.data <- tibble(id = test_matrix$id.sample_local_c,
#                     observed = pred_cubist_response$data$truth,
#                     predicted = pred_cubist_response$data$response)
# 
# ggplot(accuracy.data) +
#   geom_point(aes(x = observed, y = predicted)) +
#   geom_abline(intercept = 0, slope = 1) +
#   theme_light()
# 
# 
# ## ----accuracy_metrics, message=FALSE, warning=FALSE-------------------------------
# library("yardstick")
# 
# # Calculating metrics
# accuracy.metrics <- accuracy.data %>%
#   summarise(n = n(),
#               rmse = rmse_vec(truth = observed, estimate = predicted),
#               bias = msd_vec(truth = observed, estimate = predicted),
#               rsq = rsq_vec(truth = observed, estimate = predicted),
#               ccc = ccc_vec(truth = observed, estimate = predicted, bias = T),
#               rpiq = rpiq_vec(truth = observed, estimate = predicted))
# 
# accuracy.metrics
# 
# 
# ## ----uncertainty, message=FALSE, warning=FALSE------------------------------------
# # Runing 10 CV with the Cubist tuned model
# cv_results <- resample(task = task_neospectra,
#                        learner = lrn_cubist_tuned,
#                        resampling = rsmp("cv", folds = 10))
# 
# # Extracting the observed and predicted values within the folds
# cv_results <- lapply(1:length(cv_results$predictions("test")), function(i){
#     as.data.table(cv_results$predictions("test")[[i]]) %>%
#       mutate(fold = i)})
# 
# cv_results <- Reduce(rbind, cv_results)
# head(cv_results)
# 
# # Quick visualization
# ggplot(cv_results) +
#   geom_point(aes(x = truth, y = response)) +
#   geom_abline(intercept = 0, slope = 1) +
#   theme_light()
# 
# # Estimating the absolute error
# cv_results <- cv_results %>%
#   mutate(error = abs(truth-response))
# 
# # Create a matrix with the same predictors as train_matrix and
# # variable error as the new target
# # you can View to see if the join worked well
# error_matrix <- train_matrix %>%
#   mutate(row_ids = row_number()) %>%
#   left_join(cv_results, by = "row_ids")
# 
# # Simplifying with only necessary data
# error_matrix <- error_matrix %>%
#   select(id.sample_local_c, error, starts_with("pc_"))
# 
# task_error <- as_task_regr(error_matrix, 
#                           target = "error",
#                           id = "error_model")
# 
# # Set row names role to the id.sample_local_c column 
# task_error$set_col_roles("id.sample_local_c", roles = "name")
# task_error
# 
# # Retrain the same model structure with the task_error
# lrn_cubist_error = lrn("regr.cubist")
# lrn_cubist_error$param_set$values = lrn_cubist_tuned$param_set$values
# 
# # Retrain
# lrn_cubist_error <- lrn_cubist_tuned$train(task_error)
# 
# # Get the plain error predictions
# error_predictions <- predict(lrn_cubist_error, newdata = as.data.table(task_error))
# 
# # Calculating the conformity scores (alpha) of the calibration set
# error_matrix <- error_matrix %>%
#     select(-starts_with("pc_")) %>%
#     mutate(pred_error = error_predictions) %>%
#     mutate(alpha_scores = error/pred_error)
# 
# error_matrix
# 
# # Sample correction and final conformity score
# target_coverage = 0.95
# n_calibration <- nrow(error_matrix)
# corrected_quantile <- ((n_calibration + 1)*(target_coverage))/n_calibration
# alpha_corrected <- quantile(error_matrix$alpha_scores, corrected_quantile)
# 
# # Predicting error on test set and calculating 95% prediction interval
# # Predicting on test set
# pred_cubist_error = lrn_cubist_error$predict_newdata(test_matrix)
# pred_cubist_error
# 
# # Final output
# final.outputs <- accuracy.data %>%
#   mutate(error = pred_cubist_error$response) %>%
#   mutate(PI = error*alpha_corrected) %>%
#   mutate(lower_PI = predicted-PI,
#          upper_PI = predicted+PI) %>%
#   select(-error, -PI)
# 
# # Back-transforming all results because we have used log1p earlier
# final.outputs <- final.outputs %>%
#   mutate(observed = expm1(observed),
#          predicted = expm1(predicted),
#          lower_PI = expm1(lower_PI),
#          upper_PI = expm1(upper_PI))
# 
# # Calculating coverage statistics
# 
# final.outputs <- final.outputs %>%
#   mutate(covered = ifelse(observed >= lower_PI & observed <= upper_PI, TRUE, FALSE),
#          wis = ifelse(observed < lower_PI,
#                         (upper_PI-lower_PI)+(2/(1-target_coverage))*(lower_PI-observed),
#                         ifelse(observed > upper_PI,
#                               (upper_PI-lower_PI)+(2/(1-target_coverage))*(observed-upper_PI),
#                                upper_PI-lower_PI)))
# 
# 
# ## ----final_outputs, message=FALSE, warning=FALSE----------------------------------
# # Coverage statistics
# final.outputs %>%
#   count(covered) %>%
#   mutate(perc = n/sum(n)*100)
# 
# # Final accuracy metrics in the original range
# final.accuracy <- final.outputs %>%
#   summarise(n = n(),
#             rmse = rmse_vec(truth = observed, estimate = predicted),
#             bias = msd_vec(truth = observed, estimate = predicted),
#             rsq = rsq_trad_vec(truth = observed, estimate = predicted),
#             ccc = ccc_vec(truth = observed, estimate = predicted, bias = T),
#             rpiq = rpiq_vec(truth = observed, estimate = predicted),
#             coverage = sum(covered)/n*100,
#             wis = mean(wis))
# 
# final.accuracy
# 
# # Final accuracy plot
# perfomance.annotation <- paste0("Rsq = ", round(final.accuracy[[1,"rsq"]], 2),
#                                 "\nRMSE = ", round(final.accuracy[[1,"rmse"]], 2), " wt%",
#                                 "\nCoverage = ", round(final.accuracy[[1,"coverage"]], 2), "%",
#                                 "\nWIS = ", round(final.accuracy[[1,"wis"]], 2), " wt%")
# 
# p.final <- ggplot(final.outputs) +
#   geom_pointrange(aes(x = observed, y = predicted,
#                       ymin = lower_PI, ymax = upper_PI,
#                       color = covered)) +
#   geom_abline(intercept = 0, slope = 1) +
#   geom_text(aes(x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2),
#                 label = perfomance.annotation, size = 3) +
#   labs(title = "Soil Organic Carbon (wt%) test prediction",
#        x = "Observed", y = "Predicted", color = "PI95% covered:") +
#   theme_light() + theme(legend.position = "bottom")
# 
# # Making sure it we have a square plot
# r.max <- max(layer_scales(p.final)$x$range$range)
# r.min <- min(layer_scales(p.final)$x$range$range)
# 
# s.max <- max(layer_scales(p.final)$y$range$range)
# s.min <- min(layer_scales(p.final)$y$range$range)
# 
# t.max <- round(max(r.max,s.max),1)
# t.min <- round(min(r.min,s.min),1)
# 
# p.final <- p.final + coord_equal(xlim = c(t.min,t.max), ylim = c(t.min,t.max))
# p.final
# 
