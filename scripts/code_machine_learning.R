## ----libraries, message=FALSE, warning=FALSE------------------------------------------------------------------------------
library("tidyverse")
library("qs")         # Loading serialized files
library("prospectr")  # prep preprocessing
library("tidymodels")
library("probably")   # Conformal prediction intervals
library("rules")
library("Cubist")      # Cubist engine
library("kernlab")    # SVM engine
library("ranger")     # Random forest engine
library("yardstick")  # Performance metrics
tidymodels_prefer()   # Resolve common function name conflicts


## ----neospectra_read, eval=TRUE-------------------------------------------------------------------------------------------
## Internet configuration for downloading large datasets
options(timeout = 10000)

## Reading serialized files using qs::qread_url()
neospectra.soil <- qread_url("https://storage.googleapis.com/soilspec4gg-public/neospectra_soillab_v1.2.qs")
dim(neospectra.soil)

neospectra.site <- qread_url("https://storage.googleapis.com/soilspec4gg-public/neospectra_soilsite_v1.2.qs")
dim(neospectra.site)

neospectra.nir <- qread_url("https://storage.googleapis.com/soilspec4gg-public/neospectra_nir_v1.2.qs")
dim(neospectra.nir)


## ----neospectra_filter, eval=TRUE-----------------------------------------------------------------------------------------
# How many samples from each country?
# Most samples are from USA, others from African countries
# We will use this column to define the train/test split
neospectra.site %>%
  count(location.country_iso.3166_txt)

# Selecting relevant site data
neospectra.site <- neospectra.site %>%
  select(id.sample_local_c, location.country_iso.3166_txt)

# Selecting relevant soil data
neospectra.soil <- neospectra.soil %>%
  select(id.sample_local_c, oc_usda.c729_w.pct)

# Selecting relevant NIR data and taking average across repeats
neospectra.nir <- neospectra.nir %>%
  select(id.sample_local_c, starts_with("scan_nir")) %>%
  group_by(id.sample_local_c) %>%
  summarise(across(everything(), mean)) %>%
  ungroup()

# Renaming spectral column headers to plain numeric wavelength values
old.names <- neospectra.nir %>%
  select(starts_with("scan_nir")) %>%
  names()

new.names <- gsub("scan_nir.|_ref", "", old.names)

neospectra.nir <- neospectra.nir %>%
  rename_with(~new.names, all_of(old.names))

spectral.column.names <- new.names

neospectra.nir[1:5, 1:5]


## ----neospectra_join, eval=TRUE-------------------------------------------------------------------------------------------
neospectra <- left_join(neospectra.site,
                        neospectra.soil,
                        by = "id.sample_local_c") %>%
  left_join(., neospectra.nir, by = "id.sample_local_c") %>%
  filter(!is.na(oc_usda.c729_w.pct))

neospectra[1:5, 1:5]


## ----neospectra_prep, eval=TRUE-------------------------------------------------------------------------------------------
neospectra.nir.prep <- neospectra %>%
  select(all_of(spectral.column.names)) %>%
  as.matrix() %>%
  savitzkyGolay(X = ., p = 2, w = 11, m = 1, delta.wav = 2) %>%
  prospectr::standardNormalVariate(X = .) %>%
  as_tibble() %>%
  bind_cols({neospectra %>%
      select(id.sample_local_c,
             location.country_iso.3166_txt,
             oc_usda.c729_w.pct)}, .)

# Quick visualization of preprocessed spectra
set.seed(42)
neospectra.nir.prep %>%
  sample_n(100) %>%
  pivot_longer(any_of(spectral.column.names),
               names_to = "wavelength",
               values_to = "reflectance") %>%
  mutate(wavelength = as.numeric(wavelength),
         reflectance = as.numeric(reflectance)) %>%
  ggplot(data = .) +
  geom_line(aes(x = wavelength, y = reflectance,
                group = id.sample_local_c),
            alpha = 0.25, linewidth = 0.25) +
  theme_light()


## ----neospectra_split, eval=TRUE------------------------------------------------------------------------------------------
neospectra.train <- neospectra.nir.prep %>%
  filter(location.country_iso.3166_txt == "USA")

neospectra.test <- neospectra.nir.prep %>%
  filter(location.country_iso.3166_txt != "USA")

nrow(neospectra.train)
nrow(neospectra.test)


## ----compression, eval=TRUE-----------------------------------------------------------------------------------------------
# Define and fit the PCA recipe on training data
pca.recipe <- recipe(neospectra.train) %>%
  update_role(id.sample_local_c, new_role = "id variable") %>%
  update_role(location.country_iso.3166_txt, new_role = "id variable") %>%
  update_role(oc_usda.c729_w.pct, new_role = "outcome") %>%
  update_role(any_of(spectral.column.names), new_role = "predictor") %>%
  step_normalize(all_predictors()) %>%
  step_pca(all_predictors(), threshold = 0.99, id = "pca")

pca.prep <- prep(pca.recipe, training = neospectra.train)

# How many components were retained?
pca.prep

# Explained variance per component from tidy()
pca.expvar <- tidy(pca.prep, id = "pca", type = "variance") %>%
  filter(terms == "percent variance") %>%
  pull(value)

# Extract training and test PC scores via bake()
train.pca.scores <- bake(pca.prep, new_data = neospectra.train)
test.pca.scores  <- bake(pca.prep, new_data = neospectra.test)

# PC column names produced by step_pca
pc.names <- train.pca.scores %>%
  select(starts_with("PC")) %>%
  names()

# PCA plot coloured by log1p(SOC)
p.pca <- ggplot(train.pca.scores) +
  geom_point(aes(x = PC01, y = PC02,
                 color = log1p(oc_usda.c729_w.pct)),
             alpha = 0.5, size = 0.5) +
  scale_colour_gradient(low = "gold", high = "darkred") +
  labs(title = "Neospectra PCA compression",
       x = paste0("PC1 (", round(pca.expvar[1], 2), "%)"),
       y = paste0("PC2 (", round(pca.expvar[2], 2), "%)")) +
  theme_light() +
  theme(legend.position = "bottom")

# Adding test samples to PCA plot
p.pca +
  geom_point(data = test.pca.scores,
             aes(x = PC01, y = PC02),
             size = 0.75) +
  labs(subtitle = "Black dots represent testing points")


## ----backtransform, message=FALSE, error=FALSE----------------------------------------------------------------------------
# Extract loadings, centering means, and scaling SDs from the prepped recipe
n.pad <- nchar(length(pc.names))

pca.loadings <- tidy(pca.prep, id = "pca", type = "coef") %>%
  mutate(component = paste0("PC", formatC(as.integer(gsub("PC", "", component)),
                                          width = n.pad, flag = "0"))) %>%
  pivot_wider(names_from = component, values_from = value) %>%
  select(terms, all_of(pc.names)) %>%
  column_to_rownames("terms") %>%
  as.matrix()

norm.stats <- tidy(pca.prep, number = 1) %>%  # step_normalize is step 1
  select(terms, statistic, value) %>%
  pivot_wider(names_from = statistic, values_from = value)

train.center <- setNames(norm.stats$mean, norm.stats$terms)
train.scale  <- setNames(norm.stats$sd,   norm.stats$terms)

# Get PC scores as matrices (columns ordered to match loadings)
train.scores.mat <- train.pca.scores %>%
  select(all_of(pc.names)) %>%
  as.matrix()

test.scores.mat <- test.pca.scores %>%
  select(all_of(pc.names)) %>%
  as.matrix()

# Back-transform: rotate scores to spectral space, then rescale and recenter
train.bt <- train.scores.mat %*% t(pca.loadings)
test.bt  <- test.scores.mat  %*% t(pca.loadings)

train.bt <- sweep(train.bt, MARGIN = 2, FUN = "*", STATS = train.scale)
test.bt  <- sweep(test.bt,  MARGIN = 2, FUN = "*", STATS = train.scale)

train.bt <- sweep(train.bt, MARGIN = 2, FUN = "+", STATS = train.center)
test.bt  <- sweep(test.bt,  MARGIN = 2, FUN = "+", STATS = train.center)

# Original spectra matrices (rows ordered to match scores)
neospectra.train.spectra <- neospectra.train %>%
  select(any_of(spectral.column.names)) %>%
  as.matrix()

neospectra.test.spectra <- neospectra.test %>%
  select(any_of(spectral.column.names)) %>%
  as.matrix()

# Q-statistic: sum of squared differences between original and back-transformed
test.q.stats <- apply((neospectra.test.spectra - test.bt)^2,
                      MARGIN = 1, sum)


## ----q_critical, message=FALSE, error=FALSE-------------------------------------------------------------------------------
# Critical value from the training set (1% significance level)
E     <- cov(neospectra.train.spectra - train.bt)
teta1 <- sum(diag(E))^1
teta2 <- sum(diag(E))^2
teta3 <- sum(diag(E))^3
h0    <- 1 - ((2*teta1*teta3) / (3*teta2^2))
Ca    <- 2.57 # 1% significance level
Qa    <- teta1*(1 - (teta2*h0*((1-h0)/teta1^2)) +
                  ((sqrt(Ca*(2*teta2*h0^2)))/teta1))^(1/h0)
Qa


## ----q_flag, message=FALSE, error=FALSE-----------------------------------------------------------------------------------
# Flag test samples and visualize on PCA plot
test.pca.scores <- test.pca.scores %>%
  mutate(q_stats = test.q.stats,
         represented = q_stats <= Qa)

p.pca +
  geom_point(data = test.pca.scores,
             aes(x = PC01, y = PC02, fill = represented),
             shape = 21, size = 1.5) +
  labs(fill = "Represented?")


## ----log1p----------------------------------------------------------------------------------------------------------------
# Distribution of raw SOC values in the training set
hist(train.pca.scores$oc_usda.c729_w.pct, breaks = 50,
     main = "SOC distribution (raw)", xlab = "SOC (wt%)")

train.data <- train.pca.scores %>%
  mutate(oc_usda.c729_w.pct = log1p(oc_usda.c729_w.pct))

test.data <- test.pca.scores %>%
  mutate(oc_usda.c729_w.pct = log1p(oc_usda.c729_w.pct))

# Distribution after transformation
hist(train.data$oc_usda.c729_w.pct, breaks = 50,
     main = "SOC distribution (log1p)", xlab = "log1p(SOC)")


## ----recipe---------------------------------------------------------------------------------------------------------------
base.recipe <- train.data %>%
  recipe() %>%
  update_role(id.sample_local_c, new_role = "id variable") %>%
  update_role(location.country_iso.3166_txt, new_role = "id variable") %>%
  update_role(oc_usda.c729_w.pct, new_role = "outcome") %>%
  update_role(all_of(pc.names), new_role = "predictor") %>%
  step_normalize(all_predictors())

base.recipe


## ----resamples------------------------------------------------------------------------------------------------------------
set.seed(42)
train.folds <- vfold_cv(train.data, v = 10)

ctrl.resamples <- control_resamples(save_pred = TRUE,
                                    save_workflow = TRUE,
                                    extract = function(x) x)

ctrl.grid <- control_grid(save_pred = TRUE,
                          save_workflow = TRUE)


## ----model_specs----------------------------------------------------------------------------------------------------------
# Linear regression (no hyperparameters to tune)
spec.lm <- linear_reg() %>%
  set_engine("lm") %>%
  set_mode("regression")

# Random forest
spec.rf <- rand_forest(mtry  = tune(),
                       min_n = tune(),
                       trees = 500) %>%
  set_engine("ranger") %>%
  set_mode("regression")

# Cubist
spec.cubist <- cubist_rules(committees = tune(),
                            neighbors  = 0) %>%
  set_engine("Cubist") %>%
  set_mode("regression")

# Support vector machine (radial basis function kernel)
spec.svm <- svm_rbf(cost      = tune(),
                    rbf_sigma = tune()) %>%
  set_engine("kernlab") %>%
  set_mode("regression")


## ----grids----------------------------------------------------------------------------------------------------------------
grid.rf <- grid_regular(
  mtry(range = c(2L, floor(length(pc.names) / 3))),
  min_n(range = c(2L, 20L)),
  levels = 4)

grid.cubist <- grid_regular(
  committees(range = c(1L, 20L)),
  levels = 3)

grid.svm <- grid_regular(
  cost(range      = c(-1, 2)),
  rbf_sigma(range = c(-4, -1)),
  levels = 4)


## ----workflows_tuning, message=FALSE, warning=FALSE-----------------------------------------------------------------------
# Linear regression — fit_resamples (no HPO needed)
wf.lm <- workflow() %>%
  add_recipe(base.recipe) %>%
  add_model(spec.lm)

set.seed(42)
res.lm <- wf.lm %>%
  fit_resamples(resamples = train.folds,
                control = ctrl.resamples)

# Random forest — tune_grid
wf.rf <- workflow() %>%
  add_recipe(base.recipe) %>%
  add_model(spec.rf)

set.seed(42)
res.rf <- wf.rf %>%
  tune_grid(resamples = train.folds,
            grid = grid.rf,
            control = ctrl.grid)

# Cubist — tune_grid
wf.cubist <- workflow() %>%
  add_recipe(base.recipe) %>%
  add_model(spec.cubist)

set.seed(42)
res.cubist <- wf.cubist %>%
  tune_grid(resamples = train.folds,
            grid = grid.cubist,
            control = ctrl.grid)

# SVM — tune_grid
wf.svm <- workflow() %>%
  add_recipe(base.recipe) %>%
  add_model(spec.svm)

set.seed(42)
res.svm <- wf.svm %>%
  tune_grid(resamples = train.folds,
            grid = grid.svm,
            control = ctrl.grid)


## ----select_best, message=FALSE, warning=FALSE----------------------------------------------------------------------------
best.rf <- select_by_one_std_err(res.rf,
                                 mtry, min_n,
                                 metric = "rmse")

best.rf

best.cubist <- select_by_one_std_err(res.cubist,
                                     committees,
                                     metric = "rmse")

best.cubist

best.svm <- select_by_one_std_err(res.svm,
                                  cost, rbf_sigma,
                                  metric = "rmse")

best.svm


## ----compare_models, message=FALSE, warning=FALSE-------------------------------------------------------------------------
cv.metrics <- bind_rows(
  {res.lm %>%
    collect_metrics() %>%
    filter(.metric == "rmse") %>%
    summarise(model = "Linear Regression", mean_rmse = mean(mean))},
  {res.rf %>%
    collect_metrics() %>%
    filter(.metric == "rmse",
           mtry  == best.rf$mtry,
           min_n == best.rf$min_n) %>%
    summarise(model = "Random Forest", mean_rmse = mean(mean))},
  {res.cubist %>%
    collect_metrics() %>%
    filter(.metric == "rmse",
           committees == best.cubist$committees) %>%
    summarise(model = "Cubist", mean_rmse = mean(mean))},
  {res.svm %>%
    collect_metrics() %>%
    filter(.metric == "rmse",
           cost      == best.svm$cost,
           rbf_sigma == best.svm$rbf_sigma) %>%
    summarise(model = "Support Vector Machine", mean_rmse = mean(mean))})

cv.metrics %>% arrange(mean_rmse)


## ----final_fits, message=FALSE, warning=FALSE-----------------------------------------------------------------------------
final.fit.lm <- wf.lm %>%
  fit(data = train.data)

final.fit.rf <- wf.rf %>%
  finalize_workflow(best.rf) %>%
  fit(data = train.data)

final.fit.cubist <- wf.cubist %>%
  finalize_workflow(best.cubist) %>%
  fit(data = train.data)

final.fit.svm <- wf.svm %>%
  finalize_workflow(best.svm) %>%
  fit(data = train.data)


## ----best_overall---------------------------------------------------------------------------------------------------------
# Replace with the model that had the lowest CV RMSE in your run
best.final.fit <- final.fit.svm
best.res <- res.svm
best.params <- best.svm


## ----test_predictions, message=FALSE, warning=FALSE-----------------------------------------------------------------------
test.results <- predict(best.final.fit, new_data = test.data) %>%
  bind_cols(test.data %>%
              select(id.sample_local_c,
                     location.country_iso.3166_txt,
                     oc_usda.c729_w.pct)) %>%
  rename(predicted = .pred,
         observed  = oc_usda.c729_w.pct) %>%
  mutate(observed  = expm1(observed),
         predicted = expm1(predicted))

test.results[1:5,]


## ----performance, message=FALSE, warning=FALSE----------------------------------------------------------------------------
test.performance <- test.results %>%
  summarise(n    = n(),
            rmse = rmse_vec(truth = observed, estimate = predicted),
            bias = msd_vec(truth = observed, estimate = predicted),
            rsq  = rsq_trad_vec(truth = observed, estimate = predicted),
            ccc  = ccc_vec(truth = observed, estimate = predicted, bias = TRUE),
            rpd = rpiq_vec(truth = observed, estimate = predicted),
            rpiq = rpiq_vec(truth = observed, estimate = predicted))

test.performance


## ----scatterplot, message=FALSE, warning=FALSE----------------------------------------------------------------------------
performance.annotation <- paste0("Rsq = ",  round(test.performance[[1, "rsq"]],  2),
                                 "\nRMSE = ", round(test.performance[[1, "rmse"]], 2), " wt%")

p.final <- ggplot(test.results) +
  geom_point(aes(x = observed, y = predicted)) +
  geom_abline(intercept = 0, slope = 1) +
  annotate(geom = "text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2,
            label = performance.annotation, size = 3) +
  labs(title = "Soil Organic Carbon (wt%) — test prediction",
       x = "Observed", y = "Predicted") +
  theme_light()

r.max <- max(layer_scales(p.final)$x$range$range)
r.min <- min(layer_scales(p.final)$x$range$range)
s.max <- max(layer_scales(p.final)$y$range$range)
s.min <- min(layer_scales(p.final)$y$range$range)
t.max <- round(max(r.max, s.max), 1)
t.min <- round(min(r.min, s.min), 1)

p.final + coord_equal(xlim = c(t.min, t.max), ylim = c(t.min, t.max))


## ----conformal, message=FALSE, warning=FALSE------------------------------------------------------------------------------
# Finalized workflow for the best model
final.wf.conformal <- wf.svm %>%
  finalize_workflow(best.svm)

# Resample with extract to populate .extracts for int_conformal_cv()
ctrl.conformal <- control_resamples(save_pred = TRUE,
                                    extract = function(x) x)

set.seed(42)
res.conformal <- final.wf.conformal %>%
  fit_resamples(resamples = train.folds,
                control   = ctrl.conformal)

# Build the conformal object from the CV resampling results of the best model
conformal.object <- int_conformal_cv(res.conformal)
conformal.object


## ----conformal_predict, message=FALSE, warning=FALSE----------------------------------------------------------------------
# Prediction intervals on the log1p scale, then back-transformed
test.intervals <- predict(conformal.object,
                          new_data = test.data,
                          level    = 0.95) %>%
  bind_cols(test.data %>%
              select(id.sample_local_c,
                     oc_usda.c729_w.pct)) %>%
  rename(observed = oc_usda.c729_w.pct) %>%
  mutate(observed    = expm1(observed),
         .pred       = expm1(.pred),
         .pred_lower = expm1(.pred_lower),
         .pred_upper = expm1(.pred_upper),
         covered     = observed >= .pred_lower & observed <= .pred_upper)

test.intervals[1:5,]


## ----conformal_coverage, message=FALSE, warning=FALSE---------------------------------------------------------------------
# Coverage statistics
test.intervals %>%
  count(covered) %>%
  mutate(percentage = n / sum(n) * 100)

# Mean interval width
test.intervals %>%
  mutate(interval_width = .pred_upper - .pred_lower) %>%
  summarise(mean_width   = mean(interval_width),
            median_width = median(interval_width))


## ----conformal_plot, message=FALSE, warning=FALSE-------------------------------------------------------------------------
conformal.annotation <- paste0(
  "Rsq = ",       round(test.performance[[1, "rsq"]],  2),
  "\nRMSE = ",    round(test.performance[[1, "rmse"]], 2), " wt%",
  "\nCoverage = ", round(mean(test.intervals$covered) * 100, 1), "%")

p.conformal <- ggplot(test.intervals) +
  geom_pointrange(aes(x = observed,
                      y = .pred,
                      ymin = .pred_lower,
                      ymax = .pred_upper,
                      color = covered),
                  alpha = 0.6, linewidth = 0.4) +
  geom_abline(intercept = 0, slope = 1) +
  annotate(geom = "text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2,
            label = conformal.annotation, size = 3) +
  scale_color_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato")) +
  labs(title = "Soil Organic Carbon (wt%) and prediction intervals (PI95%)",
       x = "Observed", y = "Predicted",
       color = "PI95% covered:") +
  theme_light() +
  theme(legend.position = "bottom")

r.max <- max(layer_scales(p.conformal)$x$range$range)
r.min <- min(layer_scales(p.conformal)$x$range$range)
s.max <- max(layer_scales(p.conformal)$y$range$range)
s.min <- min(layer_scales(p.conformal)$y$range$range)
t.max <- round(max(r.max, s.max), 1)
t.min <- round(min(r.min, s.min), 1)

p.conformal + coord_equal(xlim = c(t.min, t.max), ylim = c(t.min, t.max))

