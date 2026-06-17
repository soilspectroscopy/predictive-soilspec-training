## ----libraries, message=FALSE, warning=FALSE---------------------------------------
library("tidyverse") # Regular data wrangling
library("mdatools") # Chemometrics
library("yardstick") # Additional performance metrics


## ----files, message=FALSE, warning=FALSE-------------------------------------------
## Internet configuration for downloading large datasets
options(timeout = 10000)

neospectra.soil <- read_csv("https://storage.googleapis.com/soilspec4gg-public/neospectra_soillab_v1.2.csv.gz")
dim(neospectra.soil)

neospectra.site <- read_csv("https://storage.googleapis.com/soilspec4gg-public/neospectra_soilsite_v1.2.csv.gz")
dim(neospectra.site)

neospectra.nir <- read_csv("https://storage.googleapis.com/soilspec4gg-public/neospectra_nir_v1.2.csv.gz")
dim(neospectra.nir)


## ----preparation, message=FALSE, warning=FALSE-------------------------------------
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

# Inspecting and renaming spectral column headers to numeric values
neospectra.nir %>%
  select(starts_with("scan_nir")) %>%
  names() %>%
  head()

old.names <- neospectra.nir %>%
  select(starts_with("scan_nir")) %>%
  names()

new.names <- gsub("scan_nir.|_ref", "", old.names)

neospectra.nir <- neospectra.nir %>%
  rename_with(~new.names, all_of(old.names))

spectral.column.names <- new.names

head(spectral.column.names)
tail(spectral.column.names)


## ----join, message=FALSE, warning=FALSE--------------------------------------------
# Joining data
neospectra <- left_join(neospectra.site,
                        neospectra.soil,
                        by = "id.sample_local_c") %>%
  left_join(., neospectra.nir, by = "id.sample_local_c")

# Filtering out samples without SOC values
neospectra <- neospectra %>%
  filter(!is.na(oc_usda.c729_w.pct))


## ----visualization, message=FALSE, warning=FALSE-----------------------------------
# Spectral visualization of a random subset
set.seed(42)
neospectra %>%
  sample_n(100) %>%
  pivot_longer(any_of(spectral.column.names),
               names_to = "wavelength",
               values_to = "reflectance") %>%
  mutate(wavelength = as.numeric(wavelength),
         reflectance = as.numeric(reflectance)) %>%
  ggplot(data = .) +
  geom_line(aes(x = wavelength, y = reflectance, group = id.sample_local_c),
            alpha = 0.25, linewidth = 0.25) +
  theme_light()


## ----preprocessing, message=FALSE, warning=FALSE-----------------------------------
# Preprocessing with SG 1st Der. and SNV
neospectra.prep <- neospectra %>%
  select(all_of(spectral.column.names)) %>%
  as.matrix() %>%
  prep.savgol(width = 11, porder = 1, dorder = 1) %>%
  prep.snv(.) %>%
  as_tibble() %>%
  bind_cols({neospectra %>%
      select(-all_of(spectral.column.names))}, .)

# Visualization of preprocessed spectra
set.seed(42)
neospectra.prep %>%
  sample_n(100) %>%
  pivot_longer(any_of(spectral.column.names),
               names_to = "wavelength",
               values_to = "reflectance") %>%
  mutate(wavelength = as.numeric(wavelength),
         reflectance = as.numeric(reflectance)) %>%
  ggplot(data = .) +
  geom_line(aes(x = wavelength, y = reflectance, group = id.sample_local_c),
            alpha = 0.25, linewidth = 0.25) +
  theme_light()


## ----split, message=FALSE, warning=FALSE-------------------------------------------
# Train split: USA samples
# Predictors and outcome are separated; log1p is applied to SOC to reduce skewness
neospectra.train <- neospectra.prep %>%
  filter(location.country_iso.3166_txt == "USA")

neospectra.train.predictors <- neospectra.train %>%
  select(all_of(spectral.column.names)) %>%
  as.matrix()

neospectra.train.outcome <- neospectra.train %>%
  select(oc_usda.c729_w.pct) %>%
  mutate(oc_usda.c729_w.pct = log1p(oc_usda.c729_w.pct)) %>%
  as.matrix()

# Test split: non-USA samples
neospectra.test <- neospectra.prep %>%
  filter(location.country_iso.3166_txt != "USA")

neospectra.test.predictors <- neospectra.test %>%
  select(all_of(spectral.column.names)) %>%
  as.matrix()

neospectra.test.outcome <- neospectra.test %>%
  select(oc_usda.c729_w.pct) %>%
  mutate(oc_usda.c729_w.pct = log1p(oc_usda.c729_w.pct)) %>%
  as.matrix()


## ----plsr, message=FALSE, warning=FALSE--------------------------------------------
set.seed(42)
pls.model.soc <- pls(x = neospectra.train.predictors,
                     y = neospectra.train.outcome,
                     ncomp = 20,
                     x.test = neospectra.test.predictors,
                     y.test = neospectra.test.outcome,
                     center = TRUE, scale = FALSE,
                     cv = 10, lim.type = "ddmoments", cv.scope = 'global',
                     info = "SOC prediction model")


## ----plsr_split, message=FALSE, warning=FALSE--------------------------------------
# Alternatively: fit calibration model first, then predict on test set
set.seed(42)
pls.model.soc.calibration <- pls(x = neospectra.train.predictors,
                                 y = neospectra.train.outcome,
                                 ncomp = 20,
                                 center = TRUE, scale = FALSE,
                                 cv = 10, lim.type = "ddmoments", cv.scope = 'global',
                                 info = "SOC prediction model")

pls.model.soc.predictions <- predict(pls.model.soc.calibration,
                                     x = neospectra.test.predictors,
                                     y = neospectra.test.outcome,
                                     cv = FALSE)


## ----model_summary-----------------------------------------------------------------
summary(pls.model.soc)


## ----model_viz---------------------------------------------------------------------
# Visualize representation limits, model coefficients, performance, and fit
plot(pls.model.soc, show.legend = TRUE)


## ----select_ncomp------------------------------------------------------------------
# Performance is more consistent across train, CV, and test with 11 components
pls.model.soc <- selectCompNum(pls.model.soc, 11)
summary(pls.model.soc)
plot(pls.model.soc, show.legend = FALSE)


## ----classification, message=FALSE, warning=FALSE----------------------------------
# Sample categorization based on Q and T2 limits
outlier.detection <- categorize(pls.model.soc.calibration,
                                pls.model.soc.predictions,
                                ncomp = 11)
head(outlier.detection)

plotResiduals(pls.model.soc.predictions,
              cgroup = outlier.detection,
              ncomp = 11)


## ----performance_components--------------------------------------------------------
plotRMSE(pls.model.soc, show.labels = TRUE)


## ----scatterplot-------------------------------------------------------------------
# Predictions scatterplot
plotPredictions(pls.model.soc, show.line = TRUE, pch = 20, cex = 0.25)
abline(a = 0, b = 1)


## ----residuals---------------------------------------------------------------------
# Inspection plot for outcome residuals
plotYResiduals(pls.model.soc, cex = 0.5, show.label = TRUE)


## ----vip---------------------------------------------------------------------------
# Variable Importance in Projection (VIP) scores
plotVIPScores(pls.model.soc)
vip <- vipscores(pls.model.soc, ncomp = 11)
head(vip)


## ----x_variance--------------------------------------------------------------------
# Cumulative variance retained from predictors (spectra)
plotXCumVariance(pls.model.soc, type = 'h', show.labels = TRUE, legend.position = 'bottomright')


## ----y_variance--------------------------------------------------------------------
# Cumulative variance retained from outcome
plotYCumVariance(pls.model.soc, type = 'b', show.labels = TRUE, legend.position = 'bottomright')


## ----x_scores----------------------------------------------------------------------
# Score plots for compressed spectra
plotXScores(pls.model.soc, comp = c(1, 2), show.legend = TRUE, cex = 0.5, legend.position = 'topleft')
# plotXScores(pls.model.soc, comp = c(1, 3), show.legend = TRUE, cex = 0.5, legend.position = 'bottomright')


## ----loadings----------------------------------------------------------------------
# Loadings plot for the first 5 components
plotXLoadings(pls.model.soc, comp = c(1, 2, 3, 4, 5), type = 'l')


## ----internal----------------------------------------------------------------------
# Exploring internal info from model object
pls.model.soc$T2lim
pls.model.soc$Qlim

pls.model.soc$res$cal
pls.model.soc$res$cal$xdecomp

pls.model.soc$res$test
pls.model.soc$res$test$xdecomp


## ----predictions-------------------------------------------------------------------
# Extracting test predictions
test.predictions <- pls.model.soc$res$test$y.pred
dim(test.predictions)

# Selecting predictions at component 11
test.predictions <- tibble(predicted = test.predictions[,11,1])

# Back-transforming predicted values from log1p to original scale
neospectra.test.results <- neospectra.test %>%
  select(-all_of(spectral.column.names)) %>%
  bind_cols(test.predictions) %>%
  rename(observed = oc_usda.c729_w.pct) %>%
  mutate(predicted = expm1(predicted))


## ----yardstick---------------------------------------------------------------------
# Calculating performance metrics
neospectra.test.performance <- neospectra.test.results %>%
  summarise(n = n(),
            rmse = rmse_vec(truth = observed, estimate = predicted),
            bias = msd_vec(truth = observed, estimate = predicted),
            rsq = rsq_trad_vec(truth = observed, estimate = predicted),
            ccc = ccc_vec(truth = observed, estimate = predicted, bias = TRUE),
            rpd = rpd_vec(truth = observed, estimate = predicted),
            rpiq = rpiq_vec(truth = observed, estimate = predicted))

neospectra.test.performance


## ----ggplot------------------------------------------------------------------------
# Final accuracy plot
performance.annotation <- paste0("Rsq = ", round(neospectra.test.performance[[1,"rsq"]], 2),
                                "\nRMSE = ", round(neospectra.test.performance[[1,"rmse"]], 2), " wt%")

p.final <- ggplot(neospectra.test.results) +
  geom_point(aes(x = observed, y = predicted)) +
  geom_abline(intercept = 0, slope = 1) +
  annotate(geom = "text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2,
            label = performance.annotation, size = 3) +
  labs(title = "Soil Organic Carbon (wt%) test prediction",
       x = "Observed", y = "Predicted") +
  theme_light() + theme(legend.position = "bottom")

# Ensuring a square plot with equal axis limits
r.max <- max(layer_scales(p.final)$x$range$range)
r.min <- min(layer_scales(p.final)$x$range$range)

s.max <- max(layer_scales(p.final)$y$range$range)
s.min <- min(layer_scales(p.final)$y$range$range)

t.max <- round(max(r.max, s.max), 1)
t.min <- round(min(r.min, s.min), 1)

p.final <- p.final + coord_equal(xlim = c(t.min, t.max), ylim = c(t.min, t.max))
p.final

