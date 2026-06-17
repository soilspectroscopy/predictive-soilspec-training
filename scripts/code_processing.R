## ----working_directory-------------------------------------------------------------
my.wd <- "~/projects/local-files/soilspec_training"
# Or within an RStudio project
# my.wd <- getwd()


## ----setup, message=FALSE, warning=FALSE-------------------------------------------
library("tidyverse")
library("asdreader")
library("opusreader2")
library("prospectr")


## ----asdreader, message=FALSE, warning=FALSE---------------------------------------
# Downloading an .asd file
visnir.spectra.url <- "https://github.com/soilspectroscopy/ossl-models/raw/main/sample-data/101453MD01.asd"

visnir.spectra.path <- file.path(my.wd, "file1.asd")

download.file(url = visnir.spectra.url,
              destfile = visnir.spectra.path,
              mode = "wb")

# Reading asd file
visnir.spectra <- asdreader::get_spectra(visnir.spectra.path)

# Inspecting the file
class(visnir.spectra)
dim(visnir.spectra)
visnir.spectra[1,1:5]

# Spectral range
range(as.numeric(colnames(visnir.spectra)))


## ----opusreader2, message=FALSE, warning=FALSE-------------------------------------
# Downloading an .0 file
mir.spectra.url <- "https://github.com/soilspectroscopy/ossl-models/raw/main/sample-data/235157XS01.0"

mir.spectra.path <- file.path(my.wd, "file2.0")

download.file(url = mir.spectra.url,
              destfile = mir.spectra.path,
              mode = "wb")

# Reading opus file
mir.spectra <- opusreader2::read_opus_single(dsn = mir.spectra.path)

# Inspecting the file
class(mir.spectra)
names(mir.spectra)

# Spectra is stored in file$ab$data
class(mir.spectra$ab$data)
dim(mir.spectra$ab$data)

# Spectral range
range(as.numeric(colnames(mir.spectra$ab$data)))


## ----read_csv, message=FALSE, warning=FALSE----------------------------------------
# Downloading a csv output from Neospectra
nir.spectra.url <- "https://github.com/soilspectroscopy/ossl-models/raw/main/sample-data/sample_neospectra_data.csv"

nir.spectra.path <- file.path(my.wd, "file3.csv")

download.file(url = nir.spectra.url,
              destfile = nir.spectra.path,
              mode = "wb")

# Reading csv file
nir.spectra <- readr::read_csv(nir.spectra.path)

# Inspecting the file
class(nir.spectra)
nir.spectra[1:5,1:5]

# Spectral range after removing first column
range(as.numeric(colnames(nir.spectra[,-1])))


## ----nir_reflectance---------------------------------------------------------------
# Original data
nir.spectra[1:5,1:5]

# Spectra column names
spectra.column.names <- nir.spectra %>%
  select(-sample_id) %>%
  names()

# Transforming reflectance (%) to reflectance factor (decimal, 0-1)
# Also, rounding to 5 decimal places of precision
nir.spectra.rf <- nir.spectra %>%
  mutate(across(all_of(spectra.column.names), ~round(.x/100, 5)))

nir.spectra.rf[1:5,1:5]


## ----visualization-----------------------------------------------------------------
## Pivot to long format
nir.spectra.rf.long <- nir.spectra.rf %>%
  pivot_longer(all_of(spectra.column.names),
               names_to = "wavelength",
               values_to = "reflectance") %>%
  mutate(wavelength = as.numeric(wavelength),
         reflectance = as.numeric(reflectance))

head(nir.spectra.rf.long)

## Visualization
ggplot(data = nir.spectra.rf.long) +
  geom_line(aes(x = wavelength, y = reflectance,
                group = sample_id),
            alpha = 0.5, linewidth = 0.5) +
  theme_light()


## ----interpolation_columns---------------------------------------------------------
# Old columns, reversed, and as numeric
old.wavelength <- as.numeric(rev(spectra.column.names))
head(old.wavelength)

# New columns, increasing order, spaced 2 nm
new.wavelength <- seq(1350, 2550, by = 2)
head(new.wavelength)


## ----interpolation-----------------------------------------------------------------
# Selecting old spectra in increasing order
# Parse to matrix (input of prospectr::resample)
# Resample
# Parse to tibble
# Bind to original sample ids
nir.spectra.rf.int <- nir.spectra.rf %>%
  select(all_of(rev(spectra.column.names))) %>%
  as.matrix() %>%
  prospectr::resample(X = .,
                      wav = old.wavelength,
                      new.wav = new.wavelength,
                      interpol = "spline") %>%
  as_tibble() %>%
  bind_cols({nir.spectra.rf %>%
      select(sample_id)}, .)

nir.spectra.rf.int[1:5,1:5]


## ----visualization_resample--------------------------------------------------------
new.spectra.column.names <- as.character(new.wavelength)

nir.spectra.rf.int %>%
  pivot_longer(all_of(new.spectra.column.names),
               names_to = "wavelength",
               values_to = "reflectance") %>%
  mutate(wavelength = as.numeric(wavelength),
         reflectance = as.numeric(reflectance)) %>%
  ggplot(data = .) +
  geom_line(aes(x = wavelength, y = reflectance, group = sample_id),
            alpha = 0.5, linewidth = 0.5) +
  theme_light()


## ----sg----------------------------------------------------------------------------
# Select spectra columns
# Parse to matrix (input of prospectr::savitzkyGolay)
# Apply preprocessing
# Parse to tibble
# Bind the spectra to id column
nir.spectra.sg <- nir.spectra.rf.int %>%
  select(all_of(new.spectra.column.names)) %>%
  as.matrix() %>%
  savitzkyGolay(X = ., p = 2, w = 11, m = 1, delta.wav = 2) %>%
  as_tibble() %>%
  bind_cols({nir.spectra.rf %>%
      select(sample_id)}, .)

nir.spectra.sg[1:5,1:5]


## ----sg_visualization--------------------------------------------------------------
nir.spectra.sg %>%
  pivot_longer(any_of(new.spectra.column.names),
               names_to = "wavelength",
               values_to = "reflectance") %>%
  mutate(wavelength = as.numeric(wavelength),
         reflectance = as.numeric(reflectance)) %>%
  ggplot(data = .) +
  geom_line(aes(x = wavelength, y = reflectance, group = sample_id),
            alpha = 0.5, linewidth = 0.5) +
  theme_light()


## ----snv---------------------------------------------------------------------------
# Select spectra columns
# Parse to matrix (input of prospectr::standardNormalVariate)
# Apply preprocessing
# Parse to tibble
# Bind the spectra to id column
nir.spectra.snv <- nir.spectra.rf.int %>%
  select(all_of(new.spectra.column.names)) %>%
  as.matrix() %>%
  prospectr::standardNormalVariate(X = .) %>%
  as_tibble() %>%
  bind_cols({nir.spectra.rf %>%
      select(sample_id)}, .)

nir.spectra.snv[1:5,1:5]


## ----snv_visualization-------------------------------------------------------------
nir.spectra.snv %>%
  pivot_longer(any_of(new.spectra.column.names),
               names_to = "wavelength",
               values_to = "reflectance") %>%
  mutate(wavelength = as.numeric(wavelength),
         reflectance = as.numeric(reflectance)) %>%
  ggplot(data = .) +
  geom_line(aes(x = wavelength, y = reflectance, group = sample_id),
            alpha = 0.5, linewidth = 0.5) +
  theme_light()

