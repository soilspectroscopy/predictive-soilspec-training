
library("knitr")

purl("processing.qmd", "scripts/code_processing.R")
purl("machinelearning.qmd", "scripts/code_machine_learning.R")
purl("chemometrics.qmd", "scripts/code_chemometrics.R")
