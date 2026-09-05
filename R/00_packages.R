# 00_packages.R -------------------------------------------------------------
# Package bootstrap. Run once; every other script sources this.

cran <- "https://cloud.r-project.org"

required <- c(
  "dplyr", "tidyr", "readr", "stringr", "purrr", "tibble", "forcats",
  "ggplot2", "scales", "agricolae", "lme4", "broom", "broom.mixed", "knitr"
)

missing <- setdiff(required, rownames(installed.packages()))
if (length(missing)) {
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = cran)
}

invisible(lapply(required, function(p) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}))

# Deterministic output for anything that resamples.
set.seed(20250905)

# Project-relative paths, so scripts run from the repo root.
paths <- list(
  raw     = "data/raw",
  derived = "data/derived",
  figures = "outputs/figures",
  tables  = "outputs/tables"
)
invisible(lapply(paths, dir.create, showWarnings = FALSE, recursive = TRUE))
