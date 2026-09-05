# run_all.R -----------------------------------------------------------------
# Reproduce every table and figure from the raw sheet. Run from the repo root:
#   Rscript run_all.R

scripts <- c(
  "R/01_tidy_data.R",
  "R/02_anova_duncan.R",
  "R/03_mixed_models.R",
  "R/04_flowering.R",
  "R/05_figures.R",
  "R/06_ode_model.R"
)

for (s in scripts) {
  message("\n--- ", s, " ---")
  source(s, echo = FALSE)
}

message("\nDone. Tables in outputs/tables, figures in outputs/figures.")
