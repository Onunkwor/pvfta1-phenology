# 03_mixed_models.R ---------------------------------------------------------
# Separate genotype from the layers of the screen-house layout.
#
# The design is a CRD: 20 accessions, 2 pots each, 2 plants kept per pot after
# thinning. So a week-6 measurement sits inside three nested sources of
# variation - accession, pot within accession, plant within pot - and the
# recording session cuts across all of them.
#
#   value ~ reading + (1 | accession) + (1 | accession:pot)
#
# reading is fixed because there are only two sessions and they differ
# systematically (see 02). accession is random here so its variance can be
# read on the same scale as the pot and residual terms, which is what
# repeatability needs.

source("R/00_packages.R")

week6 <- read_csv(file.path(paths$derived, "week6_wide.csv"),
                  show_col_types = FALSE)

traits <- c(
  PH = "Plant height (cm)",
  NL = "Number of leaves",
  NB = "Number of branches",
  LL = "Leaf length (cm)",
  LB = "Leaf breadth (cm)",
  SD = "Stem diameter"
)

fit_mixed <- function(tr) {
  d <- week6 |>
    select(accession, pot, reading, value = all_of(tr)) |>
    filter(!is.na(value)) |>
    mutate(across(c(accession, pot, reading), factor))

  fit <- suppressMessages(suppressWarnings(
    lme4::lmer(value ~ reading + (1 | accession) + (1 | accession:pot),
               data = d, REML = TRUE,
               control = lme4::lmerControl(check.conv.singular = "ignore"))
  ))

  vc <- as.data.frame(lme4::VarCorr(fit))
  get_var <- function(grp) {
    v <- vc$vcov[vc$grp == grp]
    if (length(v) == 0) 0 else v
  }

  v_acc <- get_var("accession")
  v_pot <- get_var("accession:pot")
  v_res <- get_var("Residual")
  total <- v_acc + v_pot + v_res

  tibble(
    trait        = tr,
    trait_label  = traits[[tr]],
    n            = nrow(d),
    var_accession = v_acc,
    var_pot       = v_pot,
    var_residual  = v_res,
    # Share of total variance attributable to accession identity: the
    # repeatability of a single plant measurement, and an upper bound on
    # broad-sense heritability under these conditions.
    repeatability = v_acc / total,
    pct_pot       = v_pot / total,
    pct_residual  = v_res / total,
    singular      = lme4::isSingular(fit)
  )
}

variance_components <- map_dfr(names(traits), fit_mixed)
write_csv(variance_components, file.path(paths$tables, "variance_components.csv"))

cat("\n== Variance partition, week-6 traits ==\n")
print(variance_components |>
        transmute(
          trait_label,
          n,
          accession = sprintf("%.1f%%", 100 * repeatability),
          pot       = sprintf("%.1f%%", 100 * pct_pot),
          residual  = sprintf("%.1f%%", 100 * pct_residual),
          singular
        ) |>
        as.data.frame(), row.names = FALSE)
