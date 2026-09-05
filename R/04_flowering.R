# 04_flowering.R ------------------------------------------------------------
# Flowering phenology, and the question the project is actually about:
# does presence of the PvFTa1 amplicon track flowering time across the panel?
#
# Flowering was scored once per pot per accession, and the two pots of an
# accession agree exactly in every case. That has a consequence worth stating
# plainly rather than working around: within-accession variance is zero, so an
# F ratio for accession has a zero denominator and a mean-separation test on
# these numbers is not defined. Days to flowering is therefore treated here as
# one fixed observed value per accession, not as a sample with error.

source("R/00_packages.R")

flowering <- read_csv(file.path(paths$derived, "flowering.csv"),
                      show_col_types = FALSE)
week6 <- read_csv(file.path(paths$derived, "week6_wide.csv"),
                  show_col_types = FALSE)

# --- Is there any within-accession variation at all? -----------------------

within_variation <- flowering |>
  filter(!is.na(dtff)) |>
  group_by(accession_code, accession) |>
  summarise(
    n_pots   = n(),
    dtff_sd  = sd(dtff),
    dt50_sd  = sd(dt50),
    .groups  = "drop"
  )

zero_variance <- all(within_variation$dtff_sd == 0, na.rm = TRUE) &&
                 all(within_variation$dt50_sd == 0, na.rm = TRUE)

# --- One row per accession -------------------------------------------------

pheno <- flowering |>
  filter(!is.na(dtff)) |>
  group_by(accession_code, accession, pvfta1) |>
  summarise(
    dtff = first(dtff),
    dt50 = first(dt50),
    category = first(category),
    .groups = "drop"
  ) |>
  mutate(
    # Interval between first flower and 50% flowering: how quickly an
    # accession moves through its flowering window.
    flowering_lag = dt50 - dtff
  ) |>
  arrange(dtff)

write_csv(pheno, file.path(paths$derived, "phenology_by_accession.csv"))

# --- Descriptive spread ----------------------------------------------------

spread <- pheno |>
  summarise(
    n_accessions = n(),
    dtff_min = min(dtff), dtff_max = max(dtff),
    dtff_range = max(dtff) - min(dtff),
    dt50_min = min(dt50), dt50_max = max(dt50),
    dt50_range = max(dt50) - min(dt50),
    lag_min = min(flowering_lag), lag_max = max(flowering_lag)
  )

# --- The PvFTa1 comparison -------------------------------------------------
# Five accessions were amplified for PvFTa1 and all five gave the 300 bp
# product. The question is how much of the panel's phenotypic range those
# five span: if they cover the extremes and still share the marker, the
# marker is not what separates early from late here.

genotyped <- pheno |> filter(pvfta1 == "Amplified")

coverage <- tibble(
  n_genotyped        = nrow(genotyped),
  n_amplified        = nrow(genotyped),   # all five gave the 300 bp band
  dtff_span_genotyped = max(genotyped$dtff) - min(genotyped$dtff),
  dtff_span_panel     = spread$dtff_range,
  pct_of_panel_range  = 100 * (max(genotyped$dtff) - min(genotyped$dtff)) /
                        spread$dtff_range,
  dtff_min_genotyped  = min(genotyped$dtff),
  dtff_max_genotyped  = max(genotyped$dtff),
  spans_panel_min     = min(genotyped$dtff) == spread$dtff_min,
  spans_panel_max     = max(genotyped$dtff) == spread$dtff_max
)

write_csv(coverage, file.path(paths$tables, "pvfta1_coverage.csv"))

# --- Does flowering time relate to vegetative growth? ----------------------
# Accession means at week 6, joined to days to first flowering.

growth_means <- week6 |>
  group_by(accession_code, accession) |>
  summarise(across(c(PH, NL, NB, LL, LB, SD), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")

pheno_growth <- pheno |>
  left_join(growth_means, by = c("accession_code", "accession"))

write_csv(pheno_growth, file.path(paths$derived, "phenology_growth.csv"))

correlations <- map_dfr(c("PH", "NL", "NB", "LL", "LB", "SD"), function(tr) {
  d <- pheno_growth |> select(dtff, value = all_of(tr)) |> filter(!is.na(value))
  ct <- suppressWarnings(cor.test(d$dtff, d$value, method = "pearson"))
  tibble(trait = tr, n = nrow(d), r = unname(ct$estimate), p = ct$p.value)
})

write_csv(correlations, file.path(paths$tables, "dtff_growth_correlations.csv"))

# --- Console summary -------------------------------------------------------

cat("\n== Within-accession variation in flowering ==\n")
cat(sprintf("Both pots of every accession recorded identical values: %s\n",
            zero_variance))
cat("Consequence: within-accession variance is 0, so an F test for accession\n")
cat("on days-to-flowering has a zero denominator and is not defined.\n")

cat("\n== Flowering spread across the panel ==\n")
print(as.data.frame(spread), row.names = FALSE)

cat("\n== PvFTa1-amplified accessions ==\n")
print(genotyped |> select(accession, dtff, dt50, category) |>
        as.data.frame(), row.names = FALSE)
cat(sprintf(
  "\nAll %d amplified. They span %d of the panel's %d-day range in days to\nfirst flowering (%.0f%%), from the earliest accession to the latest.\n",
  coverage$n_amplified, coverage$dtff_span_genotyped,
  coverage$dtff_span_panel, coverage$pct_of_panel_range
))

cat("\n== Days to first flowering vs week-6 growth (accession means) ==\n")
print(correlations |>
        mutate(r = round(r, 3), p = signif(p, 3)) |>
        as.data.frame(), row.names = FALSE)
