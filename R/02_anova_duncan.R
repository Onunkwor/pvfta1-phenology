# 02_anova_duncan.R ---------------------------------------------------------
# The classical analysis the thesis reports: one-way ANOVA across accessions
# for each week-6 growth trait, followed by Duncan's multiple range test at
# alpha = 0.05.
#
# Two models are fitted per trait:
#   naive     value ~ accession                 (what the thesis reports)
#   adjusted  value ~ accession + reading       (adds the recording session)
#
# The second is fitted because the two recording sessions are not
# interchangeable: session 2 runs systematically higher on every size trait
# and uses a different scale for stem diameter. Leaving that in the residual
# inflates the error term and costs power.

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

# --- Session offset --------------------------------------------------------
# Quantify the between-session difference before deciding to adjust for it.
session_offset <- map_dfr(names(traits), function(tr) {
  d <- week6 |>
    select(reading, value = all_of(tr)) |>
    filter(!is.na(value))
  if (n_distinct(d$reading) < 2) return(NULL)
  tt <- t.test(value ~ reading, data = d)
  tibble(
    trait      = tr,
    trait_label = traits[[tr]],
    mean_r1    = unname(tt$estimate[1]),
    mean_r2    = unname(tt$estimate[2]),
    ratio_r2_r1 = unname(tt$estimate[2] / tt$estimate[1]),
    p_value    = tt$p.value
  )
})

write_csv(session_offset, file.path(paths$tables, "session_offset.csv"))

# --- ANOVA + Duncan --------------------------------------------------------

fit_trait <- function(tr, adjust_reading) {
  d <- week6 |>
    select(accession, reading, value = all_of(tr)) |>
    filter(!is.na(value)) |>
    mutate(accession = factor(accession), reading = factor(reading))

  # Duncan's test needs residual df; a trait with one observation per
  # accession cannot supply any.
  if (nrow(d) - n_distinct(d$accession) < 2) return(NULL)

  form <- if (adjust_reading) value ~ accession + reading else value ~ accession
  fit  <- aov(form, data = d)

  aov_tbl <- broom::tidy(fit) |>
    mutate(trait = tr, trait_label = traits[[tr]],
           model = if (adjust_reading) "adjusted" else "naive")

  duncan <- agricolae::duncan.test(fit, "accession", group = TRUE)
  grp <- duncan$groups |>
    tibble::rownames_to_column("accession") |>
    rename(mean = 2, duncan_group = groups) |>
    mutate(trait = tr, trait_label = traits[[tr]],
           model = if (adjust_reading) "adjusted" else "naive")

  list(anova = aov_tbl, groups = grp)
}

results <- map(names(traits), function(tr) {
  list(
    naive    = fit_trait(tr, adjust_reading = FALSE),
    adjusted = fit_trait(tr, adjust_reading = TRUE)
  )
})
names(results) <- names(traits)

anova_table <- map_dfr(results, function(r) {
  bind_rows(r$naive$anova, r$adjusted$anova)
}) |>
  select(trait, trait_label, model, term, df, sumsq, meansq, statistic, p.value)

duncan_table <- map_dfr(results, function(r) {
  bind_rows(r$naive$groups, r$adjusted$groups)
}) |>
  select(trait, trait_label, model, accession, mean, duncan_group) |>
  arrange(trait, model, desc(mean))

write_csv(anova_table,  file.path(paths$tables, "anova_week6.csv"))
write_csv(duncan_table, file.path(paths$tables, "duncan_week6.csv"))

# --- Console summary -------------------------------------------------------

cat("\n== Recording-session offset (week 6) ==\n")
print(session_offset |>
        mutate(across(where(is.numeric), ~ round(.x, 3))) |>
        as.data.frame(), row.names = FALSE)

cat("\n== Accession effect on week-6 traits ==\n")
print(anova_table |>
        filter(term == "accession") |>
        mutate(p.value = signif(p.value, 3),
               statistic = round(statistic, 2)) |>
        select(trait_label, model, df, F = statistic, p = p.value) |>
        as.data.frame(), row.names = FALSE)
