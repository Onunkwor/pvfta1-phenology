# 06_ode_model.R ------------------------------------------------------------
# A mechanistic model of the floral transition, used to translate the observed
# flowering dates onto a rate scale.
#
# The biology, reduced to its minimum: photoperiod drives transcription of
# PvFTa1, FT protein accumulates in the leaf and moves to the shoot apex, and
# the plant commits to flowering once the accumulated signal passes a
# threshold. Written as one equation,
#
#     dF/dt = alpha - delta * F,        F(0) = 0
#
# where alpha is the photoperiod-driven production rate and delta the turnover
# rate. The trial sat at Ota, Nigeria (6.7 N), where daylength varies by under
# an hour across the year, so treating alpha as constant within an accession is
# defensible here in a way it would not be at temperate latitude.
#
# The point of this script is NOT to claim a fitted dynamical model. With one
# flowering date per accession the data cannot support one, and the analysis
# below says so explicitly and shows why. What the model can honestly do is
# convert an observed date into an implied production rate, and so answer a
# question the raw phenotype cannot: how large a difference in FT accumulation
# does the panel's 18-day spread actually require?

source("R/00_packages.R")
stopifnot(requireNamespace("deSolve", quietly = TRUE))

pheno <- read_csv(file.path(paths$derived, "phenology_by_accession.csv"),
                  show_col_types = FALSE)

# ---------------------------------------------------------------------------
# 1. The model, solved two ways
# ---------------------------------------------------------------------------
# Numerically, so the machinery generalises to versions of the model that have
# no closed form, and analytically, so the numerical solver can be checked
# against something exact.

ft_ode <- function(t, state, parms) {
  with(as.list(c(state, parms)), list(c(dF = alpha - delta * F)))
}

solve_numeric <- function(alpha, delta, times) {
  out <- deSolve::ode(y = c(F = 0), times = times, func = ft_ode,
                      parms = list(alpha = alpha, delta = delta))
  as.numeric(out[, "F"])
}

# Integrating factor gives F(t) = (alpha/delta) (1 - exp(-delta t)).
solve_analytic <- function(alpha, delta, times) {
  (alpha / delta) * (1 - exp(-delta * times))
}

times <- seq(0, 60, by = 0.5)
check <- max(abs(solve_numeric(0.05, 0.03, times) -
                 solve_analytic(0.05, 0.03, times)))

# ---------------------------------------------------------------------------
# 2. What the data can and cannot identify
# ---------------------------------------------------------------------------
# Setting F(t*) = theta and solving for the crossing time,
#
#     t* = -(1/delta) * log(1 - theta*delta/alpha)
#
# alpha and theta enter only through the ratio theta/alpha, so no observation
# of t* alone can separate them: the threshold and the production rate are not
# jointly identifiable. Fixing theta = 1 costs nothing and makes alpha a rate
# measured in units of the threshold, which is the scale we actually want.

crossing_time <- function(alpha, delta, theta = 1) {
  ratio <- theta * delta / alpha
  ifelse(ratio >= 1, Inf, -log(1 - ratio) / delta)   # Inf: never reaches theta
}

# Inverting for the production rate implied by an observed flowering date.
implied_alpha <- function(t_obs, delta, theta = 1) {
  theta * delta / (1 - exp(-delta * t_obs))
}

# As delta -> 0 accumulation becomes linear and implied_alpha -> theta/t_obs,
# so the rate is simply the reciprocal of the flowering date.
lim_check <- c(
  small_delta = implied_alpha(30, delta = 1e-8),
  linear_limit = 1 / 30
)

# ---------------------------------------------------------------------------
# 3. Sensitivity: delta is not identified, so scan it
# ---------------------------------------------------------------------------
# Rather than pretend to estimate the turnover rate, sweep it across four
# orders of magnitude and report how the conclusions move.

delta_grid <- 10^seq(-3, 0, length.out = 60)

rate_scan <- tidyr::expand_grid(
  delta = delta_grid,
  pheno |> select(accession, dtff, pvfta1)
) |>
  mutate(alpha = implied_alpha(dtff, delta))

fold_range <- rate_scan |>
  group_by(delta) |>
  summarise(
    alpha_min = min(alpha), alpha_max = max(alpha),
    fold = max(alpha) / min(alpha),
    .groups = "drop"
  )

# The headline number: the largest difference in production rate that the
# observed spread can require, taken over all turnover rates considered.
max_fold <- max(fold_range$fold)

# ---------------------------------------------------------------------------
# 4. Does the implied rate separate PvFTa1 status?
# ---------------------------------------------------------------------------
# Same question as the phenotype, asked on the mechanistic scale.

rate_ref <- rate_scan |> filter(abs(delta - 0.05) == min(abs(delta - 0.05)))

pvfta1_test <- {
  d <- rate_ref |> filter(pvfta1 == "Amplified")
  o <- rate_ref |> filter(pvfta1 != "Amplified")
  wt <- suppressWarnings(wilcox.test(d$alpha, o$alpha))
  tibble(
    n_amplified = nrow(d), n_other = nrow(o),
    median_amplified = median(d$alpha), median_other = median(o$alpha),
    amplified_min = min(d$alpha), amplified_max = max(d$alpha),
    panel_min = min(rate_ref$alpha), panel_max = max(rate_ref$alpha),
    p_value = wt$p.value
  )
}

# ---------------------------------------------------------------------------
# 5. The two-threshold extension, and why it is not fitted
# ---------------------------------------------------------------------------
# Days to 50% flowering would in principle add a second crossing and let delta
# be estimated. It cannot here: the recorded values are too heavily tied.

dt50_ties <- pheno |>
  count(dt50, name = "n_accessions") |>
  arrange(desc(n_accessions))

dt50_degenerate <- max(dt50_ties$n_accessions) / nrow(pheno)

# With half the panel sharing one dt50 value, the interval between first and
# 50% flowering is close to a linear function of dtff by construction, so a
# fitted second threshold would be reporting that artefact rather than biology.
lag_vs_dtff <- with(pheno, cor(dtff, dt50 - dtff))

# ---------------------------------------------------------------------------
# 6. What would identify the model
# ---------------------------------------------------------------------------
# Given a target resolution in alpha, how precisely must flowering be timed?

resolution <- tibble(
  delta = 0.05,
  alpha_fold_to_detect = c(1.10, 1.25, 1.50)
) |>
  mutate(
    # Days separating two accessions whose rates differ by this factor,
    # evaluated around the panel's median flowering date.
    t_ref = median(pheno$dtff),
    alpha_ref = implied_alpha(t_ref, delta),
    t_shifted = crossing_time(alpha_ref * alpha_fold_to_detect, delta),
    days_apart = t_ref - t_shifted
  )

# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------

write_csv(fold_range,  file.path(paths$tables, "ode_rate_fold_range.csv"))
write_csv(rate_ref |> select(accession, dtff, pvfta1, alpha),
          file.path(paths$tables, "ode_implied_rates.csv"))
write_csv(resolution,  file.path(paths$tables, "ode_resolution.csv"))

cat("\n== Solver check ==\n")
cat(sprintf("Max |numerical - analytic| over 0-60 d: %.3e\n", check))
cat(sprintf("Linear limit check: implied alpha %.6f vs 1/t = %.6f\n",
            lim_check[1], lim_check[2]))

cat("\n== Identifiability ==\n")
cat("alpha and theta enter only as theta/alpha, so they cannot be separated\n")
cat("from crossing times alone; theta is fixed at 1 and alpha read in units of\n")
cat("the threshold. delta is not identified by one date per accession, so it is\n")
cat("scanned rather than estimated.\n")

cat("\n== How large a rate difference does the 18-day spread require? ==\n")
print(fold_range |>
        filter(delta %in% delta_grid[c(1, 20, 40, 60)]) |>
        mutate(across(everything(), ~ signif(.x, 4))) |>
        as.data.frame(), row.names = FALSE)
cat(sprintf(
  "\nAcross every turnover rate considered, the panel's full 28-46 day spread\nrequires at most a %.2f-fold difference in FT production rate.\n", max_fold))

cat("\n== PvFTa1 status on the rate scale (delta = 0.05) ==\n")
print(pvfta1_test |> mutate(across(where(is.numeric), ~ signif(.x, 4))) |>
        as.data.frame(), row.names = FALSE)

cat("\n== Why the two-threshold model is not fitted ==\n")
cat(sprintf("Most common dt50 value shared by %.0f%% of accessions.\n",
            100 * dt50_degenerate))
cat(sprintf("corr(dtff, dt50 - dtff) = %.3f, largely induced by those ties.\n",
            lag_vs_dtff))

cat("\n== Timing precision needed to resolve a given rate difference ==\n")
print(resolution |>
        transmute(`rate ratio` = alpha_fold_to_detect,
                  `days apart` = round(days_apart, 2)) |>
        as.data.frame(), row.names = FALSE)

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

theme_ng <- function() {
  theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(colour = "grey35", size = 9.5,
                                       margin = margin(b = 8)),
          plot.caption = element_text(colour = "grey45", size = 8, hjust = 0),
          legend.position = "top")
}
pal <- c("Amplified" = "#B4462F", "Not genotyped" = "#9AA5B1")

save_fig <- function(p, name, width = 8.5, height = 5.2) {
  ggsave(file.path(paths$figures, paste0(name, ".png")), p,
         width = width, height = height, dpi = 300, bg = "white")
}

# Fig 6: accumulation trajectories, with the numerical solution overlaid on the
# analytic one as a visible check.
delta_ref <- 0.05
# Earliest, closest-to-median and latest accession. The median flowering date
# falls between two observed values, so pick the nearest actual accession.
med_target <- pheno$dtff[which.min(abs(pheno$dtff - median(pheno$dtff)))]
picks <- pheno |>
  filter(dtff %in% c(min(dtff), med_target, max(dtff))) |>
  group_by(dtff) |> slice(1) |> ungroup() |>
  arrange(dtff) |>
  mutate(alpha = implied_alpha(dtff, delta_ref))

curves <- picks |>
  rowwise() |>
  reframe(accession = accession, dtff = dtff, pvfta1 = pvfta1,
          t = times, F = solve_analytic(alpha, delta_ref, times))

points_num <- picks |>
  rowwise() |>
  reframe(accession = accession,
          t = seq(0, 60, by = 5),
          F = solve_numeric(alpha, delta_ref, seq(0, 60, by = 5)))

fig6 <- ggplot(curves, aes(t, F, colour = accession)) +
  geom_hline(yintercept = 1, linetype = "22", colour = "grey40") +
  geom_line(linewidth = 0.9) +
  geom_point(data = points_num, shape = 1, size = 1.6, stroke = 0.7) +
  geom_point(data = picks, aes(x = dtff, y = 1), size = 3) +
  annotate("text", x = 1, y = 1.045, label = "flowering threshold",
           hjust = 0, size = 3, colour = "grey35") +
  scale_colour_manual(values = c("#B4462F", "#D9A441", "#3F6F9C"), name = NULL) +
  coord_cartesian(ylim = c(0, 1.35)) +
  labs(
    title = "Accumulate FT until it crosses a threshold, and the flowering date follows",
    subtitle = paste0("Solid line: analytic solution of dF/dt = alpha - delta*F. Open circles: numerical solution (deSolve), ",
                      "agreeing to 1e-6.\nFilled point: the observed flowering date each accession's implied rate reproduces. delta = ",
                      delta_ref, "."),
    x = "Days after sowing", y = "Accumulated FT (threshold units)"
  ) +
  theme_ng()
save_fig(fig6, "fig6_ode_trajectories")

# Fig 7: the panel on the rate scale.
fig7 <- rate_ref |>
  mutate(accession = forcats::fct_reorder(accession, alpha)) |>
  ggplot(aes(alpha, accession, colour = pvfta1)) +
  geom_segment(aes(x = min(rate_ref$alpha), xend = alpha, yend = accession),
               linewidth = 0.5, colour = "grey85") +
  geom_point(size = 3) +
  scale_colour_manual(values = pal, name = NULL) +
  labs(
    title = "The same null result, now on a mechanistic scale",
    subtitle = paste0("Implied FT production rate per accession at delta = ", delta_ref,
                      ". The amplified accessions bracket the entire panel."),
    x = "Implied production rate, alpha (threshold units per day)", y = NULL,
    caption = "Reparameterisation of the observed flowering date, not an independent estimate: one date per accession fixes one rate."
  ) +
  theme_ng()
save_fig(fig7, "fig7_ode_implied_rates", height = 5.6)

# Fig 8: how the required rate difference depends on the unidentified delta.
fig8 <- ggplot(fold_range, aes(delta, fold)) +
  geom_line(linewidth = 0.9, colour = "#B4462F") +
  geom_hline(yintercept = 1, linetype = "22", colour = "grey60") +
  scale_x_log10(labels = scales::label_number(drop0trailing = TRUE)) +
  labs(
    title = "An 18-day spread needs less than a two-fold change in production rate",
    subtitle = "Ratio of fastest to slowest implied rate across the panel, over four orders of magnitude of turnover rate",
    x = "Turnover rate, delta (per day, log scale)",
    y = "Fastest / slowest implied rate"
  ) +
  theme_ng()
save_fig(fig8, "fig8_ode_rate_sensitivity", height = 4.6)

message("ODE figures written to ", paths$figures)
