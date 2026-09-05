# 05_figures.R --------------------------------------------------------------
# Figures. Each one is meant to carry a single claim from the README.

source("R/00_packages.R")

growth_long <- read_csv(file.path(paths$derived, "growth_long.csv"),
                        show_col_types = FALSE)
week6   <- read_csv(file.path(paths$derived, "week6_wide.csv"),
                    show_col_types = FALSE)
pheno   <- read_csv(file.path(paths$derived, "phenology_by_accession.csv"),
                    show_col_types = FALSE)
duncan  <- read_csv(file.path(paths$tables, "duncan_week6.csv"),
                    show_col_types = FALSE)
vc      <- read_csv(file.path(paths$tables, "variance_components.csv"),
                    show_col_types = FALSE)

theme_ng <- function() {
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = "grey92"),
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(colour = "grey35", size = 9.5,
                                   margin = margin(b = 8)),
      plot.caption  = element_text(colour = "grey45", size = 8, hjust = 0),
      legend.position = "top",
      legend.title  = element_text(size = 9),
      strip.text    = element_text(face = "bold", size = 9)
    )
}

pal <- c("Amplified" = "#B4462F", "Not genotyped" = "#9AA5B1")

save_fig <- function(plot, name, width = 8, height = 5.5) {
  ggsave(file.path(paths$figures, paste0(name, ".png")),
         plot, width = width, height = height, dpi = 300, bg = "white")
  invisible(plot)
}

# --- Fig 1: flowering time across the panel, PvFTa1 status marked -----------
# The headline figure: the five accessions carrying the amplicon are not
# clustered at the early end, they run the length of the axis.

fig1 <- pheno |>
  mutate(accession = forcats::fct_reorder(accession, dtff)) |>
  ggplot(aes(x = dtff, y = accession, colour = pvfta1)) +
  geom_segment(aes(x = dtff, xend = dt50, yend = accession),
               linewidth = 1.4, alpha = 0.55) +
  geom_point(size = 2.6) +
  geom_point(aes(x = dt50), size = 2.6, shape = 21, fill = "white",
             stroke = 1.1) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_x_continuous(breaks = seq(28, 52, 4)) +
  labs(
    title = "PvFTa1 amplifies across the whole flowering range, not just the early end",
    subtitle = paste(
      "Filled point: days to first flowering. Open point: days to 50% flowering.",
      "All five genotyped accessions gave the 300 bp product."
    ),
    x = "Days after sowing", y = NULL,
    caption = "Screen-house trial, Bells University of Technology, Nov 2024 - Jul 2025. One value per accession; both pots agreed exactly."
  ) +
  theme_ng()

save_fig(fig1, "fig1_flowering_range", width = 8.5, height = 6)

# --- Fig 2: growth trajectories, weeks 2-6 ---------------------------------

fig2 <- growth_long |>
  filter(trait == "PH", !is.na(value)) |>
  group_by(accession, week, pvfta1) |>
  summarise(mean_ph = mean(value), .groups = "drop") |>
  ggplot(aes(x = week, y = mean_ph, group = accession, colour = pvfta1)) +
  geom_line(alpha = 0.75, linewidth = 0.7) +
  geom_point(size = 1.3) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_x_continuous(breaks = 2:6) +
  labs(
    title = "Plant height, weeks 2 to 6",
    subtitle = "Accession means pooled across pots, plants and both recording sessions",
    x = "Week after sowing", y = "Plant height (cm)"
  ) +
  theme_ng()

save_fig(fig2, "fig2_height_trajectories")

# --- Fig 3: Duncan groups for week-6 traits --------------------------------

fig3_data <- duncan |>
  filter(model == "adjusted") |>
  group_by(trait_label) |>
  mutate(accession = forcats::fct_reorder(accession, mean)) |>
  ungroup()

fig3 <- ggplot(fig3_data, aes(x = mean, y = accession)) +
  geom_col(fill = "#7C8DA6", width = 0.72) +
  geom_text(aes(label = duncan_group), hjust = -0.25, size = 2.5,
            colour = "grey25") +
  facet_wrap(~ trait_label, scales = "free", ncol = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Week-6 accession means with Duncan groups",
    subtitle = "Duncan's multiple range test, alpha = 0.05, adjusted for recording session. Accessions sharing a letter do not differ.",
    x = NULL, y = NULL
  ) +
  theme_ng() +
  theme(axis.text.y = element_text(size = 6.5))

save_fig(fig3, "fig3_duncan_groups", width = 11, height = 9)

# --- Fig 4: where the variance actually sits -------------------------------

fig4 <- vc |>
  select(trait_label, accession = repeatability, pot = pct_pot,
         residual = pct_residual) |>
  pivot_longer(-trait_label, names_to = "source", values_to = "share") |>
  mutate(
    source = factor(source, levels = c("accession", "pot", "residual"),
                    labels = c("Accession (genotype)", "Pot within accession",
                               "Plant / session (residual)")),
    trait_label = forcats::fct_reorder(
      trait_label, share * (source == "Accession (genotype)"), .fun = max)
  ) |>
  ggplot(aes(x = share, y = trait_label, fill = source)) +
  geom_col(width = 0.68) +
  scale_x_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = c("#B4462F", "#D9A441", "#C9CED6"), name = NULL) +
  labs(
    title = "How much of each trait is genotype?",
    subtitle = "Variance partition from value ~ reading + (1|accession) + (1|accession:pot)",
    x = "Share of total variance", y = NULL
  ) +
  theme_ng()

save_fig(fig4, "fig4_variance_components", width = 8.5, height = 5)

# --- Fig 5: flowering time against vegetative growth -----------------------

growth_means <- week6 |>
  group_by(accession, pvfta1) |>
  summarise(PH = mean(PH, na.rm = TRUE), NL = mean(NL, na.rm = TRUE),
            .groups = "drop")

fig5 <- pheno |>
  left_join(growth_means, by = c("accession", "pvfta1")) |>
  pivot_longer(c(PH, NL), names_to = "trait", values_to = "value") |>
  mutate(trait = recode(trait, PH = "Plant height (cm)",
                        NL = "Number of leaves")) |>
  ggplot(aes(x = dtff, y = value)) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey45",
              fill = "grey88", linewidth = 0.6, formula = y ~ x) +
  geom_point(aes(colour = pvfta1), size = 2.6) +
  facet_wrap(~ trait, scales = "free_y") +
  scale_colour_manual(values = pal, name = NULL) +
  labs(
    title = "Flowering time carries little information about vegetative size",
    subtitle = "Accession means at week 6 against days to first flowering",
    x = "Days to first flowering", y = NULL
  ) +
  theme_ng()

save_fig(fig5, "fig5_flowering_vs_growth", width = 9, height = 4.8)

message("Figures written to ", paths$figures)
