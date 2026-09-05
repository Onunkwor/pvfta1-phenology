# 01_tidy_data.R ------------------------------------------------------------
# Turn the screen-house recording sheet into tidy, analysis-ready tables.
#
# The exported sheet holds two stacked blocks separated by a "Second reading"
# marker row. Each block carries the same 40 rows: 20 accessions x 2 pots
# (suffix a/b). Two plants were retained per pot after thinning at week 2,
# which is where the 80 plants reported in the thesis come from.

source("R/00_packages.R")

raw_path <- file.path(paths$raw, "morphological_data_export.csv")
raw_lines <- read_lines(raw_path)

# Locate the block boundaries rather than hard-coding row numbers, so the
# pipeline survives a re-export of the sheet.
marker <- grep("Second reading", raw_lines)
stopifnot(length(marker) == 1L)

header_rows <- grep("^Cultivar_Replicate", raw_lines)
stopifnot(length(header_rows) == 2L)

read_block <- function(start, end) {
  read_csv(
    I(raw_lines[c(start, seq(start + 1L, end))]),
    show_col_types = FALSE,
    name_repair    = "unique_quiet"
  )
}

block1 <- read_block(header_rows[1], marker - 2L)              # drop trailing blank
block2 <- read_block(header_rows[2], length(raw_lines))

# ---------------------------------------------------------------------------
# Reference tables
# ---------------------------------------------------------------------------

# Display names follow Table 2 of the thesis so figures can be cross-read
# against it.
accession_names <- tibble(
  accession_code = 1:20,
  accession = c(
    "Kura1", "Danhassan2", "Kumbotso3", "Garu4", "Danzabura5",
    "Bichi6", "Sumailia7", "Gani8", "Rano9", "Ruwan10",
    "Dawakin11", "Dawanao12", "Round13", "Kwanar14", "Wudil15",
    "Darki16", "Lahadin17", "Dawaki18", "Yankaba19", "Sabuwar20"
  )
)

# The sheet spells several LGAs inconsistently ("Samalia" for Sumaila,
# "Dawakin toea" for Dawakin Tofa, "Dawakin Kadu" for Dawakin Kudu).
# Normalise against Table 1 of the thesis.
lga_fixes <- c(
  "Samalia Local Gov"         = "Sumaila",
  "Kura Local Gov Kano"       = "Kura",
  "Kumbotso Local Gov Kano"   = "Kumbotso",
  "Bichi Local Gov Kano"      = "Bichi",
  "Rano Local Gov"            = "Rano",
  "Dawakin toea Local Gov"    = "Dawakin Tofa",
  "Kiru Local Gov"            = "Kiru",
  "Wudil Local Gov"           = "Wudil",
  "Dawakin Kadu"              = "Dawakin Kudu",
  "Nasarawa Local Kano"       = "Nasarawa"
)

# Accessions whose PvFTa1 amplicon was scored on the gel (Plate 2, section 4.5).
pvfta1_amplified <- c(2L, 3L, 8L, 12L, 19L)

# ---------------------------------------------------------------------------
# Growth traits -> long format
# ---------------------------------------------------------------------------

trait_labels <- c(
  PH = "Plant height (cm)",
  NL = "Number of leaves",
  NB = "Number of branches",
  LL = "Leaf length (cm)",
  LB = "Leaf breadth (cm)",
  SD = "Stem diameter"
)

tidy_growth <- function(df, reading) {
  df |>
    rename(pot_id = Cultivar_Replicate) |>
    mutate(lga_raw = .data[[grep("^Local Gov", names(df), value = TRUE)[1]]]) |>
    select(pot_id, lga_raw, matches("^Wk[2-6](PH|Ph|NL|NB|LL|LB|SD)$")) |>
    filter(!is.na(pot_id), str_detect(pot_id, "^[0-9]+[ab]$")) |>
    pivot_longer(
      cols      = -c(pot_id, lga_raw),
      names_to  = c("week", "trait"),
      names_pattern = "^Wk([2-6])(PH|Ph|NL|NB|LL|LB|SD)$",
      values_to = "value"
    ) |>
    mutate(
      reading = reading,
      week    = as.integer(week),
      trait   = toupper(trait),
      value   = suppressWarnings(as.numeric(value))
    )
}

growth_long <- bind_rows(
  tidy_growth(block1, 1L),
  tidy_growth(block2, 2L)
) |>
  mutate(
    accession_code = as.integer(str_extract(pot_id, "^[0-9]+")),
    pot            = str_extract(pot_id, "[ab]$"),
    lga            = unname(lga_fixes[str_trim(lga_raw)])
  ) |>
  left_join(accession_names, by = "accession_code") |>
  mutate(
    # A plant that was never there was recorded as a row of zeros. Zero is a
    # real value for branch count but impossible for a height, a leaf count,
    # a leaf dimension or a stem diameter, so those become missing.
    value = if_else(trait != "NB" & value == 0, NA_real_, value),
    trait_label = unname(trait_labels[trait]),
    pvfta1 = if_else(accession_code %in% pvfta1_amplified,
                     "Amplified", "Not genotyped"),
    plant_id = paste0(pot_id, "-r", reading)
  ) |>
  select(accession_code, accession, lga, pot_id, pot, reading, plant_id,
         week, trait, trait_label, value, pvfta1) |>
  arrange(accession_code, pot, reading, week, trait)

stopifnot(
  n_distinct(growth_long$accession_code) == 20L,
  n_distinct(growth_long$plant_id) == 80L
)

# ---------------------------------------------------------------------------
# Flowering traits (recorded once per pot, in the first block only)
# ---------------------------------------------------------------------------

flowering <- block1 |>
  rename(pot_id = Cultivar_Replicate) |>
  select(
    pot_id,
    dtff = Days_To_First_Flowering,
    dt50 = `Days_To_50%_Flowering`,
    category_raw = starts_with("Flowring Category")
  ) |>
  filter(str_detect(pot_id, "^[0-9]+[ab]$")) |>
  mutate(
    accession_code = as.integer(str_extract(pot_id, "^[0-9]+")),
    pot   = str_extract(pot_id, "[ab]$"),
    dtff  = suppressWarnings(as.numeric(dtff)),
    dt50  = suppressWarnings(as.numeric(dt50)),
    # The sheet's own class bands, restated from the column header:
    # early 38-40, intermediate 41-45, late 46-51 days to 50% flowering.
    category = str_to_title(str_trim(category_raw))
  ) |>
  left_join(accession_names, by = "accession_code") |>
  mutate(pvfta1 = if_else(accession_code %in% pvfta1_amplified,
                          "Amplified", "Not genotyped")) |>
  select(accession_code, accession, pot_id, pot, dtff, dt50, category, pvfta1) |>
  arrange(accession_code, pot)

# ---------------------------------------------------------------------------
# Week-6 wide table, the frame the classical analysis runs on
# ---------------------------------------------------------------------------

week6 <- growth_long |>
  filter(week == 6L) |>
  select(accession_code, accession, lga, pot_id, pot, reading, plant_id,
         trait, value, pvfta1) |>
  pivot_wider(names_from = trait, values_from = value)

write_csv(growth_long, file.path(paths$derived, "growth_long.csv"))
write_csv(flowering,   file.path(paths$derived, "flowering.csv"))
write_csv(week6,       file.path(paths$derived, "week6_wide.csv"))

message(sprintf(
  "Tidied %d growth observations across %d plants; %d flowering records.",
  nrow(growth_long), n_distinct(growth_long$plant_id), nrow(flowering)
))
