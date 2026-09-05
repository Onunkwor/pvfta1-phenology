# Data provenance

## Source

`raw/morphological_data_export.csv` is an unedited CSV export of the Google
Sheet `morphological_data`, the recording sheet kept during the screen-house
trial at the Department of Biological Sciences, Bells University of
Technology, Ota, Ogun State, between November 2024 and July 2025.

Nothing in `raw/` is modified. Every correction, unit fix and exclusion is
applied in `R/01_tidy_data.R` and lands in `derived/`, so the cleaning is
reviewable rather than baked into the inputs.

## Field trial

Twenty *Phaseolus vulgaris* accessions were bought from markets across ten
Local Government Areas of Kano State, Nigeria. Forty pots (20 accessions x 2
pots, suffixed `a` and `b`) were filled with ~6 kg sandy loam and arranged in
rows of four at 50 cm between rows and 30 cm within rows. Four seeds were sown
per pot and thinned to two plants at two weeks. Completely randomised design,
80 plants in total.

Weekly measurements ran from week 2 to week 6: plant height, number of leaves,
number of branches, leaf length, leaf breadth, stem diameter. Days to first
flowering and days to 50% flowering were scored once per pot.

## Sheet layout

The export holds two stacked blocks of 40 rows each, separated by a row
reading `Second reading`:

| Block | Rows | Contents |
|-------|------|----------|
| 1 | 1-40 | Weeks 2-6 growth traits, flowering dates, pod and seed traits |
| 2 | 44-83 | Weeks 2-6 growth traits only |

`R/01_tidy_data.R` finds these boundaries by searching for the marker and the
two header rows rather than hard-coding line numbers.

## Known issues in the sheet

These are recorded here rather than silently repaired, because each one
changes how a result should be read.

**1. The two blocks are not interchangeable.** Pooled across all accessions,
block 2 runs about three times higher than block 1 for plant height and about
1.2 to 1.8 times higher for leaf traits, while branch counts are effectively
identical. Two plants thinned from the same pot cannot differ that
systematically. Block 2 is therefore treated as a separate recording session,
and `reading` enters every model as a fixed effect. See the open question in
the top-level README.

**2. Stem diameter changes scale between blocks.** Block 1 values sit around
0.5 and block 2 around 3.7, a ratio close to 7. The two blocks are not on a
common scale, and stem-diameter results should not be read across them.

**3. Absent plants were recorded as rows of zeros.** Eight rows of block 2 are
all-zero. Zero is a legitimate branch count but not a legitimate height, leaf
count, leaf dimension or stem diameter, so a zero in any trait other than
`NB` is converted to missing.

**4. Four accessions have no flowering record.** Garu4, Round13, Wudil15 and
Lahadin17 have empty flowering columns. Table 2 of the thesis reports days to
flowering for all twenty, so the values for these four come from a source
outside this sheet. This pipeline analyses the sixteen the sheet supports.

**5. Both pots of an accession always agree exactly on flowering.** Every
accession with a flowering record has identical values in pots `a` and `b`.
Within-accession variance is exactly zero, which is why days to flowering is
handled as a single fixed value per accession rather than through ANOVA.

**6. Place names are spelt inconsistently.** `Samalia` for Sumaila,
`Dawakin toea` for Dawakin Tofa, `Dawakin Kadu` for Dawakin Kudu. Normalised
against Table 1 of the thesis in `01_tidy_data.R`.

## Derived files

| File | One row per |
|------|-------------|
| `derived/growth_long.csv` | accession x pot x reading x week x trait |
| `derived/week6_wide.csv` | accession x pot x reading, week-6 traits in columns |
| `derived/flowering.csv` | accession x pot |
| `derived/phenology_by_accession.csv` | accession |
| `derived/phenology_growth.csv` | accession, phenology joined to week-6 means |
