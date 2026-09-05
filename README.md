# PvFTa1 and flowering time in twenty Kano common bean accessions

A reproducible reanalysis, in R, of a screen-house trial of *Phaseolus
vulgaris* accessions collected from markets across ten Local Government Areas
of Kano State, Nigeria. The trial ran at the Department of Biological
Sciences, Bells University of Technology, Ota, from November 2024 to July
2025.

The wet-lab work, the trial and the original scoring are mine. This repository
rebuilds the statistics from the raw recording sheet, so that every number
below can be traced to a line of code and a row of data.

## The question

*PvFTa1* is the common bean homologue of *Arabidopsis* FLOWERING LOCUS T, and
is treated in the breeding literature as a candidate for early flowering. Five
of the twenty accessions were amplified for a 300 bp *PvFTa1* fragment. The
question this repository answers is narrow and answerable:

> Does carrying the *PvFTa1* amplicon track flowering time across this panel?

## The finding

**No. The locus is more conserved than the phenotypic spread predicts.**

Days to first flowering ranges over 18 days across the panel, from 28 days
(Kumbotso3) to 46 days (Dawanao12). The five accessions selected for
genotyping span **that entire 18-day range** — they include both the earliest
and the latest accession in the trial — and **all five gave the 300 bp
product**. Presence of the amplicon separates nothing.

| Accession | Days to first flowering | Days to 50% flowering | Class | PvFTa1 |
|---|---|---|---|---|
| Kumbotso3 | 28 | 38 | Early | 300 bp |
| Yankaba19 | 30 | 38 | Early | 300 bp |
| Danhassan2 | 36 | 42 | Intermediate | 300 bp |
| Gani8 | 36 | 38 | Early | 300 bp |
| Dawanao12 | 46 | 51 | Late | 300 bp |

A presence/absence assay at this locus therefore has no discriminating power
in this germplasm. Either the locus is fixed across these landraces, or the
variation that matters is sequence-level rather than presence-level. The
useful next step is to sequence the amplicon and look for polymorphism, not to
screen more accessions for the band.

Two further results support reading the marker result this way.

**Flowering time is the only sharply structured trait in the dataset.** Both
pots of every accession recorded exactly the same days to flowering — the
within-accession variance is zero for all sixteen accessions with a record.
Flowering time is crisp and repeatable here.

**Vegetative traits are not.** Fitting
`value ~ reading + (1|accession) + (1|accession:pot)` to the week-6
measurements puts accession identity at under 10% of total variance for every
trait:

| Trait | Accession | Pot within accession | Plant / session |
|---|---|---|---|
| Number of branches | 9.6% | 29.6% | 60.7% |
| Plant height | 8.0% | 0.0% | 92.0% |
| Leaf breadth | 7.5% | 0.0% | 92.5% |
| Number of leaves | 6.8% | 13.6% | 79.6% |
| Leaf length | 0.8% | 0.0% | 99.2% |
| Stem diameter | 0.0% | 0.0% | 100.0% |

The one-way ANOVA agrees: no week-6 growth trait shows a significant accession
effect (F between 0.12 and 1.77, all p > 0.05). This does **not** show that the
accessions are identical. It shows that with two pots per accession, heavy
missing data and two recording sessions that differ systematically, this trial
cannot resolve vegetative differences between them. Flowering time survives
that noise; plant size does not.

So the panel's genuinely informative phenotype is flowering time, flowering
time varies over 18 days, and the marker is present in every accession tested
across that whole range.

## Figures

| | |
|---|---|
| `fig1_flowering_range.png` | Flowering window per accession, PvFTa1 status marked. The headline figure. |
| `fig2_height_trajectories.png` | Plant height, weeks 2-6 |
| `fig3_duncan_groups.png` | Week-6 accession means with Duncan groups |
| `fig4_variance_components.png` | Where the variance sits for each trait |
| `fig5_flowering_vs_growth.png` | Flowering time against vegetative size |

## Methods

- **Design.** Completely randomised: 20 accessions x 2 pots, thinned to 2
  plants per pot at week 2. 80 plants.
- **ANOVA and Duncan's multiple range test** (`agricolae`) on week-6 traits,
  alpha = 0.05, reproducing the original analysis. Fitted twice: once as
  reported (`value ~ accession`) and once adjusting for recording session
  (`value ~ accession + reading`).
- **Mixed models** (`lme4`) to partition variance across accession, pot within
  accession, and the plant/session residual, and to estimate the repeatability
  of a single plant measurement.
- **Flowering** is treated as one fixed value per accession, because the
  within-accession variance is exactly zero and an F ratio would have a zero
  denominator.

## Reproducing

```
Rscript R/00_packages.R    # once, installs dependencies
Rscript run_all.R          # rebuilds every table and figure
```

Written against R 4.6.1. `lme4` needs `cmake` on the path to build `nloptr`
from source (`brew install cmake` on macOS).

```
data/raw/         unedited export of the recording sheet
data/derived/     tidied tables, rebuilt by 01_tidy_data.R
R/                numbered scripts, run in order
outputs/tables/   CSV results
outputs/figures/  PNG figures
```

`data/README.md` documents the sheet layout and every cleaning decision.

## Departures from the original write-up

The thesis analysed these data by hand. Rebuilding from the sheet surfaced
several things worth stating plainly.

**Accession effects on growth traits do not reach significance when
recomputed.** The thesis describes significant variability in plant height,
leaf number and stem diameter at week 6. Recomputed from the sheet, no growth
trait clears p = 0.05. The Duncan letters in `outputs/tables/duncan_week6.csv`
are reported for completeness, but a mean separation test after a
non-significant F should not be read as evidence of real differences.

**Four accessions have no flowering record in the sheet.** Garu4, Round13,
Wudil15 and Lahadin17 have empty flowering columns; the thesis reports values
for all twenty. This analysis uses the sixteen the sheet supports.

**The two recording sessions are not comparable.** Block 2 runs about three
times higher than block 1 for plant height, and stem diameter differs by a
factor of about seven, which is a change of scale rather than of plants.
Stem-diameter results should not be read across sessions at all.

**The genotyped set is described inconsistently in the thesis.** Section 3.4
lists accessions 3a, 19a, 2a, 12a and 10a; section 4.5 and the Plate 2 gel
lanes list Kumbotso3, Danhassan2, Gani8, Dawanao12 and Yankaba19 — Gani8 in
place of Ruwan10. This repository follows the gel, which is the primary
record, and the discussion chapter agrees with it. Note that the substitution
does not change the finding: both candidate sets span the panel's full
flowering range.

## Open question

Whether the second block of the sheet is a second plant per pot or a later
recording round is not resolvable from the file alone. The models treat
`reading` as a fixed effect either way, so the variance partition and the
ANOVA hold under both readings; only the label on the residual term changes
("plant within pot" versus "measurement occasion"). The flowering result does
not touch the second block at all.

## Licence

Code MIT. Data available for reuse with attribution.
