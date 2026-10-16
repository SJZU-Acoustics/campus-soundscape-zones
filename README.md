R code for reproducing the statistical analyses, figures and tables for the manuscript "Soundscape appropriateness, comfort and richness across the functional zones of a campus landscape".

## Requirements

- R 4.5+ (developed and verified on R 4.5.3)

- CRAN packages: `tidyverse`, `readxl`, `lme4`, `lmerTest`, `pbkrtest`, `emmeans`, `performance`, `broom`, `broom.mixed`, `ordinal`, `patchwork`, `ggtext`, `cowplot`, `ragg`

- Install:

  ```r
  install.packages(c("tidyverse", "readxl", "lme4", "lmerTest", "pbkrtest", "emmeans",
                     "performance", "broom", "broom.mixed", "ordinal", "patchwork",
                     "ggtext", "cowplot", "ragg"))
  ```

- Python 3 (standard library only), once, to fetch the optional OpenStreetMap layer of Figure 1.

- The figures are set in Helvetica. Where it is not installed, the PNGs use a substitute font and differ from the published figures in typeface only.

- No non-standard hardware is required. A complete run took 6.5 minutes on an Apple M3 MacBook Air; the permutation test of A09b takes about a quarter of it, and the cluster bootstraps of A13, A14b, A18 and A21 and the index screens of A17 most of the rest.

## Data

The analysis reads a single input: the Zenodo workbook
`Campus_soundwalk_summer_2026_perception_acoustics_visual.xlsx` (CC BY 4.0).

1. Download it from Zenodo: https://doi.org/10.5281/zenodo.23420131.
2. Place the `.xlsx` file in the `data/` folder (see `data/README.md`).

The workbook holds the 990 ratings that 30 students gave at 208 stops on eight sites of one campus: source salience, the eight ICSAQ descriptors with Comfort and Richness, appropriateness, noisiness, and the stops' sound levels, psychoacoustic metrics and panorama composition. It also holds the evaluator, session and site tables, and for each of the 143 calibrated recordings the ecoacoustic indices and the detection scores of a pretrained audio tagger.

The recordings and panoramas themselves are not part of the deposit. Every variable used by the models is deposited, so nothing here needs the media.

Figure 1a also draws an OpenStreetMap base layer, which is fetched separately (see `data/osm/README.md`).

## File structure

- `run_all.R` — master script: runs every analysis module, then builds every display item of the paper.
- `code/load_data.R` — single data entry: maps each analysis table to its workbook sheet and reads it through a CSV cache, so column types match the working pipeline.
- `code/helpers.R` — shared data entry, variable sets (outcomes, sources and source groups, acoustic and visual variables) and model helpers used by every module.
- `code/a*.R` — the analysis modules:
  - A01 inventory of the sample; A02 the two acoustic windows (XL2 measurement and 30-s region);
  - A03 zone profiles; A04 acoustic indicators, A05 perceived sources and A06 visual composition as predictors of the four outcomes, with A06b re-run without session W04;
  - A07 day type and time slot; A08 person-level variables; A09 variance partitioning and A09b the evaluator random-slope permutation test;
  - A10b the zone × Comfort test on appropriateness; A11 site-held-out prediction by feature family; A12 robustness of the five-feature model, including cumulative-link models;
  - A13 the natural–mechanical exchange ratio; A14b zone-specific level slopes and neutral levels; A15 Richness and source diversity; A16 and A16b the ICSAQ quadrants;
  - A17 ecoacoustic indices; A18 detected sources; A19 the structure of the source groups; A20 sensitivity checks of the headline results; A21 the other listeners' source reports.
- `code/style.R` — plot style and the reader for the module result tables.
- `code/build_figures.R` — Figures 1–6.
- `code/build_tables.R` — Tables 1–2 and Supplementary Tables S1–S34 (LaTeX fragments).
- `code/fetch_osm.py` — fetches the OpenStreetMap layer of Figure 1a.

## Usage

From the repository root:

```bash
python3 code/fetch_osm.py   # once, for Figure 1
Rscript run_all.R
```

Outputs are written to:

- `output/figures/` — Figures 1–6 (PNG, 600 dpi)
- `output/tables/` — Tables 1–2 (LaTeX fragments)
- `output/tables/si/` — Supplementary Tables S1–S34 (LaTeX fragments; some tables combine two fragments)
- `output/<analysis>/` — the full result tables of each analysis module

To keep the console log alongside the outputs:

```bash
Rscript run_all.R 2>&1 | tee output/run_log.txt
```

`output/` is produced at run time and is safe to delete. It also holds `_workbook_cache/`, the extracted workbook sheets; delete it to force a fresh read.

## Verification

Run against the deposited workbook, this pipeline reproduces the paper's display items: Tables 1–2 and all 34 Supplementary Tables (the 49 LaTeX fragments) and Figures 1, 2, 3 and 6 are **byte-identical** to the versions in the paper (Figure 1 with the OpenStreetMap layer used for the paper). Figures 4 and 5 match apart from anti-aliasing: 694 and 37 of about 15 million pixels differ, by at most 2 of 255 in one colour channel. The bootstraps and the permutation test draw from fixed seeds set in the scripts and run serially, so they regenerate exactly.

Workbook cells hold 16 significant digits, so some computed values in the deposit differ from the working data in their last digit, and the unrounded columns of the module result tables differ from the working pipeline in their last digits. Across the 26,614 numeric result cells, the two largest relative differences, about 3e-4, are on the difference of two nearly equal R² values (1.0604e-4 against 1.0600e-4) and on a p-value of 1.5e-5. All other cells differ by less than 3e-5, and estimates and test statistics by less than 2e-5: an input perturbation of about 1e-16 grows through model fitting and the numerical derivatives behind the degrees of freedom. Every quantity reported in the paper is printed to far fewer digits, and none of them moves.

## Notes

- Each module runs in its own R session, reads the deposited tables through `code/helpers.R` and writes only into its analysis folder under `output/`. A06b reads the A06 screen and therefore runs after it; `build_tables.R` reads the Figure 3a estimates that `build_figures.R` writes. The other modules can be run individually and in any order.
- The rating models are linear mixed models with crossed random intercepts for evaluator and stop, fitted by REML with Kenward–Roger degrees of freedom for marginal means and trends (Satterthwaite in A03); model comparisons refit by ML.
- The primary acoustic layer is the 30-s region of the 666 ratings at stops with a calibrated recording; the 990-rating models use the XL2 level.
- Benjamini–Hochberg false-discovery-rate control is applied within the declared families of each module; single pre-specified tests report their own p.
- Site-held-out prediction (A11, A17–A21) refits every model on seven sites and predicts the eighth; the stop-mean R² compares observed and predicted stop means.

## License

Code in this repository is released under the MIT License (see `LICENSE`). The input data are archived separately under CC BY 4.0 at Zenodo (https://doi.org/10.5281/zenodo.23420131). The OpenStreetMap layer of Figure 1a is © OpenStreetMap contributors under the Open Database License.
