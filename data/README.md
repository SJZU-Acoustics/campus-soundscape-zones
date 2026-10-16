# Data folder

## The Zenodo workbook (required)

Download `Campus_soundwalk_summer_2026_perception_acoustics_visual.xlsx` from Zenodo (CC BY 4.0; the DOI is given in the
paper's Data availability statement) and place it in **this folder**:

```
data/Campus_soundwalk_summer_2026_perception_acoustics_visual.xlsx
```

The workbook is not tracked in this repository; Zenodo is its citable home.

Sheets used by the analysis:

| Sheet | What it holds |
|---|---|
| `records` | one row per rating (990): the evaluator, session, site, zone and stop; the salience of 24 sound sources; the eight ICSAQ descriptors with Comfort and Richness; the ratings of appropriateness, noisiness and the visual and sound environments; and the stop's sound levels, psychoacoustic metrics and panorama composition |
| `people` | the 30 evaluators: sex, age band, education, walking group, rater type, WHO-5 and noise-sensitivity items and scores |
| `environments` | the 208 stops: session, site, visit time, XL2 levels, the B&K metrics of the 143 calibrated recordings, and panorama composition |
| `sessions` | the 14 walking sessions |
| `sites` | the eight survey sites with their functional zones and WGS84 coordinates |
| `ecoacoustic_indices` | 68 ecoacoustic indices and levels (scikit-maad) for each of the 143 recordings, on the 30-s region and on the whole recording |
| `detected_sources` | PANNs CNN14 detection scores for the 24 sources and three source groups, for the same two windows |

The workbook also holds `audioset_source_map` (the assignment of the 527 AudioSet classes to the 24 sources, which
defines the detection scores), a sheet summary, a variable dictionary and a value-code sheet.

The recordings, the raw sound-level-meter files and the panoramas are not part of the deposit. Every variable that the
analyses use is deposited, so nothing here needs the media.

## `osm/` (optional)

The OpenStreetMap base layer of Figure 1a; see `osm/README.md`.
