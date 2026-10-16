# OpenStreetMap base layer (optional; Figure 1 only)

Figure 1a draws the campus buildings, roads, water and green space from OpenStreetMap.
The layer is not part of this repository or of the Zenodo deposit: it is third-party data,
© OpenStreetMap contributors, available under the Open Database License (ODbL).

To draw Figure 1, fetch the layer once from the repository root (Python 3, standard library only):

```bash
python3 code/fetch_osm.py
```

This queries the Overpass API for the campus bounding box and writes

```
data/osm/osm_campus_features.csv
```

Without the file, `run_all.R` skips Figure 1 with a message and builds everything else.
OpenStreetMap is edited continually, so a fresh fetch can differ in detail from the base map
in the published figure; the survey points and panel b come from the deposited data.
