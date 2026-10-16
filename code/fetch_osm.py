"""Fetch the campus base-map features from the OpenStreetMap Overpass API and save
them as one CSV of polylines (data/osm/osm_campus_features.csv) for Figure 1a.
Data (c) OpenStreetMap contributors, ODbL. Run once from the repository root:
    python3 code/fetch_osm.py"""
import csv, json, pathlib, sys, urllib.request, urllib.parse
OUT = pathlib.Path(__file__).resolve().parents[1] / "data" / "osm" / "osm_campus_features.csv"
S, W, N, E = 41.7355, 123.5020, 41.7475, 123.5180   # campus bbox (WGS84)
Q = f"""[out:json][timeout:90];
(
  way["highway"]({S},{W},{N},{E});
  way["building"]({S},{W},{N},{E});
  way["natural"="water"]({S},{W},{N},{E});
  way["waterway"]({S},{W},{N},{E});
  way["leisure"]({S},{W},{N},{E});
  way["landuse"]({S},{W},{N},{E});
  way["amenity"="university"]({S},{W},{N},{E});
  relation["natural"="water"]({S},{W},{N},{E});
  relation["amenity"="university"]({S},{W},{N},{E});
);
out geom;"""
import time
MIRRORS = ["https://overpass.kumi.systems/api/interpreter", "https://overpass-api.de/api/interpreter", "https://lz4.overpass-api.de/api/interpreter"]
j = None
for url in MIRRORS:
    try:
        req = urllib.request.Request(url,
                                     data=urllib.parse.urlencode({"data": Q}).encode(),
                                     headers={"User-Agent": "p35-campus-map/1.0"})
        j = json.loads(urllib.request.urlopen(req, timeout=150).read()); print("fetched from", url); break
    except Exception as e:
        print("mirror failed:", url, e); time.sleep(3)
if j is None: sys.exit("all Overpass mirrors failed")
rows = []
def cls(t):
    if "building" in t: return "building"
    if t.get("natural") == "water" or "waterway" in t: return "water"
    if t.get("highway") in ("motorway", "trunk", "primary", "secondary", "tertiary"): return "road_major"
    if "highway" in t: return "road_minor"
    if t.get("leisure") in ("pitch", "sports_centre", "stadium", "track"): return "sport"
    if t.get("leisure") in ("park", "garden") or t.get("landuse") in ("grass", "farmland", "meadow", "forest", "recreation_ground"): return "green"
    if t.get("amenity") == "university": return "campus"
    return None
for el in j["elements"]:
    t = el.get("tags", {}); c = cls(t)
    if c is None: continue
    if el["type"] == "way":
        for i, g in enumerate(el.get("geometry", [])):
            rows.append((el["id"], 0, i, c, t.get("highway", ""), t.get("name", ""), g["lon"], g["lat"]))
    elif el["type"] == "relation":
        for m in el.get("members", []):
            if m["type"] != "way" or m.get("role") not in ("outer", ""): continue
            for i, g in enumerate(m.get("geometry", [])):
                rows.append((el["id"], m["ref"], i, c, "", t.get("name", ""), g["lon"], g["lat"]))
OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["osm_id", "member_ref", "seq", "class", "highway", "name", "lon", "lat"]); w.writerows(rows)
from collections import Counter
print("elements", len(j["elements"]), "rows", len(rows)); print(Counter(r[3] for r in rows))
