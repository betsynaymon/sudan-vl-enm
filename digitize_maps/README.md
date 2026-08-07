# Map Digitization

Extracts georeferenced occurrence points from published map figures for inclusion in the VL occurrence dataset. Each map is configured in `map_config.json` with two geographic reference points; the script uses an interactive click-based workflow to collect point locations from the image and transform pixel coordinates to latitude/longitude.

## Source figures not included

The `maps/` directory holds PNG extracts of figures from published journal articles (Mueller et al. 2012, Gerstl et al. 2006, van den Bogaart et al. 2013, Elnaiem et al. 2003, Hassan et al. 2020). These are copyrighted and are not redistributed in this repository, so the workflow below cannot be re-run directly from a clone.

What this directory provides is a record of how the digitized coordinates were produced: the script (`digitize_map.py`) and the per-map configuration (`map_config.json`) are included. The coordinates themselves are already incorporated into `data/raw/compiled_vl_presences.csv`, the file the analysis pipeline consumes — the digitization does not need to be re-run. To reproduce it, obtain the figures from the cited papers and place them in `maps/` following the paths in `map_config.json`.

## Folder structure

```
digitize_maps/
├── digitize_map.py       # Interactive digitization script
├── map_config.json       # Map entries: reference points, citations, years
├── maps/                 # Source figures — not committed
└── output/               # Generated CSVs — not committed; feed compiled_vl_presences.csv
```

`maps/` and `output/` are gitignored; only the script and config are tracked.

## Requirements

- Python 3
- numpy, matplotlib, Pillow

## How it works

The script converts pixel coordinates to geographic coordinates using a two-point scale-and-translate transform — independent x and y scaling plus translation, with no rotation or warping. For each map, you click two known locations (defined in the config), then click all the data points you want to extract. The transform assumes the map projection is approximately linear over the area shown, which is reasonable at the sub-national scale of these figures.

## Usage

**List configured maps and their digitization status:**

```bash
python digitize_map.py --list
```

This prints all entries from `map_config.json` with a ✓ next to any that already have output CSVs.

**Digitize a single map:**

```bash
python digitize_map.py mueller_et_al_2012_fig2
```

This opens the map image in a matplotlib window and walks through two steps:

1. **Click the two reference points** in the order shown in the terminal. These are the known geographic anchors (e.g., a named town) that calibrate the pixel-to-coordinate transform. The script prints a sanity check confirming the reference points map back correctly.

2. **Click all data points** (village locations, case markers, etc.). Left-click to add points; right-click or press Enter when finished.

The output CSV is saved to `output/digitized_<map_key>.csv`.

**Merge all individual CSVs into one file:**

```bash
python digitize_map.py --merge
```

This combines all per-map CSVs in `output/` into `output/digitized_all_maps.csv` with renumbered `point_id` values.

## Output columns

| Column | Description |
|--------|-------------|
| `point_id` | Sequential ID (per map, or renumbered in merged file) |
| `pixel_x` | Clicked x-coordinate in the source image |
| `pixel_y` | Clicked y-coordinate in the source image |
| `longitude` | Transformed longitude (decimal degrees) |
| `latitude` | Transformed latitude (decimal degrees) |
| `source` | Publication citation (from config) |
| `year` | Year of data collection (from config) |
| `map_key` | Config key identifying the source map |

## Adding a new map

1. Save the map image to `maps/`
2. Add an entry to `map_config.json` with the image filename, two reference points (identifiable locations with known coordinates), the source citation, and the data year
3. Run `python digitize_map.py <new_key>`

## Notes

- The digitized points in `output/` are the authoritative coordinates used in the compiled occurrence dataset (`data/raw/compiled_vl_presences.csv`). This directory documents provenance — the config and code that produced those rows — but is not part of the main analysis pipeline and does not need to be re-run.
- Pixel coordinates are retained in the output to allow verification against the source image.
- The transform is a simple scale-and-translate, which is sufficient for the map scales involved but would not be appropriate for large-area or heavily projected maps.