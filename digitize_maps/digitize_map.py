 """
digitize_map.py — Config-driven map digitisation for VL occurrence points
=========================================================================
Usage:
    python digitize_map.py <map_key>           Digitise a single map
    python digitize_map.py --list              Show all available maps
    python digitize_map.py --merge             Merge all individual CSVs into one

Config lives in digitize_maps/map_config.json
Each entry needs: image filename, two reference points with lat/lon, source, year.

Output: one CSV per map in output/ (digitized_<map_key>.csv) with columns:
    point_id, pixel_x, pixel_y, longitude, latitude, source, year
"""

import sys
import os
import json
import csv
import glob
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

# ------------ Paths -------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CONFIG_FILE = os.path.join(SCRIPT_DIR, "map_config.json")
MAPS_DIR = os.path.join(SCRIPT_DIR, "maps")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "output")


def load_config():
    """Load map_config.json from the same directory as this script."""
    if not os.path.exists(CONFIG_FILE):
        print(f"Error: Config file not found at {CONFIG_FILE}")
        print("Create map_config.json with your map entries.")
        sys.exit(1)
    with open(CONFIG_FILE, "r") as f:
        return json.load(f)


def compute_affine_transform(px_refs, geo_refs):
    """
    Build scale-and-translate transform from two reference points.
    px_refs:  [[px_x1, px_y1], [px_x2, px_y2]]  (pixel coordinates)
    geo_refs: [[lon1, lat1], [lon2, lat2]]        (geographic coordinates)
    Returns a function: (px_x, px_y) -> (lon, lat)
    """
    px1, px2 = np.array(px_refs[0]), np.array(px_refs[1])
    geo1, geo2 = np.array(geo_refs[0]), np.array(geo_refs[1])

    # Scale: degrees per pixel in x and y
    dpx = px2 - px1
    dgeo = geo2 - geo1

    scale_x = dgeo[0] / dpx[0] if dpx[0] != 0 else 0
    scale_y = dgeo[1] / dpx[1] if dpx[1] != 0 else 0

    # Translate: offset so ref1 maps correctly
    offset_x = geo1[0] - scale_x * px1[0]
    offset_y = geo1[1] - scale_y * px1[1]

    def transform(px_x, px_y):
        lon = scale_x * px_x + offset_x
        lat = scale_y * px_y + offset_y
        return lon, lat

    return transform


def digitize_one_map(map_key, config_entry):
    """Run the interactive digitisation workflow for a single map."""

    image_file = config_entry["image"]
    ref_points = config_entry["ref_points"]
    source = config_entry["source"]
    year = config_entry.get("year", "")
    notes = config_entry.get("notes", "")

    # Resolve image path (check maps/ subdir, then script dir, then cwd)
    img_path = os.path.join(MAPS_DIR, image_file)
    if not os.path.exists(img_path):
        img_path = os.path.join(SCRIPT_DIR, image_file)
    if not os.path.exists(img_path):
        img_path = os.path.join(os.getcwd(), image_file)
    if not os.path.exists(img_path):
        print(f"Error: Cannot find image '{image_file}'")
        print(f"  Looked in: {MAPS_DIR}")
        print(f"         and: {SCRIPT_DIR}")
        print(f"         and: {os.getcwd()}")
        return

    img = Image.open(img_path)

    ref1 = ref_points[0]
    ref2 = ref_points[1]

    # ── Step 1: Click reference points ──
    print("=" * 60)
    print(f"MAP: {source}")
    if notes:
        print(f"Notes: {notes}")
    print()
    print("STEP 1: Click the TWO reference points IN ORDER")
    print(f"  Click 1: {ref1['name']} ({ref1['lat']}N, {ref1['lon']}E)")
    print(f"  Click 2: {ref2['name']} ({ref2['lat']}N, {ref2['lon']}E)")
    print("=" * 60)

    fig, ax = plt.subplots(1, 1, figsize=(12, 16))
    ax.imshow(img)
    ax.set_title(
        f"{source}\n"
        f"CLICK 2 REFERENCE POINTS:\n"
        f"1) {ref1['name']}  2) {ref2['name']}",
        fontsize=13,
    )

    ref_pixels = plt.ginput(n=2, timeout=0)
    plt.close()

    if len(ref_pixels) != 2:
        print("Error: Need exactly 2 reference points")
        return

    print(f"  Ref 1 pixel: ({ref_pixels[0][0]:.1f}, {ref_pixels[0][1]:.1f})")
    print(f"  Ref 2 pixel: ({ref_pixels[1][0]:.1f}, {ref_pixels[1][1]:.1f})")

    # Build transform
    px_refs = [list(ref_pixels[0]), list(ref_pixels[1])]
    geo_refs = [[ref1["lon"], ref1["lat"]], [ref2["lon"], ref2["lat"]]]
    transform = compute_affine_transform(px_refs, geo_refs)

    # Sanity check
    for i, (px, geo) in enumerate(zip(px_refs, geo_refs)):
        computed = transform(px[0], px[1])
        print(
            f"  Ref {i+1} check: expected ({geo[0]:.4f}, {geo[1]:.4f}), "
            f"got ({computed[0]:.4f}, {computed[1]:.4f})"
        )

    # ── Step 2: Click data points ──
    print()
    print("=" * 60)
    print("STEP 2: Click ALL village/case location circles")
    print("  Left-click to add points")
    print("  Right-click or press Enter when DONE")
    print("=" * 60)

    fig, ax = plt.subplots(1, 1, figsize=(12, 16))
    ax.imshow(img)

    # Show reference points as red crosses
    for px in px_refs:
        ax.plot(px[0], px[1], "r+", markersize=15, markeredgewidth=2)

    ax.set_title(
        f"{source}\n"
        "CLICK ALL POINTS — right-click or Enter when done",
        fontsize=13,
    )

    data_pixels = plt.ginput(n=-1, timeout=0)
    plt.close()

    print(f"\n  Collected {len(data_pixels)} points")

    if len(data_pixels) == 0:
        print("  No points collected — skipping.")
        return

    # ── Step 3: Transform and save ──
    results = []
    for i, (px_x, px_y) in enumerate(data_pixels):
        lon, lat = transform(px_x, px_y)
        results.append(
            {
                "point_id": i + 1,
                "pixel_x": round(px_x, 1),
                "pixel_y": round(px_y, 1),
                "longitude": round(lon, 5),
                "latitude": round(lat, 5),
                "source": source,
                "year": year,
                "map_key": map_key,
            }
        )
        print(
            f"  Point {i+1}: pixel({px_x:.0f}, {px_y:.0f}) -> "
            f"({lat:.4f}N, {lon:.4f}E)"
        )

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    output_csv = os.path.join(OUTPUT_DIR, f"digitized_{map_key}.csv")
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    print(f"\n  Saved {len(results)} points to {output_csv}")
    print("  Done!\n")


def merge_all():
    """Merge all digitized_*.csv files into one combined file."""
    pattern = os.path.join(OUTPUT_DIR, "digitized_*.csv")
    csv_files = sorted(glob.glob(pattern))

    # Exclude the merged file itself if it exists
    merged_name = os.path.join(OUTPUT_DIR, "digitized_all_maps.csv")
    csv_files = [f for f in csv_files if f != merged_name]

    if not csv_files:
        print("No digitized_*.csv files found to merge.")
        return

    all_rows = []
    for filepath in csv_files:
        with open(filepath, "r") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            all_rows.extend(rows)
            print(f"  {os.path.basename(filepath)}: {len(rows)} points")

    if not all_rows:
        print("No data rows found.")
        return

    # Re-number point_id across the merged file
    for i, row in enumerate(all_rows):
        row["point_id"] = i + 1

    with open(merged_name, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=all_rows[0].keys())
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"\n  Merged {len(all_rows)} total points from {len(csv_files)} maps")
    print(f"  Saved to {merged_name}")


def list_maps(config):
    """Print all available map keys and their details."""
    print(f"\nAvailable maps ({len(config)} entries):\n")
    for key, entry in config.items():
        ref1 = entry["ref_points"][0]
        ref2 = entry["ref_points"][1]
        status = "✓" if os.path.exists(
            os.path.join(OUTPUT_DIR, f"digitized_{key}.csv")
        ) else " "
        print(f"  [{status}] {key}")
        print(f"      Source: {entry['source']}  |  Image: {entry['image']}")
        print(f"      Ref 1: {ref1['name']} ({ref1['lat']}, {ref1['lon']})")
        print(f"      Ref 2: {ref2['name']} ({ref2['lat']}, {ref2['lon']})")
        if entry.get("notes"):
            print(f"      Notes: {entry['notes']}")
        print()


def main():
    config = load_config()

    if len(sys.argv) < 2:
        print("Usage:")
        print("  python digitize_map.py <map_key>    Digitise one map")
        print("  python digitize_map.py --list        Show available maps")
        print("  python digitize_map.py --merge       Merge all CSVs")
        print()
        list_maps(config)
        return

    arg = sys.argv[1]

    if arg == "--list":
        list_maps(config)
        return

    if arg == "--merge":
        merge_all()
        return

    # Digitise a specific map
    map_key = arg
    if map_key not in config:
        print(f"Error: '{map_key}' not found in config.")
        print("Available keys:")
        for k in config:
            print(f"  {k}")
        return

    digitize_one_map(map_key, config[map_key])


if __name__ == "__main__":
    main()
