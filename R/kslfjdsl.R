# Quick diagnostic: compare background points between old and new repos
# Adjust paths to match your local setup

old_bg <- read.csv("~/repos/sudan-enm-v2/data/processed/background_points.csv")
new_bg <- read.csv("~/repos/sudan-vl-enm/data/processed/background_points.csv")

cat("Old:", nrow(old_bg), "rows | New:", nrow(new_bg), "rows\n\n")

# Spatial locations
coords_match <- all.equal(old_bg[, c("longitude", "latitude")],
                          new_bg[, c("longitude", "latitude")])
cat("Coordinates match:", coords_match, "\n")

# Year assignments
years_match <- identical(old_bg$year, new_bg$year)
cat("Year assignments match:", years_match, "\n")

if (!years_match) {
  cat("Year mismatches:", sum(old_bg$year != new_bg$year), "of", nrow(old_bg), "\n")
  cat("\nOld year distribution:\n")
  print(table(old_bg$year))
  cat("\nNew year distribution:\n")
  print(table(new_bg$year))
}