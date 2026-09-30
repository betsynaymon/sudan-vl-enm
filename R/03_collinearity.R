# ============================================================================
# 03_collinearity.R
# Screens the full candidate covariate set for collinearity, selects the final
# predictor set on ecological grounds, and writes the retained variable names
# to disk for all downstream scripts. Screening uses the prediction-surface
# layers (static covariates + 2000–2024 long-term means).
#
# Inputs:  COV_FILES (9 candidates), DOMAIN_FILE
# Outputs: outputs/models/retained_vars.rds
#          outputs/models/correlation_matrix.rds (read by 15 for labels)
#          outputs/figures/collinearity_correlation_matrix.png
#          outputs/figures/collinearity_correlation_matrix.pdf
# ============================================================================

source(here::here("R", "params.R"))

suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
})

set.seed(SEED)

# ----------------------- Load candidate covariates --------------------------
covs <- rast(file.path(DIR_COVARIATES, COV_FILES))
names(covs) <- names(COV_FILES)

# Screen over the calibration domain: collinearity matters where the model
# is fitted, and the background covers all of Sudan. Values are read into
# memory; sampling rows is far faster than spatSample() on a file-backed stack.
in_domain <- values(rast(DOMAIN_FILE), mat = FALSE) == 1
vals <- values(covs)[in_domain, ]

cat("Non-NA domain cells per layer:\n")
print(colSums(!is.na(vals)))

# ----------------------- Sample and correlate -------------------------------
# Region >= 150 mm within the domain. The mask is not clipped to Sudan,
# so it is always intersected with the domain.
sens_dom <- values(rast(SENS_MASK_FILE), mat = FALSE)[in_domain] %in% 1
ok <- complete.cases(vals)

set.seed(SEED)
sample_vals <- as.data.frame(vals[ok, ][sample(sum(ok), COLLIN_SAMPLE_N), ])
set.seed(SEED)
wet_vals <- vals[ok & sens_dom, ][sample(sum(ok & sens_dom), COLLIN_SAMPLE_N), ]
rm(vals)
cor_wet <- cor(wet_vals)
saveRDS(cor_wet, file.path(DIR_MODELS, "correlation_matrix_150mm.rds"))
cor_mat <- cor(sample_vals, use = "complete.obs", method = "pearson")

saveRDS(cor_mat, file.path(DIR_MODELS, "correlation_matrix.rds"))

cat("Sample:", nrow(sample_vals), "cells\n")
cat("\nPearson correlation matrix:\n")
print(round(cor_mat, 2))

# Flag pairs exceeding threshold
high_pairs <- which(abs(cor_mat) >= COR_THRESHOLD & upper.tri(cor_mat), arr.ind = TRUE)
if (nrow(high_pairs) > 0) {
  cat("\nPairs exceeding |r| >=", COR_THRESHOLD, ":\n")
  for (i in seq_len(nrow(high_pairs))) {
    r <- high_pairs[i, ]
    cat("  ", rownames(cor_mat)[r[1]], " — ", colnames(cor_mat)[r[2]],
        ": ", round(cor_mat[r[1], r[2]], 2), "\n")
  }
}

# ------------------------ VIF candidate sets --------------------------------
# Candidate sets resolve the greenness-wetness group (NDVI, tree cover,
# rainfall) and the temperature-elevation group (LST day, LST night,
# elevation) differently.

candidate_sets <- list(
  A_rain         = c("slope", "river_dist", "vertisols", "lst_night", "rainfall"),
  B_ndvi         = c("slope", "river_dist", "vertisols", "lst_night", "ndvi"),
  C_rain_daytemp = c("slope", "river_dist", "vertisols", "lst_night", "lst_day", "rainfall"),
  D_rain_tree    = c("slope", "river_dist", "vertisols", "lst_night", "rainfall", "treecover")
)

for (set_name in names(candidate_sets)) {
  cat("\n---", set_name, "---\n")
  print(round(diag(solve(cor_mat[candidate_sets[[set_name]], candidate_sets[[set_name]]])), 2))
}

# ---------------------- Finalize variable set -------------------------------
# Set A, selected on ecological grounds:
#   slope       topographic steepness
#   river_dist  distance to drainage, including seasonal khors
#   vertisols   black-cotton soil
#   lst_night   night temperature, the conditions the nocturnal vector experiences
#   rainfall    the moisture-vegetation axis (most stable single measure)
# Dropped as redundant with a retained variable (evidence printed below):
#   elevation   distal proxy; acts through temperature
#   lst_day     daytime thermal window
#   ndvi        conflates crop and woodland
#   treecover   forest-built metric, near-zero in the endemic belt

retained_vars <- c("slope", "river_dist", "vertisols", "lst_night", "rainfall")

cat("\nStrongest correlation of each dropped variable with a retained one:\n")
for (v in setdiff(names(COV_FILES), retained_vars)) {
  r <- cor_mat[v, retained_vars]
  j <- which.max(abs(r))
  cat(sprintf("  %-10s %-10s r = %5.2f\n", v, retained_vars[j], r[j]))
}

cat("\nSame, within the >= 150 mm region of the domain:\n")
for (v in setdiff(names(COV_FILES), retained_vars)) {
  r <- cor_wet[v, retained_vars]
  j <- which.max(abs(r))
  cat(sprintf("  %-10s %-10s r = %5.2f\n", v, retained_vars[j], r[j]))
}

r_ret <- cor_mat[retained_vars, retained_vars]
diag(r_ret) <- NA
vif_ret <- diag(solve(cor_mat[retained_vars, retained_vars]))
cat("\nFinal variable set VIF:\n")
print(round(vif_ret, 2))
cat("Max |r| among retained:", round(max(abs(r_ret), na.rm = TRUE), 2), "\n")

# Stop before writing retained_vars.rds if the set fails either screen:
# the selection then needs reassessing, not propagating downstream.
stopifnot(
  "A retained pair exceeds COR_THRESHOLD" = max(abs(r_ret), na.rm = TRUE) < COR_THRESHOLD,
  "A retained VIF exceeds VIF_THRESHOLD"  = max(vif_ret) < VIF_THRESHOLD
)

saveRDS(retained_vars, file.path(DIR_MODELS, "retained_vars.rds"))
cat("Saved retained_vars.rds\n")

# -------------------- Correlation matrix figure -----------------------------

var_order  <- c("rainfall", "ndvi", "treecover",
                "lst_day", "lst_night", "elevation",
                "slope", "river_dist", "vertisols")
var_labels <- c(rainfall = "Rainfall", ndvi = "NDVI", treecover = "Tree cover",
                lst_day = "Day LST", lst_night = "Night LST", elevation = "Elevation",
                slope = "Slope", river_dist = "Distance to river", vertisols = "Vertisols")

cor_mat <- cor_mat[var_order, var_order]
cor_df <- as.data.frame(as.table(cor_mat))
names(cor_df) <- c("var1", "var2", "r")
cor_df$var1 <- factor(cor_df$var1, levels = var_order)
cor_df$var2 <- factor(cor_df$var2, levels = var_order)
cor_df <- cor_df[as.integer(cor_df$var1) >= as.integer(cor_df$var2), ]
cor_df$flag <- abs(cor_df$r) >= COR_THRESHOLD & cor_df$var1 != cor_df$var2

p <- ggplot(cor_df, aes(var2, var1, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_tile(data = cor_df[cor_df$flag, ], fill = NA,
            colour = "grey10", linewidth = 1.1) +
  geom_text(aes(label = sprintf("%.2f", r),
                colour = abs(r) > 0.5), size = 3.2) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey20"),
                      guide = "none") +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
                       midpoint = 0, limits = c(-1, 1),
                       breaks = seq(-1, 1, 0.2), name = "Pearson r") +
  scale_x_discrete(limits = var_order, labels = var_labels) +
  scale_y_discrete(limits = rev(var_order), labels = var_labels) +
  coord_fixed() + labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.height = unit(1.6, "cm"))

ggsave(file.path(DIR_FIGS, "collinearity_correlation_matrix.pdf"),
       p, width = 7.5, height = 6.5)
ggsave(file.path(DIR_FIGS, "collinearity_correlation_matrix.png"),
       p, width = 7.5, height = 6.5, dpi = 300)
cat("Saved correlation matrix figure\n")

cat("03_collinearity.R complete\n")