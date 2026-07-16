# ============================================================================
# 10_precision_sensitivity.R
# Drops the 8 occurrence records with >4 km positional uncertainty (van de
# Bogaart 2013, Hassan 2020), refits MaxEnt on the high-precision subset,
# and compares the suitability surface and ARP to the full-data model.
#
# Inputs:  outputs/models/training_data.rds
#          outputs/models/spatial_cv_folds.rds
#          outputs/models/retained_vars.rds
#          outputs/surfaces/maxent_suitability.tif
#          outputs/surfaces/worldpop_2025_aligned.tif
#          outputs/tables/arp_summary.csv
#          data/processed/occurrences_thinned.csv
# Outputs: outputs/models/maxent_precision_subset.rds
#          outputs/surfaces/maxent_suitability_precision_subset.tif
#          outputs/tables/precision_sensitivity.csv
# ============================================================================

source(here::here("R", "params.R"))

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(maxnet)
  library(ecospat)
  library(ggplot2)
})

# ------------------------------ Load inputs ---------------------------------

train   <- readRDS(file.path(DIR_MODELS, "training_data.rds"))
folds   <- readRDS(file.path(DIR_MODELS, "spatial_cv_folds.rds"))
suit_mx <- rast(file.path(DIR_SURFACES, "maxent_suitability.tif"))
vars    <- readRDS(file.path(DIR_MODELS, "retained_vars.rds"))
tuning  <- readRDS(file.path(DIR_MODELS, "selected_tuning.rds"))

best_classes <- tolower(tuning$fc)
cat("Tuning:", tuning$fc, "rm =", tuning$rm, "\n")

# ----------------------- Subset presences -----------------------------------

drop_sources <- c("vandebogaart_et_al_2013_fig1",
                  "vandebogaart_et_al_2013_fig2",
                  "vandebogaart_et_al_2013_fig3",
                  "vandebogaart_et_al_2013_fig4",
                  "vandebogaart_et_al_2013_fig5",
                  "hassan_et_al_2020_fig1",
                  "hassan_et_al_2020_fig2",
                  "hassan_et_al_2020_fig3")

drop_idx <- which(train$occ_clean$source %in% drop_sources)

cat("Dropping", length(drop_idx), "map-digitised points (>4 km accuracy)\n")
cat("Full model:", nrow(train$occ_env), "presences\n")

occ_env_sub   <- train$occ_env[-drop_idx, ]
occ_clean_sub <- train$occ_clean[-drop_idx, ]

cat("Precision subset:", nrow(occ_env_sub), "presences\n")

# ---------------------- Fold alignment --------------------------------------
# Build fold vector: start from full folds_ids (10,099), drop the row that
# MaxEnt's cell×year dedup removed, then drop the 8 precision rows.

occ_orig   <- read.csv(here::here("data", "processed", "occurrences_thinned.csv"))
occ_keys   <- paste(round(occ_orig$longitude, 5), round(occ_orig$latitude, 5), occ_orig$year)
train_keys <- paste(round(train$occ_clean$longitude, 5), round(train$occ_clean$latitude, 5), train$occ_clean$year)
dedup_idx  <- which(!occ_keys %in% train_keys)

# Remove dedup row from fold vector → length 10,098 (98 pres + 10,000 bg)
fold_ids <- folds$folds_ids[-(dedup_idx)]

# Now drop the 8 precision rows from the presence portion
fold_ids_sub <- fold_ids[-(drop_idx)]

# Combine subset presences + full background
df_sub <- bind_rows(
  bind_cols(occ_env_sub, occ_clean_sub[, c("longitude", "latitude")]) |>
    mutate(pa = 1),
  bind_cols(train$bg_env, train$bg_clean[, c("longitude", "latitude")]) |>
    mutate(pa = 0)
)

stopifnot(
  "Fold vector length doesn't match subset data" =
    length(fold_ids_sub) == nrow(df_sub)
)

df_sub$fold <- fold_ids_sub

cat("Combined:", nrow(df_sub), "rows (",
    sum(df_sub$pa == 1), "pres,", sum(df_sub$pa == 0), "bg)\n")
cat("NA folds:", sum(is.na(df_sub$fold)), "\n")
cat("Presences per fold:\n")
print(table(df_sub$fold[df_sub$pa == 1]))

# ------------------------ Refit MaxEnt --------------------------------------

p_sub <- as.matrix(occ_env_sub[, vars])
b_sub <- as.matrix(train$bg_env[, vars])

mod_sub <- maxnet(
  p    = c(rep(1, nrow(p_sub)), rep(0, nrow(b_sub))),
  data = as.data.frame(rbind(p_sub, b_sub)),
  f    = maxnet.formula(
    p    = c(rep(1, nrow(p_sub)), rep(0, nrow(b_sub))),
    data = as.data.frame(rbind(p_sub, b_sub)),
    classes = best_classes
  ),
  regmult = tuning$rm
)

cat("Refit model: ", sum(mod_sub$betas != 0), "non-zero /",
    length(mod_sub$betas), "total coefficients\n")

# ---------------------- Spatial CV evaluation --------------------------------

fold_cbi <- fold_auc <- numeric(4)

for (k in 1:4) {
  idx_train <- df_sub$fold != k
  idx_test  <- df_sub$fold == k

  p_tr <- as.matrix(df_sub[idx_train & df_sub$pa == 1, vars])
  b_tr <- as.matrix(df_sub[idx_train & df_sub$pa == 0, vars])

  mod_k <- maxnet(
    p    = c(rep(1, nrow(p_tr)), rep(0, nrow(b_tr))),
    data = as.data.frame(rbind(p_tr, b_tr)),
    f    = maxnet.formula(
      p    = c(rep(1, nrow(p_tr)), rep(0, nrow(b_tr))),
      data = as.data.frame(rbind(p_tr, b_tr)),
      classes = best_classes
    ),
    regmult = tuning$rm
  )

  test_data <- df_sub[idx_test, vars]
  pred_test <- predict(mod_k, newdata = test_data, type = "cloglog")[, 1]

  pres_pred <- pred_test[df_sub$pa[idx_test] == 1]
  bg_pred   <- pred_test[df_sub$pa[idx_test] == 0]

  boyce <- ecospat.boyce(fit = pred_test, obs = pres_pred,
                         nclass = 0, PEplot = FALSE)
  fold_cbi[k] <- boyce$cor

  n1 <- length(pres_pred)
  n0 <- length(bg_pred)
  fold_auc[k] <- (sum(rank(c(pres_pred, bg_pred))[1:n1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

cat("\nSpatial CV — precision subset:\n")
cat("CBI:", round(mean(fold_cbi, na.rm = TRUE), 3),
    "\u00b1", round(sd(fold_cbi, na.rm = TRUE), 3), "\n")
cat("AUC:", round(mean(fold_auc), 3),
    "\u00b1", round(sd(fold_auc), 3), "\n")

# ---------------------- Predict and compare ---------------------------------

suit_sub_path <- file.path(DIR_SURFACES, "maxent_suitability_precision_subset.tif")

if (file.exists(suit_sub_path)) {
  cat("Loading cached precision subset surface\n")
  suit_sub <- rast(suit_sub_path)
} else {
  cov_stack <- rast(file.path(DIR_COVARIATES, COV_FILES[vars]))
  names(cov_stack) <- vars

  suit_sub <- predict(cov_stack, mod_sub, type = "cloglog",
                      clamp = TRUE, na.rm = TRUE)

  writeRaster(suit_sub, suit_sub_path, overwrite = TRUE)
  cat("Computed and saved precision subset surface\n")
}

# Surface correlation
set.seed(SEED)
valid_cells <- which(!is.na(values(suit_mx)) & !is.na(values(suit_sub)))
samp_idx <- sample(valid_cells, min(50000, length(valid_cells)))
r_pearson <- round(cor(values(suit_mx)[samp_idx],
                       values(suit_sub)[samp_idx]), 3)

cat("\nSurface correlation (50k sample): Pearson =", r_pearson, "\n")

# ----------------------- ARP comparison -------------------------------------

pop_aligned <- rast(file.path(DIR_SURFACES, "worldpop_2025_aligned.tif"))
total_pop   <- global(pop_aligned, "sum", na.rm = TRUE)[[1]]

# Subset thresholds
pred_occ_sub <- predict(mod_sub, newdata = occ_env_sub[, vars],
                        type = "cloglog")[, 1]
pred_bg_sub  <- predict(mod_sub, newdata = train$bg_env[, vars],
                        type = "cloglog")[, 1]

sub_p10 <- unname(quantile(pred_occ_sub, 0.10))

candidates <- sort(unique(c(pred_occ_sub, pred_bg_sub)))
sens <- sapply(candidates, function(t) mean(pred_occ_sub >= t))
spec <- sapply(candidates, function(t) mean(pred_bg_sub < t))
sub_maxsss <- candidates[which.max(sens + spec)]

cat("Subset p10:   ", round(sub_p10, 4), "\n")
cat("Subset maxSSS:", round(sub_maxsss, 4), "\n")

# ARP estimates
sub_arp_p10      <- global(pop_aligned * (suit_sub >= sub_p10), "sum", na.rm = TRUE)[[1]]
sub_arp_maxsss   <- global(pop_aligned * (suit_sub >= sub_maxsss), "sum", na.rm = TRUE)[[1]]
sub_arp_weighted <- global(pop_aligned * suit_sub, "sum", na.rm = TRUE)[[1]]

# Load full-model ARP for comparison
mx_arp       <- read.csv(file.path(DIR_TABLES, "arp_summary.csv"))
full_weighted <- mx_arp$arp[mx_arp$metric == "risk_weighted"]
full_maxsss   <- mx_arp$arp[mx_arp$metric == "maxSSS"]
full_p10      <- mx_arp$arp[mx_arp$metric == "p10"]

cat("\n--- ARP comparison ---\n")
cat(sprintf("%-20s %14s %14s\n", "", "Full (n=98)", "Subset (n=90)"))
cat(sprintf("%-20s %14s %14s\n", "Risk-weighted",
    format(round(full_weighted), big.mark = ","),
    format(round(sub_arp_weighted), big.mark = ",")))
cat(sprintf("%-20s %14s %14s\n", "maxSSS binary",
    format(round(full_maxsss), big.mark = ","),
    format(round(sub_arp_maxsss), big.mark = ",")))
cat(sprintf("%-20s %14s %14s\n", "p10 binary",
    format(round(full_p10), big.mark = ","),
    format(round(sub_arp_p10), big.mark = ",")))

pct_shift <- round(100 * abs(sub_arp_weighted - full_weighted) / full_weighted, 1)
cat("Risk-weighted ARP shift:", pct_shift, "%\n")

# --------------------------------- Save -------------------------------------

saveRDS(mod_sub, file.path(DIR_MODELS, "maxent_precision_subset.rds"))

precision_summary <- data.frame(
  metric     = c("Presences", "Coefficients", "CBI (spatial CV)", "AUC (spatial CV)",
                  "Mean suitability", "ARP risk-weighted", "ARP maxSSS binary",
                  "ARP p10 binary", "p10 threshold", "maxSSS threshold",
                  "Surface correlation (Pearson)"),
  full_model = c(98, 13,
                 0.873, 0.849,
                 round(global(suit_mx, "mean", na.rm = TRUE)[[1]], 4),
                 round(full_weighted), round(full_maxsss), round(full_p10),
                 round(mx_arp$threshold[mx_arp$metric == "p10"], 4),
                 round(mx_arp$threshold[mx_arp$metric == "maxSSS"], 4),
                 NA),
  subset     = c(90, sum(mod_sub$betas != 0),
                 round(mean(fold_cbi, na.rm = TRUE), 3),
                 round(mean(fold_auc), 3),
                 round(global(suit_sub, "mean", na.rm = TRUE)[[1]], 4),
                 round(sub_arp_weighted), round(sub_arp_maxsss), round(sub_arp_p10),
                 round(sub_p10, 4), round(sub_maxsss, 4),
                 r_pearson)
)

write.csv(precision_summary, file.path(DIR_TABLES, "precision_sensitivity.csv"),
          row.names = FALSE)

cat("\nSaved:\n")
cat("  ", file.path(DIR_MODELS, "maxent_precision_subset.rds"), "\n")
cat("  ", suit_sub_path, "\n")
cat("  ", file.path(DIR_TABLES, "precision_sensitivity.csv"), "\n")

cat("\n10_precision_sensitivity.R complete\n")