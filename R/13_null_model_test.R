# ============================================================================
# 13_null_model_test.R
# Tests whether the model captures real environmental structure beyond what
# random or spatially structured presences would produce. Two variants:
# (1) Uniform null: random pseudo-presences from the ecological mask.
# (2) Accessibility-weighted null: pseudo-presences drawn proportional to
#     geographic accessibility, testing whether environmental signal exceeds
#     spatial pattern of accessibility alone.
#
# Both are cached after first run.
#
# Inputs:  outputs/models/training_data.rds
#          outputs/models/retained_vars.rds
#          data/processed/ecological_mask.tif
#          data/raw/ (covariate rasters)
#          data/raw/weiss_travel_time.tif
# Outputs: outputs/models/null_model_results.rds
#          outputs/figures/null_model_test.png
#          outputs/figures/null_model_test_accessible.png
#.         outputs/figures/fig_null_model_panel.png
# ============================================================================

source(here::here("R", "params.R"))
source(here::here("R", "plotting_theme.R"))

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(maxnet)
  library(ggplot2)
  library(ecospat)
  library(patchwork)
})

# ------------------------------ Load inputs ---------------------------------

train    <- readRDS(file.path(DIR_MODELS, "training_data.rds"))
vars     <- readRDS(file.path(DIR_MODELS, "retained_vars.rds"))
eco_mask <- rast(here::here("data", "processed", "ecological_mask.tif"))

cov_stack <- rast(file.path(DIR_COVARIATES, COV_FILES))
names(cov_stack) <- names(COV_FILES)

cat("Covariates:", paste(names(cov_stack), collapse = ", "), "\n")
cat("Presences:", nrow(train$occ_env), "\n")

# Load observed CBI for comparison
tuning <- readRDS(file.path(DIR_MODELS, "selected_tuning.rds"))
res_table <- read.csv(file.path(DIR_TABLES, "enmeval_results.csv"))
obs_row <- res_table |> filter(fc == tuning$fc, rm == tuning$rm)
obs_cbi <- obs_row$cbi.val.avg
obs_auc <- obs_row$auc.val.avg

cat("Observed CBI:", round(obs_cbi, 3), "| AUC:", round(obs_auc, 3), "\n")

# ======================== UNIFORM NULL MODEL ================================

null_path <- file.path(DIR_MODELS, "null_model_results.rds")

if (file.exists(null_path)) {
  cat("Loading cached null model results\n")
  null_cache <- readRDS(null_path)
  null_results     <- null_cache$uniform
  null_results_acc <- null_cache$accessibility
} else {

  # --------------------- Build uniform null pool ----------------------------

  mask_binary <- subst(eco_mask, 0, NA)
  mask_cells  <- which(values(eco_mask) == 1)

  mask_env <- cov_stack[mask_cells]
  complete <- complete.cases(mask_env)
  mask_cells <- mask_cells[complete]
  mask_env   <- mask_env[complete, ]
  cat("Uniform null pool:", nrow(mask_env), "cells\n")

  # ---------------------- Uniform null loop --------------------------------

  set.seed(SEED)

  occ_env <- as.matrix(train$occ_env[, vars])
  bg_env  <- as.matrix(train$bg_env[, vars])

  cat("Running", N_NULL, "uniform null iterations...\n")

  null_results <- data.frame(
    iter = integer(N_NULL),
    cbi  = numeric(N_NULL),
    auc  = numeric(N_NULL)
  )

  for (i in seq_len(N_NULL)) {
    set.seed(SEED + i)

    idx <- sample(nrow(mask_env), N_PRES, replace = FALSE)
    null_train <- mask_env[idx, ]

    p_vec <- c(rep(1, N_PRES), rep(0, nrow(train$bg_env)))
    all_data <- as.data.frame(rbind(
      as.matrix(null_train[, vars]),
      bg_env
    ))

    mod_null <- tryCatch(
      maxnet(
        p    = p_vec,
        data = all_data,
        f    = maxnet.formula(p = p_vec, data = all_data, classes = best_classes),
        regmult = tuning$rm
      ),
      error = function(e) NULL
    )

    if (is.null(mod_null)) {
      null_results$iter[i] <- i
      null_results$cbi[i]  <- NA
      null_results$auc[i]  <- NA
      next
    }

    pred_occ <- predict(mod_null, as.data.frame(occ_env), type = "cloglog")
    pred_bg  <- predict(mod_null, as.data.frame(bg_env),  type = "cloglog")

    boyce <- tryCatch(
      ecospat.boyce(fit = pred_bg, obs = pred_occ, PEplot = FALSE),
      error = function(e) list(cor = NA)
    )

    labels  <- c(rep(1, length(pred_occ)), rep(0, length(pred_bg)))
    preds   <- c(pred_occ, pred_bg)
    n1      <- sum(labels == 1)
    n0      <- sum(labels == 0)
    ranks   <- rank(preds)
    auc_val <- (sum(ranks[labels == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)

    null_results$iter[i] <- i
    null_results$cbi[i]  <- boyce$cor
    null_results$auc[i]  <- auc_val

    if (i %% 25 == 0) cat("  Completed", i, "of", N_NULL, "\n")
  }

  # ================ ACCESSIBILITY-WEIGHTED NULL MODEL =======================

  tt_raw <- rast(here::here("data", "raw", "weiss_travel_time.tif"))
  bias_r <- 1 / sqrt(1 + resample(tt_raw, eco_mask, method = "bilinear"))
  bias_r <- mask(bias_r, mask_binary)

  set.seed(SEED)
  null_pool_pts <- spatSample(bias_r, size = 50000, method = "weights",
                              na.rm = TRUE, as.points = TRUE)

  mask_env_acc <- terra::extract(cov_stack, null_pool_pts, ID = FALSE)
  complete <- complete.cases(mask_env_acc)
  mask_env_acc <- mask_env_acc[complete, ]
  cat("Accessibility-weighted null pool:", nrow(mask_env_acc), "locations\n")

  # ------------------- Accessibility null loop ------------------------------

  set.seed(SEED)

  cat("Running", N_NULL, "accessibility-weighted null iterations...\n")

  null_results_acc <- data.frame(
    iter = integer(N_NULL),
    cbi  = numeric(N_NULL),
    auc  = numeric(N_NULL)
  )

  for (i in seq_len(N_NULL)) {
    set.seed(SEED + i)

    idx <- sample(nrow(mask_env_acc), N_PRES, replace = FALSE)
    null_train <- mask_env_acc[idx, ]

    p_vec <- c(rep(1, N_PRES), rep(0, nrow(train$bg_env)))
    all_data <- as.data.frame(rbind(
      as.matrix(null_train[, vars]),
      bg_env
    ))

    mod_null <- tryCatch(
      maxnet(
        p    = p_vec,
        data = all_data,
        f    = maxnet.formula(p = p_vec, data = all_data, classes = best_classes),
        regmult = tuning$rm
      ),
      error = function(e) NULL
    )

    if (is.null(mod_null)) {
      null_results_acc$iter[i] <- i
      null_results_acc$cbi[i]  <- NA
      null_results_acc$auc[i]  <- NA
      next
    }

    pred_occ <- predict(mod_null, as.data.frame(occ_env), type = "cloglog")
    pred_bg  <- predict(mod_null, as.data.frame(bg_env),  type = "cloglog")

    boyce <- tryCatch(
      ecospat.boyce(fit = pred_bg, obs = pred_occ, PEplot = FALSE),
      error = function(e) list(cor = NA)
    )

    labels  <- c(rep(1, length(pred_occ)), rep(0, length(pred_bg)))
    preds   <- c(pred_occ, pred_bg)
    n1      <- sum(labels == 1)
    n0      <- sum(labels == 0)
    ranks   <- rank(preds)
    auc_val <- (sum(ranks[labels == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)

    null_results_acc$iter[i] <- i
    null_results_acc$cbi[i]  <- boyce$cor
    null_results_acc$auc[i]  <- auc_val

    if (i %% 25 == 0) cat("  Completed", i, "of", N_NULL, "\n")
  }

  # Save both sets together
  saveRDS(list(uniform = null_results, accessibility = null_results_acc),
          null_path)
  cat("Computed and saved null model results\n")
}

# =========================== RESULTS ========================================

cat("\n--- Uniform null distribution ---\n")
cat("CBI \u2014 mean:", round(mean(null_results$cbi, na.rm = TRUE), 3),
    " sd:", round(sd(null_results$cbi, na.rm = TRUE), 3),
    " 95th:", round(quantile(null_results$cbi, 0.95, na.rm = TRUE), 3), "\n")
cat("p-value (CBI):", round(mean(null_results$cbi >= obs_cbi, na.rm = TRUE), 3), "\n")

cat("\n--- Accessibility-weighted null distribution ---\n")
cat("CBI \u2014 mean:", round(mean(null_results_acc$cbi, na.rm = TRUE), 3),
    " sd:", round(sd(null_results_acc$cbi, na.rm = TRUE), 3),
    " 95th:", round(quantile(null_results_acc$cbi, 0.95, na.rm = TRUE), 3), "\n")
cat("p-value (CBI):", round(mean(null_results_acc$cbi >= obs_cbi, na.rm = TRUE), 3), "\n")

# =========================== FIGURES ========================================

p_null <- ggplot(null_results, aes(x = cbi)) +
  geom_histogram(binwidth = 0.05, fill = "grey70", colour = "white") +
  geom_vline(xintercept = obs_cbi, colour = "red", linewidth = 1,
             linetype = "dashed") +
  annotate("text", x = obs_cbi - 0.03, y = Inf, vjust = 2, hjust = 1,
           label = paste0("Observed\nCBI = ", round(obs_cbi, 3)),
           colour = "red", size = 3.5, fontface = "bold") +
  labs(
    x = "Continuous Boyce Index (CBI)",
    y = "Count",
    title = "Null model significance test (uniform)",
    subtitle = paste0("99 randomisations; observed CBI exceeds null 95th ",
                      "percentile (", round(quantile(null_results$cbi, 0.95,
                      na.rm = TRUE), 3), ")")
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(DIR_FIGS, "null_model_test.png"), p_null,
       width = 7, height = 5, dpi = 300, bg = "white")
cat("Saved null_model_test.png\n")

p_null_acc <- ggplot(null_results_acc, aes(x = cbi)) +
  geom_histogram(binwidth = 0.05, fill = "grey70", colour = "white") +
  geom_vline(xintercept = obs_cbi, colour = "red", linewidth = 1,
             linetype = "dashed") +
  annotate("text", x = obs_cbi - 0.03, y = Inf, vjust = 2, hjust = 1,
           label = paste0("Observed\nCBI = ", round(obs_cbi, 3)),
           colour = "red", size = 3.5, fontface = "bold") +
  labs(
    x = "Continuous Boyce Index (CBI)",
    y = "Count",
    title = "Null model significance test (accessibility-weighted)",
    subtitle = paste0("99 randomisations; p = ",
                      round(mean(null_results_acc$cbi >= obs_cbi, na.rm = TRUE), 2))
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(DIR_FIGS, "null_model_test_accessible.png"), p_null_acc,
       width = 7, height = 5, dpi = 300, bg = "white")
cat("Saved null_model_test_accessible.png\n")


# ---------------- PANEL FIGURE -------------

# Shared x-axis so both distributions are visually comparable
x_limits <- c(-1, 1.05)

p_uniform <- ggplot(null_results, aes(x = cbi)) +
  geom_histogram(binwidth = 0.05, colour = "white") +
  geom_vline(xintercept = obs_cbi, colour = "red", linewidth = 0.8,
             linetype = "dashed") +
  annotate("text", x = obs_cbi - 0.03, y = Inf, vjust = 2, hjust = 1,
           label = paste0("Observed\nCBI = ", round(obs_cbi, 3)),
           colour = "red", size = 3, fontface = "bold") +
  coord_cartesian(xlim = x_limits) +
  labs(x = "Continuous Boyce Index (CBI)", y = "Count",
       title = "(a) Uniform null")

p_accessible <- ggplot(null_results_acc, aes(x = cbi)) +
  geom_histogram(binwidth = 0.05, colour = "white") +
  geom_vline(xintercept = obs_cbi, colour = "red", linewidth = 0.8,
             linetype = "dashed") +
  annotate("text", x = obs_cbi - 0.03, y = Inf, vjust = 2, hjust = 1,
           label = paste0("Observed\nCBI = ", round(obs_cbi, 3)),
           colour = "red", size = 3, fontface = "bold") +
  coord_cartesian(xlim = x_limits) +
  labs(x = "Continuous Boyce Index (CBI)", y = "Count",
       title = "(b) Accessibility-weighted null")

p_null_panel <- p_uniform + p_accessible

ggsave(file.path(DIR_FIGS, "fig_null_model_panel.png"), p_null_panel,
         width = FIG_WIDTH_FULL, height = 8)

cat("Saved fig_null_model_panel.png\n")

cat("\n13_null_model_test.R complete\n")
