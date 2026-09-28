# ============================================================================
# 13_null_model_test.R
# Tests whether the model captures real environmental structure beyond what
# random or spatially structured presences would produce. Two variants:
# (1) Uniform null: pseudo-presences drawn at random from within the
#     ecological mask.
# (2) Accessibility-weighted null: pseudo-presences drawn in proportion to
#     accessibility (1 / sqrt(1 + travel time)), testing whether the
#     environmental signal exceeds what the spatial pattern of accessibility
#     alone would produce.
#
# Each null pool is a 50,000-point sample. Every point is assigned a year from
# the occurrence-year distribution and given year-matched covariate values,
# so pseudo-presences are measured the same way as the real presences.
# Each iteration fits MaxEnt (same feature classes and regularisation as the
# selected model) to 98 pseudo-presences against the real background, then
# scores it on the real presences.
#
# p-values use (r + 1) / (n + 1), where r = null iterations with CBI at least
# the observed CBI and n = successful iterations. Results are cached after
# the first run; delete null_model_results.rds to recompute.
#
# Inputs:  outputs/models/training_data.rds
#          outputs/models/retained_vars.rds
#          outputs/models/selected_tuning.rds
#          outputs/tables/enmeval_results.csv
#          data/raw/ecological_mask_150mm.tif
#          data/raw/ (static and year-specific covariate rasters)
#          data/raw/weiss_travel_time.tif
#          R/extract_year_matched.R
# Outputs: outputs/models/null_model_results.rds
#          outputs/figures/null_model_test.png
#          outputs/figures/null_model_test_accessible.png
#          outputs/figures/fig_null_model_panel.png
# ============================================================================

source(here::here("R", "params.R"))
source(here::here("R", "plotting_theme.R"))
source(here::here("R", "extract_year_matched.R"))

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
year_weights <- train$occ_clean |> count(year, name = "weight")
assign_years <- function(df) {
    df$year <- sample(year_weights$year, nrow(df), replace = TRUE,
                      prob = year_weights$weight)
    df
  }
  pts_to_df <- function(v) {
    as.data.frame(v, geom = "XY") |>
      rename(longitude = x, latitude = y) |>
      select(longitude, latitude)
  }


vars     <- readRDS(file.path(DIR_MODELS, "retained_vars.rds"))
eco_mask <- rast(here::here("data", "raw", "ecological_mask_150mm.tif"))

cov_stack <- rast(file.path(DIR_COVARIATES, COV_FILES))
names(cov_stack) <- names(COV_FILES)

cat("Covariates:", paste(names(cov_stack), collapse = ", "), "\n")
cat("Presences:", nrow(train$occ_env), "\n")

# Load observed CBI for comparison
tuning <- readRDS(file.path(DIR_MODELS, "selected_tuning.rds"))
best_classes <- tolower(tuning$fc)
cat("Null models use:", tuning$fc, "rm =", tuning$rm, "\n")
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

  
  # Uniform pool: 50,000 random cells in the mask, each given a year from the
  # occurrence-year distribution and year-matched values, as real presences are
  set.seed(SEED)
  pool_u_df <- spatSample(mask_binary, size = 50000, method = "random",
                          na.rm = TRUE, as.points = TRUE) |>
    pts_to_df() |>
    assign_years()
  mask_env <- extract_year_matched(pool_u_df, vars, eco_mask)
  mask_env <- mask_env[complete.cases(mask_env), ]
  cat("Uniform null pool:", nrow(mask_env), "locations\n")

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

  null_pool_df <- pts_to_df(null_pool_pts) |> assign_years()
  mask_env_acc <- extract_year_matched(null_pool_df, vars, eco_mask)
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
  
  n_ok_u <- sum(!is.na(null_results$cbi))
  n_ok_a <- sum(!is.na(null_results_acc$cbi))
  cat("Successful fits — uniform:", n_ok_u, "| accessibility:", n_ok_a, "of", N_NULL, "\n")
  if (min(n_ok_u, n_ok_a) < 0.9 * N_NULL) {
    stop("More than 10% of null fits failed — check best_classes and inputs")
  }
  
  saveRDS(list(uniform = null_results, accessibility = null_results_acc),
          null_path)
  cat("Computed and saved null model results\n")
}

# =========================== RESULTS ========================================
# Permutation p-value: (r + 1) / (n + 1), where r = null iterations whose CBI
# is at least the observed CBI, and n = null iterations that ran successfully.
# Counting the observed model as one possible outcome means the smallest
# attainable p is 1 / (n + 1), never 0.
perm_p <- function(null_cbi, obs) {
  null_cbi <- null_cbi[!is.na(null_cbi)]
  (sum(null_cbi >= obs) + 1) / (length(null_cbi) + 1)
}

p_val_uniform <- perm_p(null_results$cbi, obs_cbi)
p_val_access  <- perm_p(null_results_acc$cbi, obs_cbi)
n_uniform <- sum(!is.na(null_results$cbi))
n_access  <- sum(!is.na(null_results_acc$cbi))

cat("\n--- Uniform null distribution ---\n")
cat("CBI \u2014 mean:", round(mean(null_results$cbi, na.rm = TRUE), 3),
    " sd:", round(sd(null_results$cbi, na.rm = TRUE), 3),
    " 95th:", round(quantile(null_results$cbi, 0.95, na.rm = TRUE), 3), "\n")
cat("p-value (CBI):", round(p_val_uniform, 3), "(", n_uniform, "successful iterations )\n")

cat("\n--- Accessibility-weighted null distribution ---\n")
cat("CBI \u2014 mean:", round(mean(null_results_acc$cbi, na.rm = TRUE), 3),
    " sd:", round(sd(null_results_acc$cbi, na.rm = TRUE), 3),
    " 95th:", round(quantile(null_results_acc$cbi, 0.95, na.rm = TRUE), 3), "\n")
cat("p-value (CBI):", round(p_val_access, 3), "(", n_access, "successful iterations )\n")

p_auc_uniform <- perm_p(null_results$auc, obs_auc)
p_auc_access  <- perm_p(null_results_acc$auc, obs_auc)

cat("\n--- AUC comparison (observed AUC =", round(obs_auc, 3), ") ---\n")
cat("Uniform null       — mean:", round(mean(null_results$auc, na.rm = TRUE), 3),
    " 95th:", round(quantile(null_results$auc, 0.95, na.rm = TRUE), 3),
    " p:", round(p_auc_uniform, 3), "\n")
cat("Accessibility null — mean:", round(mean(null_results_acc$auc, na.rm = TRUE), 3),
    " 95th:", round(quantile(null_results_acc$auc, 0.95, na.rm = TRUE), 3),
    " p:", round(p_auc_access, 3), "\n")

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
        subtitle = paste0(N_NULL, " randomisations; p = ", round(p_val_uniform, 2))
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
        subtitle = paste0(N_NULL, " randomisations; p = ", round(p_val_access, 2))
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
