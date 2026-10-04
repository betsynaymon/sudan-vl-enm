# ============================================================================
# 13_null_model_test.R
# Does the model carry environmental signal beyond what random or merely
# accessible locations would produce? Null presences are drawn (1) uniformly
# over the domain, a sanity check, and (2) in proportion to (1 + travel
# time)^-b, with b set so the null presences are as close to cities as the
# real ones (median travel time). Each null presence gets a year from the
# presence-year distribution, year-matched values, and a fold from 05's
# blocks. Each iteration is cross-validated exactly as in 06: trained on null
# presences outside a fold, scored on the real presences and background
# inside it, with the selected configuration (Bohl et al. 2019). The observed
# statistic is 06's cross-validated CBI, reproduced here through the same code.
#
# p = (r + 1) / (n + 1); r = iterations with CBI >= observed, n = iterations run.
#
# Inputs:  TRAIN_FILE, TUNING_FILE, FOLDS_FILE, DOMAIN_FILE, TT_FILE,
#          retained_vars.rds, covariate rasters (COV_FILES, COV_ANNUAL)
# Outputs: outputs/tables/null_model_iterations.csv
#          outputs/tables/null_model_summary.csv
#          outputs/figures/fig_null_model_panel.png
# ============================================================================

source(here::here("R", "params.R"))
source(here::here("R", "helpers.R"))
source(here::here("R", "plotting_theme.R"))

suppressPackageStartupMessages({
  library(terra); library(dplyr); library(maxnet); library(ggplot2); library(patchwork)
})

# ------------------------------ Load inputs ---------------------------------

train   <- readRDS(TRAIN_FILE)
tuning  <- readRDS(TUNING_FILE)
vars    <- readRDS(file.path(DIR_MODELS, "retained_vars.rds"))
occ     <- train$occ_clean
bg      <- train$bg_clean
occ_env <- train$occ_env[, vars]
bg_env  <- train$bg_env[, vars]
n_pres  <- nrow(occ)
cat("Model:", tuning$fc, "rm =", tuning$rm, "| presences:", n_pres,
    "| background:", nrow(bg), "\n")

stopifnot(
  "Block folds disagree with the presences' folds (05)"  = all(fold_from_blocks(occ) == occ$fold),
  "Block folds disagree with the background's folds (05)" = all(fold_from_blocks(bg) == bg$fold)
)
# ------------------------------- Scoring ------------------------------------
# As in 06: for each fold, fit on the presences outside it plus the background
# outside it; score on the real presences and background inside it.

score_cv <- function(p_env, p_fold) {
  bind_rows(lapply(seq_len(K_FOLDS), function(k) {
    tr_p <- p_env[p_fold != k, , drop = FALSE]
    tr_b <- bg_env[bg$fold != k, ]
    m <- fit_maxnet(c(rep(1, nrow(tr_p)), rep(0, nrow(tr_b))), rbind(tr_p, tr_b),
                    tuning$fc, tuning$rm)
    eval_fold(m, test_occ = occ_env[occ$fold == k, ],
              test_bg = bg_env[bg$fold == k, ], train_occ = tr_p)
  }))
}

obs <- score_cv(occ_env, occ$fold)
stopifnot("Observed CV CBI differs from 06" = abs(mean(obs$cbi) - tuning$cbi) < 1e-8,
          "Observed CV AUC differs from 06" = abs(mean(obs$auc) - tuning$auc) < 1e-8)
obs_cbi <- tuning$cbi
obs_auc <- tuning$auc
cat("Observed CV CBI:", round(obs_cbi, 3), "| AUC:", round(obs_auc, 3), "(reproduced from 06)\n")

# --------------------------------- Pools ------------------------------------
# Candidate cells: domain cells with every long-term covariate (the
# background's frame, 04). Accessibility weight (1 + tt)^-b, with b set so the
# weighted median travel time equals the presences' (rule introduced during
# analysis, fixed before the null results; disclose). Cells without travel
# time get weight 0 (no fill; 11). b = 0.5 was the dissertation's weighting.

covs  <- domain_covs(vars)
cand  <- !any(is.na(covs))
cand  <- mask(cand, cand, maskvalues = 0)
tt    <- travel_time(cand)
cells <- which(!is.na(values(cand, mat = FALSE)))
tt_c  <- values(tt, mat = FALSE)[cells]
cat("\nCandidate cells:", length(cells), "| without travel time:", sum(is.na(tt_c)), "\n")

has_tt <- !is.na(tt_c)
occ_tt <- terra::extract(tt, as.matrix(occ[, c("longitude", "latitude")]))[, 1]
target <- median(occ_tt)
b      <- tt_exponent(tt_c, target)
cat("Presence median travel time:", round(target), "min | b =", b,
    "| weighted median at b:", round(tt_wmedian(tt_c, b)),
    "| at b = 0.5:", round(tt_wmedian(tt_c, 0.5)), "\n")

# Pools drawn with draw_points() (helpers.R).
year_w <- table(occ$year)
pools <- list(
  uniform       = draw_points(cand, cells, N_NULL_POOL, NULL, year_w, vars, SEED),
  accessibility = draw_points(cand, cells, N_NULL_POOL, ifelse(has_tt, (1 + tt_c)^(-b), 0),
                              year_w, vars, SEED + 1)
)
for (nm in names(pools)) pools[[nm]]$tt <- values(tt, mat = FALSE)[pools[[nm]]$cell]

q3 <- function(x) paste(round(quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)), collapse = " / ")
cat("\nTravel time q25 / median / q75 (min). Presences:", q3(occ_tt), "\n")
for (nm in names(pools)) {
  pl <- pools[[nm]]
  cat(nm, "pool:", nrow(pl$env), "points (", pl$dropped, "dropped ) |", q3(pl$tt),
      "| P[presence closer]:", round(auc_ties(-occ_tt, -pl$tt[!is.na(pl$tt)]), 3),
      "| per fold:", paste(table(factor(pl$fold, levels = seq_len(K_FOLDS))), collapse = " / "), "\n")
}

# ------------------------------ Null loops ----------------------------------
# Each iteration draws n_pres null presences from a pool (without replacement)
# under seed SEED + i. An iteration counts only if all four folds were scored.

run_null <- function(pl, label) {
  t0 <- Sys.time()
  bind_rows(lapply(seq_len(N_NULL), function(i) {
    set.seed(SEED + i)
    idx <- sample(nrow(pl$env), n_pres)
    r <- tryCatch(score_cv(pl$env[idx, ], pl$fold[idx]), error = function(e) NULL)
    if (i == 1) cat("\n", label, ": first iteration took ",
                    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
                    " min; expect about ", N_NULL, " times that\n", sep = "")
    if (i %% 10 == 0) cat("  ", label, i, "of", N_NULL, "\n")
    data.frame(null = label, iter = i,
               cbi = if (is.null(r) || anyNA(r$cbi)) NA_real_ else mean(r$cbi),
               auc = if (is.null(r)) NA_real_ else mean(r$auc))
  }))
}
iters <- bind_rows(run_null(pools$uniform, "uniform"),
                   run_null(pools$accessibility, "accessibility"))

# ------------------------------- Results ------------------------------------
# Reading fixed before the run: the uniform null is a sanity check (expected
# p = 1 / (N_NULL + 1)). The accessibility null is the test: p <= NULL_ALPHA
# means the model carries environmental signal beyond what equally accessible
# locations produce; otherwise that claim fails.

perm_p <- function(x, obs) (sum(x >= obs) + 1) / (length(x) + 1)

summary_df <- iters |> filter(!is.na(cbi)) |> group_by(null) |>
  summarise(n_ok = n(), cbi_mean = mean(cbi), cbi_sd = sd(cbi),
            cbi_q95 = quantile(cbi, 0.95), cbi_max = max(cbi),
            p_cbi = perm_p(cbi, obs_cbi), z_cbi = (obs_cbi - mean(cbi)) / sd(cbi),
            auc_mean = mean(auc), p_auc = perm_p(auc, obs_auc), .groups = "drop") |>
  mutate(b = ifelse(null == "accessibility", b, NA_real_),
         beaten = p_cbi <= NULL_ALPHA)
stopifnot("A null has no successful iteration"      = nrow(summary_df) == 2,
          "More than 10% of null iterations failed" = all(summary_df$n_ok >= 0.9 * N_NULL))

cat("\nObserved CV CBI", round(obs_cbi, 3), "against each null:\n")
summary_df |> mutate(across(where(is.double), ~ round(., 3))) |> print(width = Inf)

# -------------------------------- Figure ------------------------------------

p_hist <- function(nm, title) {
  ggplot(filter(iters, null == nm, !is.na(cbi)), aes(cbi)) +
    geom_histogram(binwidth = 0.05, boundary = 0, fill = "grey70", colour = "white") +
    geom_vline(xintercept = obs_cbi, colour = col_observed, linewidth = 0.8, linetype = "dashed") +
    annotate("text", x = obs_cbi, y = Inf, vjust = 1.5, hjust = 1.05, size = 3,
             colour = col_observed, label = paste0("Observed\nCBI = ", round(obs_cbi, 3))) +
    coord_cartesian(xlim = c(-1, 1)) +
    labs(x = "Cross-validated CBI", y = "Iterations", title = title,
         subtitle = paste0("p = ", format(round(summary_df$p_cbi[summary_df$null == nm], 2), nsmall = 2)))
}
fig <- p_hist("uniform", "(a) Uniform null") +
  p_hist("accessibility", "(b) Accessibility-matched null")
save_fig(file.path(DIR_FIGS, "fig_null_model_panel.png"), fig,
         width = FIG_WIDTH_FULL, height = FIG_HEIGHT_PLOT)

# --------------------------------- Save -------------------------------------

write.csv(iters,      file.path(DIR_TABLES, "null_model_iterations.csv"), row.names = FALSE)
write.csv(summary_df, file.path(DIR_TABLES, "null_model_summary.csv"),    row.names = FALSE)
cat("13_null_model_test.R complete\n")









