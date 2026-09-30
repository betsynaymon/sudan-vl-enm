# ============================================================================
# helpers.R
# Shared functions, sourced after params.R. 
# ============================================================================

# Year-matched covariate values for a point set (longitude, latitude, year).
# Static covariates come from COV_FILES; dynamic ones from the annual raster
# for each point's year. Points must lie in the domain, and by default every
# value must be present: nothing is dropped silently.
extract_year_matched <- function(pts, vars, annual = COV_ANNUAL,
                                 require_complete = TRUE) {
  xy <- as.matrix(pts[, c("longitude", "latitude")])

  in_dom <- terra::extract(terra::rast(DOMAIN_FILE), xy)[, 1]
  if (!all(in_dom %in% 1)) stop(sum(!in_dom %in% 1), " point(s) outside the study domain")

  dyn  <- intersect(vars, names(annual))
  stat <- setdiff(vars, dyn)
  out  <- as.data.frame(matrix(NA_real_, nrow(pts), length(vars),
                               dimnames = list(NULL, vars)))

  if (length(stat) > 0) {
    r <- terra::rast(file.path(DIR_COVARIATES, COV_FILES[stat]))
    names(r) <- stat
    out[, stat] <- terra::extract(r, xy)
  }
  for (yr in sort(unique(pts$year))) {
    idx <- which(pts$year == yr)
    f <- file.path(DIR_COVARIATES, sub("{year}", yr, annual[dyn], fixed = TRUE))
    if (!all(file.exists(f))) stop("Missing annual raster(s) for ", yr)
    r <- terra::rast(f)
    names(r) <- dyn
    out[idx, dyn] <- terra::extract(r, xy[idx, , drop = FALSE])
  }

  if (require_complete && !all(complete.cases(out)))
    stop(sum(!complete.cases(out)), " point(s) have a missing covariate value")
  out
}

# Fold for each point, joined by ID from FOLD_TABLE_FILE (05).
attach_fold <- function(ids, pa_value, fold_table) {
  ft <- fold_table[fold_table$pa == pa_value, ]
  if (!setequal(ids, ft$id)) stop("Fold table and data disagree: rerun 05")
  ft$fold[match(ids, ft$id)]
}

# MaxEnt fit for a feature-class string (e.g. "LQH") and regularisation
# multiplier.
fit_maxnet <- function(p, data, fc, rm) {
  maxnet::maxnet(p = p, data = data,
                 f = maxnet::maxnet.formula(p = p, data = data, classes = tolower(fc)),
                 regmult = rm)
}

# Validation metrics for one fold. AUC counts ties as half; CBI uses the
# moving-window Boyce index; omission uses the OMISSION_Q training threshold.
eval_fold <- function(mod, test_occ, test_bg, train_occ) {
  p_occ <- as.numeric(predict(mod, test_occ,  type = "cloglog"))
  p_bg  <- as.numeric(predict(mod, test_bg,   type = "cloglog"))
  p_tr  <- as.numeric(predict(mod, train_occ, type = "cloglog"))

  auc <- mean(outer(p_occ, p_bg, ">") + 0.5 * outer(p_occ, p_bg, "=="))
  cbi <- tryCatch(
    ecospat::ecospat.boyce(fit = c(p_occ, p_bg), obs = p_occ, nclass = 0,
                           window.w = "default", res = BOYCE_RES,
                           PEplot = FALSE)$cor,
    error = function(e) NA_real_)
  or10 <- mean(p_occ < quantile(p_tr, OMISSION_Q))

  data.frame(auc = auc, cbi = cbi, or_10p = or10)
}