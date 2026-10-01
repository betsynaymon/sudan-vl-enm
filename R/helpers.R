# ============================================================================
# helpers.R
# Shared functions, sourced after params.R. 
# ============================================================================

# AUC for presences vs background, ties counted as half.
auc_ties <- function(p_occ, p_bg) {
  mean(outer(p_occ, p_bg, ">") + 0.5 * outer(p_occ, p_bg, "=="))
}

# Response curves: each covariate varied across the range of `env`, the
# others held at `ref`. Binary covariates are evaluated at 0 and 1.
response_curves <- function(mod, ref, env, vars, n = 200) {
  do.call(rbind, lapply(vars, function(v) {
    x <- if (all(env[[v]] %in% c(0, 1))) c(0, 1) else
      seq(min(env[[v]]), max(env[[v]]), length.out = n)
    nd <- as.data.frame(matrix(rep(ref, each = length(x)), ncol = length(ref),
                               dimnames = list(NULL, names(ref))))
    nd[[v]] <- x
    data.frame(variable = v, value = x,
               suit = as.numeric(predict(mod, nd, type = "cloglog", clamp = TRUE)))
  }))
}

# Suitability over a grid of two covariates, the others held at `ref`.
pair_surface <- function(mod, ref, env, x, y, n = 100) {
  g <- expand.grid(seq(min(env[[x]]), max(env[[x]]), length.out = n),
                   seq(min(env[[y]]), max(env[[y]]), length.out = n))
  names(g) <- c(x, y)
  nd <- as.data.frame(matrix(rep(ref, each = nrow(g)), ncol = length(ref),
                             dimnames = list(NULL, names(ref))))
  nd[[x]] <- g[[x]]; nd[[y]] <- g[[y]]
  cbind(g, suit = as.numeric(predict(mod, nd, type = "cloglog", clamp = TRUE)))
}


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

  auc <- auc_ties(p_occ, p_bg)
  cbi <- tryCatch(
    ecospat::ecospat.boyce(fit = c(p_occ, p_bg), obs = p_occ, nclass = 0,
                           window.w = "default", res = BOYCE_RES,
                           PEplot = FALSE)$cor,
    error = function(e) NA_real_)
  or10 <- mean(p_occ < quantile(p_tr, OMISSION_Q))

  data.frame(auc = auc, cbi = cbi, or_10p = or10)
}

# p10 (OMISSION_Q quantile of presence predictions) and maxSSS (maximum
# sensitivity + specificity against the background).
thresholds_from <- function(pred_occ, pred_bg) {
  cand <- sort(unique(c(pred_occ, pred_bg)))
  sss  <- sapply(cand, function(t) mean(pred_occ >= t) + mean(pred_bg < t))
  c(p10 = unname(quantile(pred_occ, OMISSION_Q)), maxsss = cand[which.max(sss)])
}

# Thousands separator for printed counts.
fmt <- function(x) format(round(x), big.mark = ",")

# Long-term covariates for `vars`, masked to the domain: the prediction stack
# (06b, 07 and every refit).
domain_covs <- function(vars) {
  r <- terra::rast(file.path(DIR_COVARIATES, COV_FILES[vars]))
  names(r) <- vars
  terra::mask(r, terra::rast(DOMAIN_FILE), maskvalues = 0)
}

# Spatial CV for one configuration, as in 06: each fold's model is fitted on
# the other folds and scored on its own. `fold` is each row's own fold. A
# failed fit stops rather than being skipped.
cv_maxnet <- function(pa, env, fold, fc, rm) {
  do.call(rbind, lapply(seq_len(K_FOLDS), function(k) {
    tr <- fold != k
    m  <- fit_maxnet(pa[tr], env[tr, ], fc, rm)
    cbind(fold = k, n_test_pres = sum(!tr & pa == 1),
          eval_fold(m, test_occ  = env[!tr & pa == 1, ],
                       test_bg   = env[!tr & pa == 0, ],
                       train_occ = env[tr & pa == 1, ]))
  }))
}

# Fit, cross-validate and threshold one configuration on a presence and
# background set: the refit used by 10, 12 and 15.
refit_maxnet <- function(occ_env, bg_env, occ_fold, bg_fold, fc, rm) {
  pa  <- c(rep(1, nrow(occ_env)), rep(0, nrow(bg_env)))
  env <- rbind(occ_env, bg_env)
  mod <- fit_maxnet(pa, env, fc, rm)
  list(mod = mod,
       cv  = cv_maxnet(pa, env, c(occ_fold, bg_fold), fc, rm),
       thr = thresholds_from(as.numeric(predict(mod, occ_env, type = "cloglog")),
                             as.numeric(predict(mod, bg_env,  type = "cloglog"))))
}

# Population at risk: risk-weighted (population x suitability) and binary at
# each threshold, over cells with both a prediction and a population value.
arp_estimates <- function(suit_r, pop_r, thr) {
  v <- data.frame(s = terra::values(suit_r, mat = FALSE),
                  p = terra::values(pop_r,  mat = FALSE))
  v <- v[complete.cases(v), ]
  c(risk_weighted = sum(v$p * v$s),
    p10    = sum(v$p[v$s >= thr[["p10"]]]),
    maxsss = sum(v$p[v$s >= thr[["maxsss"]]]))
}

# One state per domain cell. touches = TRUE keeps the domain's boundary cells;
# cells on internal borders go to one state, never two.
state_zones <- function(template) {
  z <- terra::rasterize(terra::vect(ADM1_FILE), template, field = "NAME_1", touches = TRUE)
  terra::mask(z, terra::rast(DOMAIN_FILE), maskvalues = 0)
}

# Weiss et al. (2018) travel time to the nearest city (minutes) on the grid of
# `template`, masked to the domain: the surface 11, 12, 13 and 22 use.
# method = "near" is for diagnosing gaps only.
travel_time <- function(template, method = "bilinear") {
  tt <- terra::resample(terra::rast(TT_FILE), template, method = method)
  names(tt) <- "travel_time"
  terra::mask(tt, terra::rast(DOMAIN_FILE), maskvalues = 0)
}

# Fold for new points (longitude, latitude) from the spatial blocks of 05: the
# block containing each point, or the nearest block if none contains it.
fold_from_blocks <- function(pts) {
  blocks <- readRDS(FOLDS_FILE)$blocks
  p <- sf::st_as_sf(pts, coords = c("longitude", "latitude"), crs = 4326)
  p <- sf::st_transform(p, sf::st_crs(blocks))
  blocks$folds[sf::st_nearest_feature(p, blocks)]
}