# ============================================================================
# helpers.R
# Shared functions, sourced after params.R. 
# ============================================================================

# AUC for presences vs background, ties counted as half.
auc_ties <- function(p_occ, p_bg) {
  mean(outer(p_occ, p_bg, ">") + 0.5 * outer(p_occ, p_bg, "=="))
}

# Response curves: each covariate varied across the range of `env`, the
# others held at `ref`. Binary covariates are evaluated at 0 and 1. `pred`
# returns suitability for a model and data frame: MaxEnt cloglog (clamped) by
# default, rf_prob for the random forest.
response_curves <- function(mod, ref, env, vars, n = 200,
                            pred = function(m, d) predict(m, d, type = "cloglog", clamp = TRUE)) {
  do.call(rbind, lapply(vars, function(v) {
    x <- if (all(env[[v]] %in% c(0, 1))) c(0, 1) else
      seq(min(env[[v]]), max(env[[v]]), length.out = n)
    nd <- as.data.frame(matrix(rep(ref, each = length(x)), ncol = length(ref),
                               dimnames = list(NULL, names(ref))))
    nd[[v]] <- x
    data.frame(variable = v, value = x, suit = as.numeric(pred(mod, nd)))
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

# Validation metrics from predictions, for any algorithm: test presences,
# test background, and training presences (for the OMISSION_Q threshold).
# AUC counts ties as half; CBI is the moving-window Boyce index.
score_preds <- function(p_occ, p_bg, p_tr) {
  cbi <- tryCatch(
    ecospat::ecospat.boyce(fit = c(p_occ, p_bg), obs = p_occ, nclass = 0,
                           window.w = "default", res = BOYCE_RES,
                           PEplot = FALSE)$cor,
    error = function(e) NA_real_)
  data.frame(auc = auc_ties(p_occ, p_bg), cbi = cbi,
             or_10p = mean(p_occ < quantile(p_tr, OMISSION_Q)))
}

# MaxEnt suitability (cloglog) as a vector.
maxnet_prob <- function(mod, data) as.numeric(predict(mod, data, type = "cloglog"))

# Validation metrics for one MaxEnt fold (06, 13).
eval_fold <- function(mod, test_occ, test_bg, train_occ)
  score_preds(maxnet_prob(mod, test_occ), maxnet_prob(mod, test_bg),
              maxnet_prob(mod, train_occ))


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
# (06b, 07 and every refit). `files` defaults to COV_FILES; 15 adds the
# seasonal composites.
domain_covs <- function(vars, files = COV_FILES) {
  r <- terra::rast(file.path(DIR_COVARIATES, files[vars]))
  names(r) <- vars
  terra::mask(r, terra::rast(DOMAIN_FILE), maskvalues = 0)
}
# Spatial CV for any algorithm, as in 06: each fold's model is fitted on the
# other folds and scored on its own. fit(pa, env) returns a model and
# pred(model, env) its suitability; `fold` is each row's own fold. With
# `within` (one logical per row), metrics are also computed on the test rows
# inside it (suffix _wet); the omission threshold still uses all training
# presences, as in 06. A failed fit stops rather than being skipped.
cv_fit <- function(pa, env, fold, fit, pred, within = NULL) {
  do.call(rbind, lapply(seq_len(K_FOLDS), function(k) {
    tr   <- fold != k
    m    <- fit(pa[tr], env[tr, ])
    p_tr <- pred(m, env[tr & pa == 1, ])
    out  <- cbind(fold = k, n_test_pres = sum(!tr & pa == 1),
                  score_preds(pred(m, env[!tr & pa == 1, ]),
                              pred(m, env[!tr & pa == 0, ]), p_tr))
    if (!is.null(within)) {
      w <- score_preds(pred(m, env[!tr & pa == 1 & within, ]),
                       pred(m, env[!tr & pa == 0 & within, ]), p_tr)
      names(w) <- paste0(names(w), "_wet")
      out <- cbind(out, w)
    }
    out
  }))
}

# Spatial CV for one MaxEnt configuration.
cv_maxnet <- function(pa, env, fold, fc, rm, within = NULL)
  cv_fit(pa, env, fold, function(p, d) fit_maxnet(p, d, fc, rm), maxnet_prob, within)

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

# Weighted median of travel time `tt` (minutes) under weights (1 + tt)^-b.
tt_wmedian <- function(tt, b) {
  s <- sort(tt[!is.na(tt)]); w <- (1 + s)^(-b)
  s[which(cumsum(w) >= sum(w) / 2)[1]]
}

# Exponent b at which the weighted median of `tt` (candidate cells) equals
# `target`: the accessibility weighting that makes a weighted sample as close
# to cities as the presences (12, 13).
tt_exponent <- function(tt, target, grid = seq(0, 3, by = 0.01)) {
  s   <- sort(tt[!is.na(tt)])
  med <- sapply(grid, function(b) { w <- (1 + s)^(-b); s[which(cumsum(w) >= sum(w) / 2)[1]] })
  b   <- grid[which.min(abs(med - target))]
  if (b <= min(grid) || b >= max(grid)) stop("Calibrated exponent at the edge of the search range")
  b
}

# Weighted sample of cells (with replacement, so exactly proportional to `w`;
# NULL = uniform) as points with a year drawn from `years` (a table of
# counts), year-matched values, and folds from 05's blocks. Points missing an
# annual value are dropped and counted. Used for 13's null pools and 12's
# corrected backgrounds.
draw_points <- function(template, cells, n, w, years, vars, seed) {
  set.seed(seed)
  idx <- sample(cells, n, replace = TRUE, prob = w)
  pts <- as.data.frame(terra::xyFromCell(template, idx))
  names(pts) <- c("longitude", "latitude")
  pts$year <- as.integer(sample(names(years), n, replace = TRUE, prob = years))
  env  <- extract_year_matched(pts, vars, require_complete = FALSE)
  keep <- complete.cases(env)
  list(pts = pts[keep, ], env = env[keep, vars], fold = fold_from_blocks(pts[keep, ]),
       cell = idx[keep], dropped = sum(!keep))
}

# In-sample permutation importance: AUC drop when one covariate is shuffled
# across presences `o` and background `b` (07, 24). N_PERM permutations per
# variable; set the seed before calling.
perm_importance <- function(mod, o, b, vars, n = N_PERM) {
  base <- auc_ties(as.numeric(predict(mod, o, type = "cloglog")),
                   as.numeric(predict(mod, b, type = "cloglog")))
  t(sapply(vars, function(v) {
    drops <- replicate(n, {
      shuf <- sample(c(o[[v]], b[[v]]))
      o2 <- o; b2 <- b
      o2[[v]] <- shuf[seq_len(nrow(o))]
      b2[[v]] <- shuf[-seq_len(nrow(o))]
      base - auc_ties(as.numeric(predict(mod, o2, type = "cloglog")),
                      as.numeric(predict(mod, b2, type = "cloglog")))
    })
    c(mean = mean(drops), sd = sd(drops))
  }))
}

# Positions read off one response curve (columns variable, value, suit;
# value ascending): the peak, and the lowest (rise_XX) and highest (fall_XX)
# values where suitability is at least XX% of the peak. Other covariates are
# held at a reference, so positions are more robust than heights (07, 09, 24).
curve_features <- function(d) {
  mx <- max(d$suit)
  x_at <- function(frac, side) {
    above <- d$value[d$suit >= frac * mx]
    if (side == "rise") min(above) else max(above)
  }
  tibble::tibble(
    variable  = unique(d$variable),
    peak_x    = d$value[which.max(d$suit)],
    peak_suit = mx,
    min_suit  = min(d$suit),
    rise_10 = x_at(0.10, "rise"), rise_50 = x_at(0.50, "rise"),
    rise_90 = x_at(0.90, "rise"),
    fall_90 = x_at(0.90, "fall"), fall_50 = x_at(0.50, "fall"),
    fall_10 = x_at(0.10, "fall")
  )
}

# Random forest for presence-background data, down-sampled (Valavi et al.
# 2021): a probability forest whose every tree draws, with replacement, as
# many presences and as many background points as there are presences.
# ranger seeds each tree from SEED, so fits are reproducible.
fit_rf <- function(pa, env, mtry, keep.inbag = FALSE) {
  n_pres <- sum(pa == 1)
  ranger::ranger(x = env, y = factor(pa, levels = c(0, 1)), probability = TRUE,
                 num.trees = RF_NTREES, mtry = mtry, replace = TRUE,
                 sample.fraction = rep(n_pres / length(pa), 2),
                 keep.inbag = keep.inbag, seed = SEED, verbose = FALSE)
}

# Random-forest probability of presence; also the `fun` for terra::predict.
rf_prob <- function(mod, data, ...)
  predict(mod, data = data, verbose = FALSE)$predictions[, "1"]

# ---- Comparisons across algorithms (06, 07, 09, 17) ----

# Each presence predicted by the model fitted without its fold, and whether it
# falls below that model's OMISSION_Q threshold (from its own training
# presences): does the model predict a record without its neighbours? One row
# per presence, in presence order (06, 17).
heldout_preds <- function(pa, env, fold, fit, pred) {
  out <- data.frame(fold = fold[pa == 1], pred = NA_real_, thr = NA_real_)
  for (k in seq_len(K_FOLDS)) {
    tr <- fold != k
    m  <- fit(pa[tr], env[tr, ])
    out$pred[out$fold == k] <- pred(m, env[!tr & pa == 1, ])
    out$thr[out$fold == k]  <- quantile(pred(m, env[tr & pa == 1, ]), OMISSION_Q, names = FALSE)
  }
  out$below_threshold <- out$pred < out$thr
  out
}

# curve_features() for every continuous covariate in a set of response
# curves; vertisols is binary (07, 09, 17).
curve_table <- function(curves) {
  cont <- curves[curves$variable != "vertisols", ]
  dplyr::bind_rows(lapply(split(cont, cont$variable), curve_features))
}

# Suitability at x as a percentage of the peak of one covariate's response
# curve (09, 17: rainfall at the wettest presence).
pct_of_peak <- function(curves, var, x) {
  r <- curves[curves$variable == var, ]
  100 * stats::approx(r$value, r$suit, xout = x)$y / max(r$suit)
}

# One algorithm's row in the comparison (09, 17): cross-validated fit, its own
# thresholds, population at risk, and rank agreement with MaxEnt's surface
# (NA for MaxEnt itself).
comparator_row <- function(model, config, cv, thr, arp,
                           rho = c(domain = NA_real_, belt = NA_real_))
  data.frame(model = model, config = config,
             cv_cbi = mean(cv$cbi), cv_cbi_sd = sd(cv$cbi), cv_cbi_wet = mean(cv$cbi_wet),
             cv_auc = mean(cv$auc), cv_or_10p = mean(cv$or_10p),
             p10 = thr[["p10"]], maxsss = thr[["maxsss"]],
             arp_weighted = arp[["risk_weighted"]], arp_p10 = arp[["p10"]],
             arp_maxsss = arp[["maxsss"]],
             rho_domain = rho[["domain"]], rho_belt = rho[["belt"]])

# Each algorithm against MaxEnt (09, 17): CV CBI within one SE of MaxEnt's
# (`se`), and the risk-weighted change against the plateau spread (%).
vs_primary <- function(df, se, plateau)
  dplyr::mutate(df,
    cbi_within_se  = cv_cbi >= cv_cbi[model == "maxent"] - se,
    rw_change_pct  = 100 * (arp_weighted / arp_weighted[model == "maxent"] - 1),
    beyond_plateau = abs(rw_change_pct) > plateau)

# Gradient boosted trees for presence-background data (17): presences weight
# 1, background down-weighted to the same total (Valavi et al. 2022). Unlike
# ranger's case weights these scale each row's loss, not its chance of being
# drawn: every tree draws GBT_BAG of the rows uniformly. gbm rescales weights
# to sum to the number of rows and counts rows for the minimum leaf size, so
# only the ratio matters. Seeded, so the first n trees of a longer fit equal
# an n-tree fit (17 checks).
fit_gbt <- function(pa, env, lr, depth, n_trees) {
  set.seed(SEED)
  gbm::gbm.fit(x = env, y = pa, distribution = "bernoulli",
               w = ifelse(pa == 1, 1, sum(pa == 1) / sum(pa == 0)),
               n.trees = n_trees, interaction.depth = depth, shrinkage = lr,
               bag.fraction = GBT_BAG, n.minobsinnode = GBT_MIN_NODE,
               keep.data = FALSE, verbose = FALSE)
}

# Boosted-tree probability of presence at the model's own number of trees;
# also the `fun` for terra::predict.
gbt_prob <- function(mod, data, ...)
  predict(mod, newdata = data, n.trees = mod$n.trees, type = "response")