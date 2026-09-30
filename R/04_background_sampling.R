# ============================================================================
# 04_background_sampling.R
# Samples uniform-random background cells from the study domain, restricted
# to cells with complete values for the retained covariates, and assigns each
# a year drawn from the occurrence-year distribution for year-matched
# extraction downstream.
#
# Inputs:  OCC_FILE, DOMAIN_FILE, ADM0_FILE, SENS_MASK_FILE,
#          retained_vars.rds, COV_FILES
# Outputs: BG_FILE
#          outputs/figures/background_vs_occurrences.png
# ============================================================================

source(here::here("R", "params.R"))

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(ggplot2)
})

# ------------------------------ Load inputs ---------------------------------

occ <- read.csv(OCC_FILE)
cat("Occurrence records:", nrow(occ), "\n")

stopifnot(
  "No occurrence records found — run 02_spatial_thinning.R first" =
    nrow(occ) > 0
)

retained_vars <- readRDS(file.path(DIR_MODELS, "retained_vars.rds"))
covs   <- rast(file.path(DIR_COVARIATES, COV_FILES[retained_vars]))
domain <- rast(DOMAIN_FILE)

# Sampling surface: domain cells with every retained long-term covariate.
# Cells with a missing covariate (mostly coastline) cannot serve as background.
in_domain <- values(domain, mat = FALSE) == 1
avail     <- in_domain & complete.cases(values(covs))
cat("Domain cells:", sum(in_domain), "| available for background:", sum(avail), "\n")

# Presences must lie where background can be drawn; otherwise the model
# contrasts them with an area they are not part of.
occ_cells <- cellFromXY(domain, as.matrix(occ[, c("longitude", "latitude")]))
stopifnot(
  "A presence falls outside the covariate grid"   = !anyNA(occ_cells),
  "A presence falls outside the sampling surface" = all(avail[occ_cells])
)

# -------------------------- Sample background -------------------------------
set.seed(SEED)
bg_cells <- sample(which(avail), N_BACKGROUND)   # without replacement: one point per cell
bg <- as.data.frame(xyFromCell(domain, bg_cells))
names(bg) <- c("longitude", "latitude")
bg <- data.frame(bg_id = seq_len(nrow(bg)), bg)

stopifnot("Background count differs from N_BACKGROUND" = nrow(bg) == N_BACKGROUND)

# Background where the presences are 
sens <- values(rast(SENS_MASK_FILE), mat = FALSE)
cat("Background in the >= 150 mm part of the domain:",
    sum(sens[bg_cells] %in% 1), "of", N_BACKGROUND, "\n")

# Assign years weighted by occurrence-year distribution
year_weights <- occ %>% count(year, name = "weight")

bg$year <- sample(year_weights$year, size = nrow(bg), replace = TRUE,
                  prob = year_weights$weight)

cat("\nYear distribution comparison:\n")
occ_pct <- occ %>% count(year) %>% mutate(occ_pct = round(100 * n / sum(n), 1))
bg_pct  <- bg  %>% count(year) %>% mutate(bg_pct  = round(100 * n / sum(n), 1))
left_join(occ_pct, bg_pct, by = "year", suffix = c("_occ", "_bg")) %>%
  select(year, n_occ, occ_pct, n_bg, bg_pct) %>%
  as.data.frame() %>% print()

# ----------------------------- Save and map ---------------------------------
write.csv(bg, BG_FILE, row.names = FALSE)
cat("Saved", nrow(bg), "background points\n")

p <- ggplot() +
  geom_sf(data = sf::st_as_sf(vect(ADM0_FILE)), fill = NA, colour = "grey30") +
  geom_point(data = bg, aes(longitude, latitude),
             colour = "grey70", size = 0.3, alpha = 0.3) +
  geom_point(data = occ, aes(longitude, latitude),
             colour = "firebrick", size = 1.5) +
  coord_sf() +
  labs(title = paste0("Background (n=", nrow(bg),
                      ") and occurrences (n=", nrow(occ), ")"),
       x = "Longitude", y = "Latitude") +
  theme_minimal()

ggsave(file.path(DIR_FIGS, "background_vs_occurrences.png"), p,
       width = 8, height = 7, dpi = 300)
cat("Saved background_vs_occurrences.png\n")

cat("04_background_sampling.R complete\n")