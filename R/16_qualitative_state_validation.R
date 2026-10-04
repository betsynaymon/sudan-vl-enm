# ============================================================================
# 16_qualitative_state_validation.R
# Does predicted suitability follow the literature's ranking of states by VL
# endemicity (STATE_STATUS_FILE)? Core-endemic status rests largely on the
# literature that supplied the presences, so states with training records
# rank high partly by construction; the informative comparison is among
# states with no training presence. Reported for area-mean suitability
# (habitat) and population-weighted suitability (what the estimate uses).
#
# Inputs:  SUIT_FILE, POP_ALIGNED_FILE, TRAIN_FILE, DOMAIN_FILE, ADM1_FILE,
#          STATE_STATUS_FILE
# Outputs: outputs/tables/state_validation.csv
# ============================================================================

source(here::here("R", "params.R"))
source(here::here("R", "helpers.R"))

suppressPackageStartupMessages({
  library(terra); library(dplyr)
})

status_levels <- c("No documented", "Reported", "Core endemic")

# ------------------------------ Load inputs ---------------------------------

suit_r <- rast(SUIT_FILE)
pop    <- rast(POP_ALIGNED_FILE)
train  <- readRDS(TRAIN_FILE)
lit    <- read.csv(STATE_STATUS_FILE)
zones  <- state_zones(suit_r)
cat("States classified:", nrow(lit), "| without a source:", sum(lit$source %in% c("", NA)), "\n")

# ---------------------------- State summaries -------------------------------
# Same state assignment as 08. Population-weighted suitability = 08's
# risk-weighted estimate / state population.

s_mean <- zonal(suit_r, zones, fun = "mean",   na.rm = TRUE)
s_med  <- zonal(suit_r, zones, fun = "median", na.rm = TRUE)
s_pw   <- zonal(c(pop * suit_r, pop), zones, fun = "sum", na.rm = TRUE)
names(s_mean) <- c("state", "mean_suit")
names(s_med)  <- c("state", "median_suit")
names(s_pw)   <- c("state", "rw", "pop")
occ_state <- as.character(terra::extract(zones,
               as.matrix(train$occ_clean[, c("longitude", "latitude")]))[, 1])

stopifnot(
  "State names differ between STATE_STATUS_FILE and ADM1_FILE" = setequal(lit$state, s_mean$state),
  "Unknown status in STATE_STATUS_FILE" = all(lit$status %in% status_levels)
)

val <- lit |>
  left_join(s_mean, by = "state") |> left_join(s_med, by = "state") |>
  left_join(s_pw, by = "state") |>
  mutate(pop_weighted_suit = rw / pop,
         n_presences = as.integer(table(factor(occ_state, levels = state))),
         status = factor(status, levels = status_levels)) |>
  arrange(desc(status), desc(pop_weighted_suit))
stopifnot("Not every presence was assigned to a state" =
            sum(val$n_presences) == nrow(train$occ_clean))

# ------------------------------- Agreement ----------------------------------

rho <- function(d, col) {
  ct <- suppressWarnings(cor.test(as.integer(d$status), d[[col]],
                                  method = "spearman", exact = FALSE))
  c(states = nrow(d), rho = unname(ct$estimate), p = ct$p.value)
}
no_pres <- filter(val, n_presences == 0)
agree <- rbind(
  `All states, area mean`            = rho(val, "mean_suit"),
  `All states, population-weighted`  = rho(val, "pop_weighted_suit"),
  `No presences, area mean`          = rho(no_pres, "mean_suit"),
  `No presences, population-weighted` = rho(no_pres, "pop_weighted_suit")
)
cat("\nSpearman rho, endemic status vs suitability:\n"); print(round(agree, 3))

cat("\nBy state (status, then population-weighted suitability):\n")
val |> select(state, status, n_presences, mean_suit, median_suit, pop_weighted_suit) |>
  mutate(across(c(mean_suit, median_suit, pop_weighted_suit), ~ round(., 3))) |>
  print(row.names = FALSE, right = FALSE)

# --------------------------------- Save -------------------------------------

write.csv(select(val, state, status, evidence, source, n_presences,
                 mean_suit, median_suit, pop_weighted_suit),
          file.path(DIR_TABLES, "state_validation.csv"), row.names = FALSE)
cat("16_qualitative_state_validation.R complete\n")