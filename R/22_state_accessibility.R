# ============================================================================
# 22_state_accessibility.R
# Share of each state's population living more than TT_REMOTE_MIN minutes from
# the nearest city (Weiss et al. 2018), for the state comparison table. Same
# grid and population as 08 (travel time resampled to the covariate grid;
# WorldPop aligned in 08). People in cells without travel time (11: the Red
# Sea coast, missing in the source grid) are reported separately and excluded
# from the share.
#
# Inputs:  TT_FILE, POP_ALIGNED_FILE, SUIT_FILE (grid), DOMAIN_FILE, ADM1_FILE,
#          outputs/tables/arp_by_state.csv (08)
# Outputs: outputs/tables/state_accessibility.csv
# ============================================================================

source(here::here("R", "params.R"))
source(here::here("R", "helpers.R"))

suppressPackageStartupMessages({
  library(terra); library(dplyr)
})

# ------------------------------ Load inputs ---------------------------------

suit_r    <- rast(SUIT_FILE)
pop       <- rast(POP_ALIGNED_FILE)
tt        <- travel_time(suit_r)
zones     <- state_zones(suit_r)
states_08 <- read.csv(file.path(DIR_TABLES, "arp_by_state.csv"))
stopifnot("Population grid differs from the surface" =
            compareGeom(pop, suit_r, stopOnError = FALSE))

# ------------------------------- By state -----------------------------------

layers <- c(pop, pop * (tt > TT_REMOTE_MIN), pop * is.na(tt))
names(layers) <- c("total_pop", "pop_beyond", "pop_no_tt")
st <- zonal(layers, zones, fun = "sum", na.rm = TRUE)
names(st)[1] <- "state"

check <- left_join(st, select(states_08, state, total_pop_08 = total_pop), by = "state")
stopifnot("State population differs from 08" =
            isTRUE(all.equal(check$total_pop, check$total_pop_08)))

st <- st |>
  mutate(pop_with_tt = total_pop - pop_no_tt,
         pct_beyond  = round(100 * pop_beyond / pop_with_tt, 1)) |>
  arrange(desc(pct_beyond))
nat <- colSums(st[c("total_pop", "pop_beyond", "pop_no_tt", "pop_with_tt")])

cat("Share of each state's population more than", TT_REMOTE_MIN, "min from a city:\n")
st |> select(state, total_pop, pop_beyond, pct_beyond, pop_no_tt) |>
  mutate(across(c(total_pop, pop_beyond, pop_no_tt), fmt)) |>
  print(right = FALSE, row.names = FALSE)
cat("\nNational:", fmt(nat[["pop_beyond"]]), "of", fmt(nat[["pop_with_tt"]]),
    "people with a travel time (", round(100 * nat[["pop_beyond"]] / nat[["pop_with_tt"]], 1),
    "% ) | without travel time:", fmt(nat[["pop_no_tt"]]), "\n")

# --------------------------------- Save -------------------------------------

write.csv(st, file.path(DIR_TABLES, "state_accessibility.csv"), row.names = FALSE)
cat("22_state_accessibility.R complete\n")