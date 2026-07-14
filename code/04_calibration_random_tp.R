#===============================================================================
# 03_calibration_indonesia_random_tp.R
#-------------------------------------------------------------------------------
# RANDOM-SEARCH TP-PERTURBATION CALIBRATION FOR THE INDONESIA NCD / CVD MARKOV
# MODEL
#
# RELATIONSHIP TO 03_calibration_indonesia_transparent.R
# ------------------------------------------------------
# This file is a NEAR-CLONE of 03_calibration_indonesia_transparent.R. The ONLY
# algorithmic change is the per-combo SEARCH STRATEGY:
#
#   * transparent.R : a (1+1) stochastic HILL-CLIMB. Each candidate is a Gaussian
#                     perturbation of the CURRENT BEST multiplier vector, folded
#                     back into the search bounds, with periodic uniform restarts
#                     and a "no-improvement" early stop.
#   * THIS FILE      : a PURE RANDOM SEARCH (Monte Carlo). Each candidate is an
#                     INDEPENDENT, i.i.d. uniform draw of the multiplier vector,
#                     drawn WITHOUT any reference to the current best. The best
#                     candidate is merely RETAINED (argmin over evaluations); it
#                     never seeds the next draw. There is no Gaussian step, no
#                     boundary reflection, and no restart schedule -- those are
#                     hill-climb constructs with no meaning for random search.
#
# Everything else -- data loading, GBD target prep, the well-sick-dead Markov
# projection (project_combo), the objective (combo_error), the probability
# constraints (enforce_tp_constraints), the multiplier algebra (build_mtab /
# apply_multipliers), the run_combo wrapper, and the output schema -- is REUSED
# UNCHANGED from the transparent script. This keeps the diff between the two
# files confined to calibrate_one_combo* so the two algorithms can be compared
# side by side. The cancer script (calibration (1).R) is used ONLY as a design
# reference for the random-perturbation / keep-best pattern; none of its disease
# stages (lcl/rgn/dst/prc) are imported -- the CVD model stays well-sick-dead.
#
# CONTEXT (carried over from transparent.R)
# -----------------------------------------
# The production calibration was split across two opaque grid searches:
#   * 031_calibration_indonesia.R : a symmetric +/-5% grid over IR and CF
#                                   (for(i in -5:5) for(j in -5:5) -> 1 + i/100).
#   * 032_adjustments_indonesia.R : a SECOND, asymmetric grid that shrinks IR by
#                                   0-20% and CF by 0-50% and then renormalises
#                                   IR/CF against background mortality.
# Both transparent.R and this file REPLACE that 031 + 032 split with ONE search
# over a single documented SYMMETRIC multiplicative range around 1
# (SEARCH_HALFWIDTH), differing only in HOW that range is explored.
#
# IT DOES NOT MODIFY 031, 032 OR ANY OTHER EXISTING SCRIPT. It only reads the
# same inputs (tps_inpt_part*.rds) and writes the same outputs
# (adjusted_searo_part{1..10}.rds in wd_data) so it is a DROP-IN replacement
# for sourcing 03_calibration.R (i.e. 031 + 032).
#
# HOW TO USE
# ----------
# Source this AFTER 02_load_inputs_indonesia.R (so tps_inpt_part*.rds and the
# globals wd / wd_raw / wd_data / wd_temp / dx_include / cause_map / locs exist)
# INSTEAD OF sourcing 03_calibration.R. Then continue with 04/05/06/07/08.
#
# IMPORTANT -- AVOID DOUBLE-CALIBRATION:
#   05_build_baseline_indonesia.R re-applies adjustments2023_age.csv when
#   run_adjustment_model == TRUE. That CSV is produced by 032, which this script
#   does NOT run. To avoid multiplying IR/CF a second time with stale 032
#   factors, set  run_adjustment_model <- FALSE  before running 05. This script
#   already bakes the calibrated multipliers into the saved IR/CF, so no further
#   adjustment is needed downstream.
#
# OBJECTIVE FUNCTION
# ------------------
# For each location-sex-cause combo we project the cohort 2009-2019 and compare
# to GBD "Number" of Deaths (fatal) and Prevalence (non-fatal), aggregated to
# 5-year age groups. The search minimises a weighted RELATIVE squared error
# (so age groups of very different magnitudes contribute comparably, matching
# calc_combo_error_fast in the cancer template):
#
#   error = sum_{year, age.group} [
#               W_DEATHS * ((Deaths_model - Deaths_gbd)/(Deaths_gbd + EPS_REL))^2
#             + W_PREV   * ((Prev_model   - Prev_gbd  )/(Prev_gbd   + EPS_REL))^2 ]
#
# Fatal estimates are weighted 2x non-fatal (W_DEATHS = 2, W_PREV = 1), exactly
# as 031's "error = 2*RMSE_deaths + RMSE_prev". Absolute RMSE_deaths / RMSE_prev
# are ALSO reported in the diagnostics table for comparison with the old method.
#
# SEARCH SPACE & STRATEGY (THE ONLY DIFFERENCE FROM transparent.R)
# ----------------------------------------------------------------
# The free parameters are multiplicative adjustment factors on the baseline
# transition probabilities IR and CF, sampled symmetrically in
#   [1 - SEARCH_HALFWIDTH, 1 + SEARCH_HALFWIDTH].
# Two granularities (set GRANULARITY below):
#   * "combo"      (PRIMARY)     : ONE IR multiplier + ONE CF multiplier shared
#                                  by all ages in the combo (2 free params).
#   * "age_group"  (SENSITIVITY) : a separate IR and CF multiplier per 5-year
#                                  age group (2 * n_age_groups free params).
#
# STRATEGY = PURE RANDOM SEARCH. For each combo we evaluate candidate 0 = the
# baseline (all multipliers = 1), then draw N_ITER candidates, each an
# INDEPENDENT i.i.d. uniform vector:
#       m^{IR}_{g,i}, m^{CF}_{g,i}  ~  U(1 - SEARCH_HALFWIDTH, 1 + SEARCH_HALFWIDTH).
# Candidate i does NOT depend on the current best (contrast transparent.R, where
# the next candidate is a Gaussian step from the best). We keep the argmin over
# all evaluations:  theta* = argmin_i E(theta_i). Because candidate 0 is the
# baseline, the calibrated fit can never be WORSE than the uncalibrated baseline.
#
# PROBABILITY / ROW CONSTRAINTS  (enforced after EVERY adjustment)
# ----------------------------------------------------------------
# IR is used in the recursion as an annual transition PROBABILITY (well -> sick),
# and CF as a probability (sick -> dead). Therefore after each adjustment:
#   * NA -> 0; then 0 <= IR <= 1 and 0 <= CF <= 1.
#   * IR + BG.mx <= 1 and CF + BG.mx <= 1 (a sick/well row's competing risks
#     cannot exceed 1).
# We PREFER to preserve the externally-supplied background mortality (BG.mx) and
# instead cap the DISEASE transition probability:
#       IR <= 1 - BG.mx - TP_EPS ;  CF <= 1 - BG.mx - TP_EPS.
# Only as a FALLBACK -- when BG.mx alone already leaves no room
# (BG.mx >= 1 - TP_EPS) -- do we proportionally renormalise IR/CF AND BG.mx
# (the share-based shrink from 032_adjustments_indonesia.R lines 46-55). Any row
# where BG.mx was modified is flagged and the count is reported, so the rare
# departure from "preserve BG.mx" is fully auditable.
#
# OUTPUT CONTRACT
# ---------------
#   * adjusted_searo_part{1..10}.rds in wd_data -- identical schema/format/
#     chunking to 031's output: original tps_inpt rows with IR/CF (and, in the
#     rare fallback, BG.mx) updated in place. This is what 05_build_baseline
#     picks up via its "adjusted" file pattern.
#   * calibration_factors_random_tp.csv -- the calibrated per-combo (and per
#     age-group) IR/CF multipliers, analogous to adjustments2023_age.csv.
#   * calibration_diagnostics_random_tp.csv -- baseline-vs-calibrated absolute
#     RMSE_deaths / RMSE_prev, weighted error, and % improvement per
#     location-sex-cause(-age.group).
#===============================================================================

# library(data.table)
# library(foreach)
# library(doParallel)
# library(parallel)

#===============================================================================
# 0. TUNABLE PARAMETERS  (all calibration "magic numbers" live here, commented)
#===============================================================================

## --- search SPACE -----------------------------------------------------------
# Symmetric half-width of the multiplicative search range around 1. 0.5 means
# IR/CF may be scaled anywhere in [0.5, 1.5]. This single symmetric range
# REPLACES 031's +/-5% grid and 032's asymmetric 0-20% IR / 0-50% CF grid; it is
# wide enough to subsume both (their combined reach was ~ -55% .. +5%) while
# being explicit and symmetric. Widen only if combos hit the bound.
SEARCH_HALFWIDTH <- 0.50

# Granularity of the multipliers. "combo" = PRIMARY (one IR + one CF per
# location-sex-cause). "age_group" = SENSITIVITY (per 5-year age group).
GRANULARITY <- "age_group"          # "combo" (primary) | "age_group" (sensitivity)

## --- search ALGORITHM (PURE RANDOM SEARCH / Monte Carlo) --------------------
# Each candidate is an INDEPENDENT i.i.d. uniform draw of the multiplier vector
# (no Gaussian step, no reflection, no restart schedule -- those are hill-climb
# constructs and do not apply here). We keep the argmin over evaluations.
N_ITER       <- 400         # number of random candidates drawn per combo
                            # (the evaluation budget; baseline is candidate 0)
CONVERGE_TOL <- 1e-4        # stop a combo early if best weighted error < this
                            # (the only early stop; random search has no
                            #  meaningful "no-improvement" criterion, so the
                            #  NO_IMPROVE_LIMIT / PERTURB_SD / RESTART_EVERY
                            #  knobs from transparent.R are intentionally absent)
SEED         <- 42          # master seed; per-combo seed = SEED + ci*10000,
                            # per-candidate seed = combo_seed + i (i.i.d. draws)

## --- objective WEIGHTS / numerics -------------------------------------------
W_DEATHS <- 2                   # fatal weight (matches 031's "2*RMSE_deaths ...")
W_PREV   <- 1                   # non-fatal weight ("... + RMSE_prev")
EPS_REL  <- 1e-6                # denominator floor for relative error (no /0)

## --- probability-constraint numerics ----------------------------------------
# Epsilon buffer kept below 1 when capping/renormalising probabilities. 032 uses
# the same 0.005 buffer in its renormalisation; we reuse it for consistency.
TP_EPS <- 0.005

## --- calibration target window ----------------------------------------------
CAL_YEAR_START <- 2009          # GBD comparison start (031 uses year >= 2009)
CAL_YEAR_END   <- 2019          # GBD comparison end   (031 uses year <  2020)

## --- execution --------------------------------------------------------------
# run_calibration_par is defined in 00_run_model_indonesia.R; default TRUE here.
RUN_PAR   <- if (exists("run_calibration_par")) isTRUE(run_calibration_par) else TRUE
MAX_CORES <- 14                 # cap on workers (031 used 14)
N_OUT_CHUNKS <- 10              # number of adjusted_searo_part*.rds files (= 031)

#===============================================================================
# 1. LOAD INPUTS (same sources 031 reads -- this script never writes them)
#===============================================================================

## locations (written by 021_get_base_rates_indonesia.R)
locs <- readRDS(paste0(wd, "locs.rds"))
locs <- as.vector(locs$location)

## baseline transition probabilities (tps_inpt_part*.rds from 022)
## schema: age, sex, location, year, cause, ALL.mx, BG.mx.all, BG.mx,
##         PREVt0, DIS.mx.t0, Nx, IR, CF   (cause = full GBD names, ages 20-95,
##         years 2000-2019; contains disease causes only, no "All causes")
tps_files <- list.files(path = wd_data, pattern = "tps_inpt", full.names = TRUE)
b_rates   <- rbindlist(lapply(tps_files, function(f) { dt <- readRDS(f); setDT(dt); dt }),
                       use.names = TRUE, fill = TRUE)

## same defensive clamps 031 applies before calibrating
b_rates[CF >= 1, CF := 0.99]
b_rates[IR >= 1, IR := 0.99]
b_rates[CF < 0,  CF := 0]
b_rates[IR < 0,  IR := 0]

## keep a frozen copy of the INPUT for end-of-run schema / row-count validation
tps_input_cols <- copy(names(b_rates))
tps_input_nrow <- nrow(b_rates)

## incoming age-20 population (UNWPP 2024), exactly as 031 prepares it
pop20 <- fread(paste0(wd_data, "PopulationsAge20_full.csv"))
pop20 <- pop20[location %in% locs]
pop20 <- pop20[year_id >= CAL_YEAR_START & year_id <= 2050]
setnames(pop20, c("year_id", "Nx"), c("year", "Nx20"))

#===============================================================================
# 2. GBD CALIBRATION TARGETS  (Deaths + Prevalence, Number, by 5-year age group)
#    Mirrors 031 lines 237-321; produces one tidy targets table.
#===============================================================================

gbd <- readRDS(paste0(wd_raw, "GBD/", "temp_1baseline_rates_gbd23.rds"))
setDT(gbd)
if ("upper" %in% names(gbd)) gbd[, upper := NULL]
if ("lower" %in% names(gbd)) gbd[, lower := NULL]

gbd <- gbd[cause_name %in% dx_include]
setnames(gbd,
         c("sex_name", "age_name", "cause_name", "measure_name", "metric_name", "location_name"),
         c("sex",      "age",      "cause",      "measure",      "metric",      "location"))

## counts only (Number), Deaths + Prevalence only, calibration window, our locs
gbd <- gbd[metric == "Number" &
             measure %in% c("Deaths", "Prevalence") &
             location %in% locs &
             year >= CAL_YEAR_START & year <= CAL_YEAR_END]

## long -> wide so each (location, sex, cause, age.group, year) has Deaths + Prev
gbd <- dcast(gbd, location + sex + cause + age + year ~ measure, value.var = "val")

## drop "All causes" (no per-cause target; 031 also excludes it from the join)
targets <- gbd[cause != cause_map[["all"]],
               .(location, sex, cause, age, year,
                 gbdDeaths = Deaths, gbdPrev = Prevalence)]
setkey(targets, location, sex, cause)

#===============================================================================
# 3. HELPER FUNCTIONS
#===============================================================================

## ----------------------------------------------------------------------------
## Single-year age (20-95) -> 5-year GBD age-group label. Identical breaks to
## 031's age_match so model output aggregates onto GBD's age bins.
## ----------------------------------------------------------------------------
make_age_match <- function() {
      am <- data.table(age = 20:95)
      am[, age.group := fcase(
            age < 25, "20-24 years",
            age < 30, "25-29 years",
            age < 35, "30-34 years",
            age < 40, "35-39 years",
            age < 45, "40-44 years",
            age < 50, "45-49 years",
            age < 55, "50-54 years",
            age < 60, "55-59 years",
            age < 65, "60-64 years",
            age < 70, "65-69 years",
            age < 75, "70-74 years",
            age < 80, "75-79 years",
            age < 85, "80-84 years",
            age < 90, "85-89 years",
            age < 95, "90-94 years",
            default = "95+ years"
      )]
      am
}

## ----------------------------------------------------------------------------
## perturb_cvd_combo_random: draw ONE independent random multiplier vector.
##
## This is the heart of the random-search strategy. Unlike transparent.R's
## Gaussian-step-from-best, this draws every multiplier i.i.d. uniformly in
## [lo, hi] = [1 - hw, 1 + hw], with NO reference to any current/best vector.
## The returned vector has the same layout build_mtab() expects:
##   * "age_group": c(ir_mult[1..n_g], cf_mult[1..n_g])   (length 2*n_g)
##   * "combo"    : c(ir_mult, cf_mult)                    (length 2)
## Reproducibility: caller passes a per-candidate seed so each draw is
## independent yet exactly repeatable.
## ----------------------------------------------------------------------------
perturb_cvd_combo_random <- function(n_age_groups, granularity, lo, hi, seed = NULL) {
      if (!is.null(seed)) set.seed(seed)
      n_par <- if (granularity == "age_group") 2L * n_age_groups else 2L
      runif(n_par, lo, hi)            # i.i.d. uniform; independent of the best
}

## ----------------------------------------------------------------------------
## Enforce the probability / row constraints on a TP table (in place).
## PREFERENCE ORDER (per task spec):
##   (1) NA -> 0, clamp IR,CF,BG.mx into [0,1];
##   (2) PRIMARY: preserve BG.mx, cap the disease TP so IR+BG.mx<=1, CF+BG.mx<=1
##       i.e. IR,CF <= 1 - BG.mx - TP_EPS;
##   (3) FALLBACK (only when BG.mx itself leaves no room, BG.mx >= 1 - TP_EPS):
##       proportionally renormalise the disease TP AND BG.mx by their shares,
##       the share-shrink from 032_adjustments_indonesia.R lines 46-55. Rows
##       whose BG.mx was modified are flagged in `bg_modified` for audit.
## Adds an integer `bg_modified` column (0/1) used only for diagnostics.
## ----------------------------------------------------------------------------
enforce_tp_constraints <- function(dt, tp_eps = TP_EPS) {
      ## (1) NA -> 0 and basic [0,1] clamp -----------------------------------
      dt[is.na(IR),    IR := 0]
      dt[is.na(CF),    CF := 0]
      dt[is.na(BG.mx), BG.mx := 0]
      dt[IR < 0, IR := 0]; dt[IR > 1, IR := 1]
      dt[CF < 0, CF := 0]; dt[CF > 1, CF := 1]
      dt[BG.mx < 0, BG.mx := 0]

      if (!("bg_modified" %in% names(dt))) dt[, bg_modified := 0L]

      ## headroom for the disease TP given FIXED background mortality
      dt[, headroom := 1 - BG.mx - tp_eps]

      ## (2) PRIMARY: cap IR/CF into the headroom, leaving BG.mx untouched ----
      dt[headroom >= 0 & IR > headroom, IR := headroom]
      dt[headroom >= 0 & CF > headroom, CF := headroom]

      ## (3) FALLBACK: BG.mx alone leaves no room (headroom < 0). Must shrink
      ##     BG.mx. Proportional renormalisation a la 032 lines 46-55.
      ##     IR side:
      dt[headroom < 0 & (IR + BG.mx) > 1, `:=`(
            IR_new   = IR    / (IR + BG.mx) - tp_eps,
            BGmx_new = BG.mx / (IR + BG.mx) - tp_eps,
            bg_modified = 1L
      )]
      dt[!is.na(IR_new), `:=`(IR = pmax(IR_new, 0), BG.mx = pmax(BGmx_new, 0))]
      dt[, c("IR_new", "BGmx_new") := NULL]
      ##     CF side (uses possibly-updated BG.mx; for normal rows CF+BG.mx<=1
      ##     already holds after the primary cap, so this only fires in fallback)
      dt[(CF + BG.mx) > 1, `:=`(
            CF_new   = CF    / (CF + BG.mx) - tp_eps,
            BGmx_new = BG.mx / (CF + BG.mx) - tp_eps,
            bg_modified = 1L
      )]
      dt[!is.na(CF_new), `:=`(CF = pmax(CF_new, 0), BG.mx = pmax(BGmx_new, 0))]
      dt[, c("CF_new", "BGmx_new") := NULL]

      dt[, headroom := NULL]
      dt[]
}

## ----------------------------------------------------------------------------
## Build the multiplier table (age.group -> ir_mult, cf_mult) from a flat
## parameter vector, given the granularity. For "combo" the same scalar pair is
## broadcast to every age group; for "age_group" each group gets its own pair.
## ----------------------------------------------------------------------------
build_mtab <- function(par, age_groups, granularity) {
      if (granularity == "age_group") {
            n <- length(age_groups)
            data.table(age.group = age_groups,
                       ir_mult = par[1:n],
                       cf_mult = par[(n + 1):(2 * n)])
      } else {
            data.table(age.group = age_groups,
                       ir_mult = par[1],
                       cf_mult = par[2])
      }
}

## ----------------------------------------------------------------------------
## Apply the multipliers to a combo's TP rows (all years), then enforce the
## probability constraints. Returns a NEW data.table (baseline is never mutated).
## Age groups absent from mtab keep their baseline IR/CF (mult defaults to 1).
## ----------------------------------------------------------------------------
apply_multipliers <- function(combo_rates, mtab, age_match) {
      cr <- copy(combo_rates)
      cr <- merge(cr, age_match, by = "age",       all.x = TRUE)
      cr <- merge(cr, mtab,      by = "age.group", all.x = TRUE)
      cr[is.na(ir_mult), ir_mult := 1]
      cr[is.na(cf_mult), cf_mult := 1]
      cr[, IR := IR * ir_mult]
      cr[, CF := CF * cf_mult]
      cr[, c("age.group", "ir_mult", "cf_mult") := NULL]
      enforce_tp_constraints(cr)
      cr
}

## ----------------------------------------------------------------------------
## Project ONE location-sex-cause combo through the well-sick-dead recursion,
## 2009 -> 2019. This is 031's state.transition specialised to a single combo
## and rewritten with a KEYED update-join (br[upd, on=.(year,age), :=...])
## instead of 031's position-based assignment, which is fragile to row order.
##
## APPROXIMATION (documented): the full model couples causes through the shared
## population pool (all.mx2 = sum of deaths over ALL causes). Projecting a single
## cause means the pool is depleted only by THAT cause's deaths. Disease deaths
## are tiny relative to population, and targets are per-cause, so this is a
## negligible and intentional simplification that lets each combo calibrate
## independently (as in the cancer calibrate_by_combo template).
## ----------------------------------------------------------------------------
project_combo <- function(cr, pop_combo, y0 = CAL_YEAR_START, y1 = CAL_YEAR_END) {
      br <- merge(cr, pop_combo, by = c("year", "location", "sex", "age"), all.x = TRUE)
      ## incoming age-20 cohort each projection year comes from UNWPP (Nx20)
      br[age == 20 & year > y0, Nx := Nx20]
      br[, Nx20 := NULL]

      ## initial states: first calibration year (all ages) + age-20 (all years)
      br[year == y0 | age == 20, sick   := Nx * PREVt0]
      br[year == y0 | age == 20, dead   := Nx * DIS.mx.t0]
      br[year == y0 | age == 20, well   := Nx * (1 - (PREVt0 + ALL.mx))]
      br[year == y0 | age == 20, pop    := Nx]
      br[year == y0 | age == 20, all.mx := Nx * ALL.mx]

      ## same defensive caps 031 applies inside state.transition
      br[CF > 0.9, CF := 0.9]
      br[IR > 0.9, IR := 0.9]

      n_steps <- y1 - y0
      for (s in 1:n_steps) {
            yr <- y0 + s
            b2 <- br[year <= yr & year >= yr - 1]
            setorder(b2, sex, location, cause, age, year)   # ensure shift() reads year yr-1
            b2[, age2 := age + 1]

            ## sick' = sick*(1-(CF+BG.mx)) + well*IR
            b2[, sick2 := shift(sick) * (1 - (CF + BG.mx)) + shift(well) * IR,
               by = .(sex, location, cause, age)]
            b2[sick2 < 0, sick2 := 0]
            ## dead' = sick*CF
            b2[, dead2 := shift(sick) * CF, by = .(sex, location, cause, age)]
            b2[dead2 < 0, dead2 := 0]
            ## pop'  = pop - all.mx
            b2[, pop2 := shift(pop) - shift(all.mx), by = .(sex, location, cause, age)]
            b2[pop2 < 0, pop2 := 0]
            ## all.mx' = (disease deaths this combo) + background mortality of pool
            b2[, all.mx2 := sum(dead2), by = .(sex, location, year, age)]
            b2[, all.mx2 := all.mx2 + (pop2 * BG.mx.all)]
            b2[all.mx2 < 0, all.mx2 := 0]
            ## well' = pop' - all.mx' - sick'
            b2[, well2 := pop2 - all.mx2 - sick2]
            b2[well2 < 0, well2 := 0]

            upd <- b2[year == yr & age2 < 96,
                      .(age = age2, year, sick2, dead2, well2, pop2, all.mx2)]
            br[upd, on = .(year, age), `:=`(
                  sick   = i.sick2,
                  dead   = i.dead2,
                  well   = i.well2,
                  pop    = i.pop2,
                  all.mx = i.all.mx2
            )]
      }

      br[year >= y0, .(location, sex, cause, year, age, sick, dead)]
}

## ----------------------------------------------------------------------------
## Aggregate a single-age projection to 5-year age groups (model Prevalence =
## sum of sick, model Deaths = sum of dead) and join GBD targets.
## ----------------------------------------------------------------------------
proj_vs_targets <- function(proj, combo_targets, age_match) {
      m <- merge(proj, age_match, by = "age", all.x = TRUE)
      ms <- m[, .(Prevalence = sum(sick, na.rm = TRUE),
                  Deaths     = sum(dead, na.rm = TRUE)),
              by = .(location, sex, cause, year, age.group)]
      setnames(ms, "age.group", "age")
      j <- merge(combo_targets, ms,
                 by = c("location", "sex", "cause", "year", "age"), all.x = TRUE)
      j[is.na(Deaths),     Deaths := 0]
      j[is.na(Prevalence), Prevalence := 0]
      j
}

## ----------------------------------------------------------------------------
## Weighted RELATIVE squared error (the search objective). Deaths weighted 2x.
## ----------------------------------------------------------------------------
combo_error <- function(proj, combo_targets, age_match,
                        w_deaths = W_DEATHS, w_prev = W_PREV, eps = EPS_REL) {
      j <- proj_vs_targets(proj, combo_targets, age_match)
      j[, sum(
            w_deaths * ((Deaths     - gbdDeaths) / (gbdDeaths + eps))^2 +
            w_prev   * ((Prevalence - gbdPrev)   / (gbdPrev   + eps))^2,
            na.rm = TRUE)]
}

## ----------------------------------------------------------------------------
## ABSOLUTE RMSE diagnostics per location-sex-cause-age.group (for comparison
## with 031, which minimised 2*RMSE_deaths + RMSE_prev in count units).
## ----------------------------------------------------------------------------
combo_diag <- function(proj, combo_targets, age_match) {
      j <- proj_vs_targets(proj, combo_targets, age_match)
      j[, .(
            RMSE_deaths = sqrt(mean((Deaths     - gbdDeaths)^2, na.rm = TRUE)),
            RMSE_prev   = sqrt(mean((Prevalence - gbdPrev)^2,   na.rm = TRUE))
      ), by = .(location, sex, cause, age)]
}

## ----------------------------------------------------------------------------
## Calibrate ONE combo: PURE RANDOM SEARCH over the multiplier vector.
##
## Candidate 0 = baseline (all multipliers = 1), so the returned fit is never
## worse than baseline. For i = 1..n_iter we draw an INDEPENDENT i.i.d. uniform
## candidate via perturb_cvd_combo_random() -- crucially, the draw does NOT
## depend on best_par (contrast transparent.R's Gaussian step from the best).
## We retain the argmin error vector. "Keep best" is storage only; it never
## seeds the next draw. The only early stop is CONVERGE_TOL; there is no
## "no-improvement" stop because that is a hill-climb notion with no meaning
## for memoryless random search. Each candidate gets its own derived seed so
## the whole search is reproducible.
## ----------------------------------------------------------------------------
calibrate_one_combo_random <- function(combo_rates, pop_combo, combo_targets, age_match,
                                       granularity, hw, n_iter, converge_tol,
                                       w_deaths, w_prev, eps, seed) {

      lo <- 1 - hw; hi <- 1 + hw
      age_groups <- sort(unique(combo_targets$age))
      n_g   <- length(age_groups)
      n_par <- if (granularity == "age_group") 2L * n_g else 2L

      eval_par <- function(par) {
            mtab <- build_mtab(par, age_groups, granularity)
            cr2  <- apply_multipliers(combo_rates, mtab, age_match)
            proj <- project_combo(cr2, pop_combo)
            combo_error(proj, combo_targets, age_match, w_deaths, w_prev, eps)
      }

      ## candidate 0: baseline (multipliers = 1) -----------------------------
      best_par <- rep(1, n_par)
      best_err <- eval_par(best_par)
      base_err <- best_err                 # baseline error, for diagnostics
      n_eval   <- 1L

      ## candidates 1..n_iter: INDEPENDENT i.i.d. uniform draws ---------------
      for (it in 1:n_iter) {
            cand <- perturb_cvd_combo_random(n_g, granularity, lo, hi,
                                             seed = seed + it)   # i.i.d., reproducible
            err  <- eval_par(cand); n_eval <- n_eval + 1L

            if (err < best_err) {                 # keep best (storage only)
                  best_err <- err; best_par <- cand
            }
            if (best_err < converge_tol) break    # only early stop for random search
      }

      list(mtab = build_mtab(best_par, age_groups, granularity),
           best_err = best_err, base_err = base_err,
           n_eval = n_eval, n_par = n_par,
           hit_bound = any(best_par <= lo + 1e-9 | best_par >= hi - 1e-9))
}

## ----------------------------------------------------------------------------
## Drive one combo end-to-end: calibrate, bake multipliers into ALL years, and
## assemble its calibrated rows, factor records, and baseline-vs-calibrated
## diagnostics. Returns a list combined across combos after the parallel loop.
## ----------------------------------------------------------------------------
run_combo <- function(ci, combos, b_rates, pop20, targets, age_match) {
      loc <- combos$location[ci]; sx <- combos$sex[ci]; cse <- combos$cause[ci]

      cr <- b_rates[location == loc & sex == sx & cause == cse]
      pc <- pop20[location == loc & sex == sx]
      ct <- targets[location == loc & sex == sx & cause == cse]

      ## baseline (mult = 1) rows + projection, for diagnostics & as fallback
      base_rows <- enforce_tp_constraints(copy(cr))

      if (nrow(ct) == 0) {
            ## no GBD target for this combo -> keep baseline unchanged
            base_rows[, bg_modified := NULL]
            return(list(
                  rows = base_rows,
                  factors = data.table(location = loc, sex = sx, cause = cse,
                                       age.group = NA_character_, ir_mult = 1, cf_mult = 1,
                                       granularity = GRANULARITY),
                  diag = data.table(location = loc, sex = sx, cause = cse,
                                    age = NA_character_,
                                    RMSE_deaths_base = NA_real_, RMSE_prev_base = NA_real_,
                                    RMSE_deaths_cal = NA_real_,  RMSE_prev_cal = NA_real_),
                  err = data.table(location = loc, sex = sx, cause = cse,
                                   base_err = NA_real_, cal_err = NA_real_,
                                   n_eval = 0L, n_par = 0L, hit_bound = FALSE,
                                   bg_modified_rows = 0L)
            ))
      }

      fit <- calibrate_one_combo_random(cr, pc, ct, age_match,
                                        GRANULARITY, SEARCH_HALFWIDTH, N_ITER,
                                        CONVERGE_TOL, W_DEATHS, W_PREV, EPS_REL,
                                        SEED + ci * 10000L)

      cal_rows  <- apply_multipliers(cr, fit$mtab, age_match)
      bg_mod_n  <- sum(cal_rows$bg_modified)
      cal_rows[, bg_modified := NULL]

      ## diagnostics: project baseline and calibrated, compute absolute RMSE
      base_diag <- combo_diag(project_combo(base_rows, pc), ct, age_match)
      cal_diag  <- combo_diag(project_combo(cal_rows,  pc), ct, age_match)
      setnames(base_diag, c("RMSE_deaths", "RMSE_prev"),
               c("RMSE_deaths_base", "RMSE_prev_base"))
      setnames(cal_diag,  c("RMSE_deaths", "RMSE_prev"),
               c("RMSE_deaths_cal", "RMSE_prev_cal"))
      diag <- merge(base_diag, cal_diag,
                    by = c("location", "sex", "cause", "age"), all = TRUE)

      factors <- copy(fit$mtab)
      factors[, `:=`(location = loc, sex = sx, cause = cse, granularity = GRANULARITY)]
      setcolorder(factors, c("location", "sex", "cause", "age.group",
                             "ir_mult", "cf_mult", "granularity"))

      err <- data.table(location = loc, sex = sx, cause = cse,
                        base_err = fit$base_err, cal_err = fit$best_err,
                        n_eval = fit$n_eval, n_par = fit$n_par,
                        hit_bound = fit$hit_bound, bg_modified_rows = bg_mod_n)

      list(rows = cal_rows, factors = factors, diag = diag, err = err)
}

#===============================================================================
# 4. RUN CALIBRATION OVER ALL location-sex-cause COMBOS
#    Parallelised over COMBOS (not locations): Indonesia is a single location,
#    so location-level parallelism would use one core; combo-level uses all 12
#    (6 causes x 2 sexes) independent units.
#===============================================================================

age_match <- make_age_match()
combos    <- unique(b_rates[, .(location, sex, cause)])
n_combos  <- nrow(combos)

cat(sprintf("Random-search TP calibration: %d location-sex-cause combos | granularity = %s\n",
            n_combos, GRANULARITY))
cat(sprintf("Search range per multiplier: [%.2f, %.2f] | %d i.i.d. candidates/combo\n",
            1 - SEARCH_HALFWIDTH, 1 + SEARCH_HALFWIDTH, N_ITER))

worker_exports <- c(
      "combos", "b_rates", "pop20", "targets", "age_match",
      "make_age_match", "perturb_cvd_combo_random", "enforce_tp_constraints",
      "build_mtab", "apply_multipliers", "project_combo", "proj_vs_targets",
      "combo_error", "combo_diag", "calibrate_one_combo_random", "run_combo",
      "GRANULARITY", "SEARCH_HALFWIDTH", "N_ITER", "CONVERGE_TOL",
      "W_DEATHS", "W_PREV", "EPS_REL",
      "TP_EPS", "CAL_YEAR_START", "CAL_YEAR_END", "SEED"
)

if (RUN_PAR && n_combos > 1) {
      n_cores <- max(1L, min(MAX_CORES, n_combos, parallel::detectCores() - 1L))
      cat(sprintf("Running in parallel on %d cores...\n", n_cores))
      cl <- makeCluster(n_cores)
      registerDoParallel(cl)
      results <- foreach(ci = seq_len(n_combos),
                         .packages = c("data.table"),
                         .export   = worker_exports) %dopar% {
            setDTthreads(1)
            run_combo(ci, combos, b_rates, pop20, targets, age_match)
      }
      stopCluster(cl)
} else {
      cat("Running sequentially...\n")
      results <- lapply(seq_len(n_combos), function(ci) {
            res <- run_combo(ci, combos, b_rates, pop20, targets, age_match)
            cat(sprintf("  [%d/%d] %s | %s | %s : err %.3g -> %.3g\n",
                        ci, n_combos, combos$location[ci], combos$sex[ci], combos$cause[ci],
                        res$err$base_err, res$err$cal_err))
            res
      })
}

## --- consolidate ------------------------------------------------------------
calibrated   <- rbindlist(lapply(results, `[[`, "rows"),    use.names = TRUE, fill = TRUE)
factors_out  <- rbindlist(lapply(results, `[[`, "factors"), use.names = TRUE, fill = TRUE)
diag_out     <- rbindlist(lapply(results, `[[`, "diag"),    use.names = TRUE, fill = TRUE)
err_out      <- rbindlist(lapply(results, `[[`, "err"),     use.names = TRUE, fill = TRUE)

## restore the exact INPUT column set/order (drop-in schema match with 031)
setcolorder(calibrated, intersect(tps_input_cols, names(calibrated)))

#===============================================================================
# 5. WRITE OUTPUTS  (drop-in: adjusted_searo_part{1..10}.rds in wd_data)
#===============================================================================

n     <- nrow(calibrated)
chunk <- ceiling(n / N_OUT_CHUNKS)
for (i in 1:N_OUT_CHUNKS) {
      start <- (i - 1) * chunk + 1
      end   <- min(i * chunk, n)
      if (start > n) {                       # write empty tail chunk for parity
            saveRDS(calibrated[0], file = paste0(wd_data, "adjusted_searo_part", i, ".rds"))
            next
      }
      saveRDS(calibrated[start:end],
              file = paste0(wd_data, "adjusted_searo_part", i, ".rds"))
}

## calibrated multipliers (analogous to 032's adjustments2023_age.csv)
fwrite(factors_out, paste0(wd_data, "calibration_factors_random_tp.csv"))

## fit diagnostics: baseline-vs-calibrated absolute RMSE per combo x age.group,
## plus the per-combo weighted error and its % improvement over baseline.
err_pct <- copy(err_out[, .(location, sex, cause, base_err, cal_err,
                            n_eval, n_par, hit_bound, bg_modified_rows)])
err_pct[, pct_improvement := 100 * (base_err - cal_err) / pmax(base_err, EPS_REL)]
diag_full <- merge(diag_out, err_pct, by = c("location", "sex", "cause"), all.x = TRUE)
fwrite(diag_full, paste0(wd_data, "calibration_diagnostics_random_tp.csv"))

#===============================================================================
# 6. VALIDATION  (assert constraints; print baseline-vs-calibrated fit summary)
#===============================================================================

cat("\n", strrep("=", 70), "\nVALIDATION\n", strrep("=", 70), "\n", sep = "")

stopifnot(
      "IR contains NA"            = !anyNA(calibrated$IR),
      "CF contains NA"            = !anyNA(calibrated$CF),
      "BG.mx contains NA"         = !anyNA(calibrated$BG.mx),
      "IR outside [0,1]"          = calibrated[, all(IR >= 0 & IR <= 1)],
      "CF outside [0,1]"          = calibrated[, all(CF >= 0 & CF <= 1)],
      "IR + BG.mx > 1"            = calibrated[, all(IR + BG.mx <= 1 + 1e-9)],
      "CF + BG.mx > 1"            = calibrated[, all(CF + BG.mx <= 1 + 1e-9)],
      "row count != input"        = nrow(calibrated) == tps_input_nrow,
      "schema != input"           = setequal(names(calibrated), tps_input_cols)
)
cat("All probability/row constraints satisfied.\n")
cat(sprintf("Rows: %d (matches input: %d). Schema matches input: TRUE.\n",
            nrow(calibrated), tps_input_nrow))

bg_rows_total <- sum(err_out$bg_modified_rows, na.rm = TRUE)
if (bg_rows_total > 0) {
      cat(sprintf("NOTE: BG.mx was renormalised (fallback) on %d rows where ",
                  bg_rows_total),
          "BG.mx alone left no room for the disease TP. See bg_modified_rows ",
          "in calibration_diagnostics_random_tp.csv.\n", sep = "")
} else {
      cat("BG.mx preserved on ALL rows (no fallback renormalisation needed).\n")
}

## weighted-error improvement (search objective) -----------------------------
tot_base_err <- sum(err_out$base_err, na.rm = TRUE)
tot_cal_err  <- sum(err_out$cal_err,  na.rm = TRUE)
cat(sprintf("\nWeighted relative error (search objective), summed over combos:\n"))
cat(sprintf("  baseline = %.4g   calibrated = %.4g   reduction = %.1f%%\n",
            tot_base_err, tot_cal_err,
            100 * (tot_base_err - tot_cal_err) / max(tot_base_err, EPS_REL)))

## absolute RMSE improvement (comparable to 031's units) ---------------------
abs_summary <- diag_out[, .(
      RMSE_deaths_base = mean(RMSE_deaths_base, na.rm = TRUE),
      RMSE_deaths_cal  = mean(RMSE_deaths_cal,  na.rm = TRUE),
      RMSE_prev_base   = mean(RMSE_prev_base,   na.rm = TRUE),
      RMSE_prev_cal    = mean(RMSE_prev_cal,    na.rm = TRUE)
)]
cat("\nMean absolute RMSE across combo x age.group cells (counts):\n")
cat(sprintf("  Deaths     : baseline = %.1f -> calibrated = %.1f\n",
            abs_summary$RMSE_deaths_base, abs_summary$RMSE_deaths_cal))
cat(sprintf("  Prevalence : baseline = %.1f -> calibrated = %.1f\n",
            abs_summary$RMSE_prev_base, abs_summary$RMSE_prev_cal))

## per-combo table for quick scan -------------------------------------------
cat("\nPer-combo weighted error (baseline -> calibrated):\n")
print(err_out[order(cause, sex),
              .(location, sex, cause,
                base_err = round(base_err, 3),
                cal_err  = round(cal_err, 3),
                n_eval, hit_bound, bg_modified_rows)])

if (any(err_out$hit_bound, na.rm = TRUE)) {
      cat("\nWARNING: some combos hit the search bound -- consider widening ",
          "SEARCH_HALFWIDTH or switching GRANULARITY to \"age_group\".\n", sep = "")
}

cat("\nWrote:\n")
cat(sprintf("  %sadjusted_searo_part{1..%d}.rds\n", wd_data, N_OUT_CHUNKS))
cat(sprintf("  %scalibration_factors_random_tp.csv\n", wd_data))
cat(sprintf("  %scalibration_diagnostics_random_tp.csv\n", wd_data))
cat("\nReminder: set  run_adjustment_model <- FALSE  before sourcing ",
    "05_build_baseline_indonesia.R to avoid re-applying 032's old factors.\n", sep = "")

#===============================================================================
# 7. OPTIONAL: baseline-vs-calibrated comparison plot
#     The transparent script produces NO plots; this block is purely optional
#     and the pipeline does NOT depend on it. It only runs if ggplot2 is
#     available and MAKE_PLOTS is TRUE, so sourcing this file never fails for
#     lack of a plotting package. Mirrors diagnostic output 4 (baseline vs
#     calibrated absolute RMSE for Deaths and Prevalence, per combo).
#===============================================================================

MAKE_PLOTS <- if (exists("make_calibration_plots")) isTRUE(make_calibration_plots) else FALSE

if (MAKE_PLOTS && requireNamespace("ggplot2", quietly = TRUE)) {
      pdat <- melt(
            diag_out[, .(location, sex, cause, age,
                         Deaths_base = RMSE_deaths_base, Deaths_cal = RMSE_deaths_cal,
                         Prev_base   = RMSE_prev_base,   Prev_cal   = RMSE_prev_cal)],
            id.vars      = c("location", "sex", "cause", "age"),
            variable.name = "metric", value.name = "rmse"
      )
      pdat[, c("measure", "stage") := tstrsplit(metric, "_", fixed = TRUE)]
      p <- ggplot2::ggplot(
            pdat, ggplot2::aes(x = age, y = rmse, colour = stage, group = stage)) +
            ggplot2::geom_line() +
            ggplot2::facet_grid(measure ~ cause + sex, scales = "free_y") +
            ggplot2::labs(
                  title    = "Random-search calibration: baseline vs calibrated RMSE",
                  subtitle = sprintf("Indonesia | granularity = %s | %d candidates/combo",
                                     GRANULARITY, N_ITER),
                  x = "Age group", y = "Absolute RMSE (counts)", colour = "") +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      ggplot2::ggsave(paste0(wd_data, "calibration_fit_random_tp.png"),
                      p, width = 14, height = 7, dpi = 150)
      cat(sprintf("  %scalibration_fit_random_tp.png (optional plot)\n", wd_data))
} else if (MAKE_PLOTS) {
      cat("NOTE: make_calibration_plots is TRUE but ggplot2 is not installed; ",
          "skipping the optional comparison plot.\n", sep = "")
}
