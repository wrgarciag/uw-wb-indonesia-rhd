# ==============================================================================
# RHD secondary-prevention investment case: STANDARD LONG-TABLE OUTPUTS (A/B/C/D)
# scripts/07_make_outputs.R
#
# ------------------------------------------------------------------------------
# ROLE OF THIS SCRIPT
# ------------------------------------------------------------------------------
# Fed by 06's per-location runs. It reads every output/out_model/<location>.rds
# and compiles THREE tidy long tables (data.table throughout):
#
#   (1) rhd_model_long.csv/.rds  — the WELL-SICK-DEAD aggregate (contract UNCHANGED)
#         scenario, age, cause, sex, year, well, sick, newcases, dead, pop, all.mx,
#         intervention, location, eff_ir, eff_cf
#       * sick   = A + B + C + D (all RHD stages collapsed), from 06's $wsd.
#       * eff_ir = RRR multiplier on INCIDENCE, [0,1] (1 = no reduction). Always 1
#                  here (secondary prevention does not reduce incidence).
#       * eff_cf = SAP RRR multiplier on RHD-specific MORTALITY, [0,1]
#                  (= 1 - sap_rrr_rhd_death x effective_treatment_coverage). NOTE:
#                  this is the SAP mortality multiplier, NOT a progression rate.
#
#   (2) rhd_stage_model_long.csv — per-stage prevalence
#         scenario, age, sex, year, location, stage, cases, prevalence, prevalence_per_1000
#
#   (3) rhd_annual_flows.csv      — all A/B/C/D flows + volumes + surgery trace
#         scenario, age, sex, year, location, new_rhd_A, A_to_no_rhd, A_to_B,
#         B_to_A, B_to_C, C_to_B, C_to_D, D_to_C, rhd_deaths_[A..D],
#         other_deaths_[A..D], n_screened, n_diagnosed, n_on_optimal_treatment,
#         C_requiring_surgery, D_requiring_surgery, surgeries_C, surgeries_D,
#         total_surgeries
#
# ALL economic / benefit-cost / DALY content lives solely in 08.
#
#   INPUT : output/out_model/*.rds                 (from 06: $wsd + $stages)
#   OUTPUT: output/tables/rhd_model_long.csv/.rds, rhd_stage_model_long.csv,
#           rhd_annual_flows.csv
# ==============================================================================

library(data.table)

if (!exists("wd_outp")) wd_outp <- paste0(here::here("output"), "/")

IN_DIR  <- paste0(wd_outp, "out_model/")
OUT_DIR <- paste0(wd_outp, "tables/")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# EXACT aggregate output column contract (order matters)
OUT_COLS <- c("scenario", "age", "cause", "sex", "year",
              "well", "sick", "newcases", "dead", "pop", "all.mx",
              "intervention", "location", "eff_ir", "eff_cf")

# stage + flow output contracts
STAGE_COLS <- c("scenario", "age", "sex", "year", "location",
                "stage", "cases", "prevalence", "prevalence_per_1000")
FLOW_COLS  <- c("scenario", "age", "sex", "year", "location",
                "living_rhd_start",
                "new_rhd_A", "A_to_no_rhd", "A_to_B", "B_to_A", "B_to_C",
                "C_to_B", "C_to_D", "D_to_C",
                "rhd_deaths_A", "rhd_deaths_B", "rhd_deaths_C", "rhd_deaths_D",
                "other_deaths_A", "other_deaths_B", "other_deaths_C", "other_deaths_D",
                "n_screened", "n_diagnosed", "n_on_optimal_treatment",
                "C_requiring_surgery", "D_requiring_surgery",
                "surgeries_C", "surgeries_D", "total_surgeries")

# ------------------------------------------------------------------------------
# 1. LOAD per-location model outputs ($wsd aggregate + $stages stock/flow)
# ------------------------------------------------------------------------------
files <- list.files(IN_DIR, pattern = "\\.rds$", full.names = TRUE)
if (length(files) == 0)
  stop("No model outputs found in ", IN_DIR,
       ".\n  Run 06_run_prevention_model.R first.", call. = FALSE)

message(sprintf("── 07_make_outputs.R : reading %d location file(s) from %s ──",
                length(files), IN_DIR))

read_tbl <- function(which) rbindlist(lapply(files, function(f) {
  obj <- readRDS(f)
  if (is.null(obj[[which]]))
    stop("File ", basename(f), " has no $", which, " table (re-run 06).", call. = FALSE)
  as.data.table(obj[[which]])
}), use.names = TRUE, fill = TRUE)

wsd <- read_tbl("wsd")
stg <- read_tbl("stages")

# ------------------------------------------------------------------------------
# 2. AGGREGATE long table (exact contract)
# ------------------------------------------------------------------------------
missing <- setdiff(OUT_COLS, names(wsd))
if (length(missing))
  stop("Model $wsd is missing required column(s): ",
       paste(missing, collapse = ", "), call. = FALSE)

long <- wsd[, ..OUT_COLS]
setorder(long, location, scenario, sex, year, age)

# ------------------------------------------------------------------------------
# 3. STAGE long table (melt A/B/C/D -> stage; prevalence per 1,000)
# ------------------------------------------------------------------------------
stage_long <- melt(stg[, .(scenario, age, sex, year, location, pop, A, B, C, D)],
                   id.vars = c("scenario", "age", "sex", "year", "location", "pop"),
                   measure.vars = c("A", "B", "C", "D"),
                   variable.name = "stage", value.name = "cases")
stage_long[, `:=`(prevalence = fifelse(pop > 0, cases / pop, 0),
                  prevalence_per_1000 = fifelse(pop > 0, cases / pop * 1000, 0))]
stage_long[, stage := as.character(stage)]
stage_long <- stage_long[, ..STAGE_COLS]
setorder(stage_long, location, scenario, sex, year, age, stage)

# ------------------------------------------------------------------------------
# 4. FLOW table (A/B/C/D flows + program volumes + surgery trace)
# ------------------------------------------------------------------------------
flow <- stg[, ..FLOW_COLS]
setorder(flow, location, scenario, sex, year, age)

# ------------------------------------------------------------------------------
# 5. VALIDATION  (fail loudly before writing)
# ------------------------------------------------------------------------------
if (!identical(names(long), OUT_COLS))
  stop("Aggregate column order does not match the contract.", call. = FALSE)
if (!identical(names(stage_long), STAGE_COLS))
  stop("Stage column order does not match the contract.", call. = FALSE)
if (!identical(names(flow), FLOW_COLS))
  stop("Flow column order does not match the contract.", call. = FALSE)
if (anyNA(long) || anyNA(stage_long) || anyNA(flow))
  stop("A compiled table contains NA.", call. = FALSE)

# aggregate stock/flow columns non-negative
for (cc in c("well", "sick", "newcases", "dead", "pop", "all.mx"))
  if (any(long[[cc]] < -1e-6)) stop("Aggregate column '", cc, "' has negative values.", call. = FALSE)
# stage cases + every flow non-negative
if (any(stage_long$cases < -1e-6)) stop("Stage cases has negative values.", call. = FALSE)
for (cc in setdiff(FLOW_COLS, c("scenario", "age", "sex", "year", "location")))
  if (any(flow[[cc]] < -1e-6)) stop("Flow column '", cc, "' has negative values.", call. = FALSE)

# effect multipliers on [0,1]
for (cc in c("eff_ir", "eff_cf"))
  if (any(long[[cc]] < 0 | long[[cc]] > 1))
    stop("Column '", cc, "' outside [0,1].", call. = FALSE)

# well-sick-dead accounting: well + sick <= pop ; all-cause deaths >= RHD deaths
if (long[, any(well + sick > pop + 1e-3)])
  stop("well + sick exceeds pop somewhere.", call. = FALSE)
if (long[, any(all.mx + 1e-6 < dead)])
  stop("all-cause deaths < RHD deaths somewhere.", call. = FALSE)

# A + B + C + D must reconstruct the aggregate `sick` exactly
recon <- merge(
  stage_long[, .(sick_stg = sum(cases)), by = .(location, scenario, sex, age, year)],
  long[, .(location, scenario, sex, age, year, sick)],
  by = c("location", "scenario", "sex", "age", "year"))
if (recon[, max(abs(sick_stg - sick))] > 1e-6)
  stop("A+B+C+D does not reconstruct aggregate `sick`.", call. = FALSE)

# surgery volumes do not exceed the number requiring surgery
if (flow[, any(surgeries_C > C_requiring_surgery + 1e-6)] ||
    flow[, any(surgeries_D > D_requiring_surgery + 1e-6)] ||
    flow[, any(total_surgeries > C_requiring_surgery + D_requiring_surgery + 1e-6)])
  stop("Surgery volume exceeds the number requiring surgery.", call. = FALSE)

# total living RHD does not exceed population
if (long[, any(sick > pop + 1e-3)])
  stop("Total living RHD (sick) exceeds population somewhere.", call. = FALSE)

# completeness: every (location, scenario) covers the same full age x sex x year grid
grid_n <- long[, .N, by = .(location, scenario)]
if (uniqueN(grid_n$N) != 1L)
  stop("Age/sex/year grid is not identical across every location-scenario.", call. = FALSE)

# reference vs SAP: SAP averts RHD deaths (cumulative, per location)
chk <- dcast(long[, .(dead = sum(dead)), by = .(location, scenario)],
             location ~ scenario, value.var = "dead")
if ("sap" %in% names(chk) && "ref" %in% names(chk) &&
    chk[, any(sap > ref + 1e-6)])
  stop("SAP scenario has MORE cumulative RHD deaths than reference in some location.",
       call. = FALSE)

# ------------------------------------------------------------------------------
# 6. WRITE + report
# ------------------------------------------------------------------------------
fwrite(long,       paste0(OUT_DIR, "rhd_model_long.csv"))
saveRDS(long,      paste0(OUT_DIR, "rhd_model_long.rds"))
fwrite(stage_long, paste0(OUT_DIR, "rhd_stage_model_long.csv"))
saveRDS(stage_long, paste0(OUT_DIR, "rhd_stage_model_long.rds"))
fwrite(flow,       paste0(OUT_DIR, "rhd_annual_flows.csv"))
saveRDS(flow,      paste0(OUT_DIR, "rhd_annual_flows.rds"))

message("── Compiled standard long tables ──────────────────────")
message(sprintf("  rhd_model_long        : %s rows (aggregate WSD)",
                formatC(nrow(long), format = "d", big.mark = ",")))
message(sprintf("  rhd_stage_model_long  : %s rows (A/B/C/D prevalence)",
                formatC(nrow(stage_long), format = "d", big.mark = ",")))
message(sprintf("  rhd_annual_flows      : %s rows (flows + volumes + surgery trace)",
                formatC(nrow(flow), format = "d", big.mark = ",")))
message(sprintf("  locations: %s | scenarios: %s | years %d-%d | ages %d-%d",
                paste(unique(long$location), collapse = ", "),
                paste(unique(long$scenario), collapse = ", "),
                min(long$year), max(long$year), min(long$age), max(long$age)))

# concise headline: cumulative RHD deaths averted + end-year stage mix (SAP arm)
tot <- long[, .(dead = sum(dead)), by = scenario]
if (all(c("ref", "sap") %in% tot$scenario)) {
  averted <- tot[scenario == "ref", dead] - tot[scenario == "sap", dead]
  message(sprintf("  Cumulative RHD deaths averted (ref - sap): %s",
                  formatC(round(averted), format = "d", big.mark = ",")))
}
endyr <- max(stage_long$year)
mix <- stage_long[scenario == "sap" & year == endyr, .(cases = sum(cases)), by = stage]
mix[, pct := 100 * cases / sum(cases)]
message(sprintf("  SAP-arm stage mix in %d: %s", endyr,
                paste(sprintf("%s %.0f%%", mix$stage, mix$pct), collapse = " / ")))
message("── 07_make_outputs.R complete ─────────────────────────")
message("  Wrote: rhd_model_long.*, rhd_stage_model_long.csv, rhd_annual_flows.csv")
message("  Next: 08_economic_evaluation.R (benefit-cost & cost-effectiveness)")
