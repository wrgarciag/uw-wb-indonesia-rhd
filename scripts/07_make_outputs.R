# ==============================================================================
# RHD secondary-prevention investment case: STANDARD LONG-TABLE OUTPUT
# scripts/07_make_outputs.R
#
# ------------------------------------------------------------------------------
# ROLE OF THIS SCRIPT (refactored)
# ------------------------------------------------------------------------------
# Fed by the outputs of 06 (the per-location model runs). It reads every
# output/out_model/<location>.rds, takes the WELL-SICK-DEAD aggregate table (1)
# ($wsd), row-binds it across locations, and compiles ONE tidy long table with
# EXACTLY these columns, in this order:
#
#   scenario, age, cause, sex, year, well, sick, newcases, dead, pop, all.mx,
#   intervention, location, eff_ir, eff_cf
#
#   * sick   = sum of all RHD cases across the tunnel states (mild + severe + post),
#              as produced by 06's $wsd table.
#   * eff_ir = relative-risk-reduction multiplier on INCIDENCE (new RHD cases), on
#              the scale [0,1] (1 = no reduction, 0 = 100% reduction).
#   * eff_cf = relative-risk-reduction multiplier on RHD cause-specific MORTALITY,
#              same [0,1] scale.
#
# ALL economic / benefit-cost / DALY / monetary content has been REMOVED from this
# script — it now lives solely in 08_economic_evaluation.R.
#
#   INPUT : output/out_model/*.rds                 (from 06)
#   OUTPUT: output/tables/rhd_model_long.csv / .rds
# ==============================================================================

library(data.table)

if (!exists("wd_outp")) wd_outp <- paste0(here::here("output"), "/")

IN_DIR  <- paste0(wd_outp, "out_model/")
OUT_DIR <- paste0(wd_outp, "tables/")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# EXACT output column contract (order matters)
OUT_COLS <- c("scenario", "age", "cause", "sex", "year",
              "well", "sick", "newcases", "dead", "pop", "all.mx",
              "intervention", "location", "eff_ir", "eff_cf")

# ------------------------------------------------------------------------------
# 1. LOAD per-location model outputs (table 1 = $wsd)
# ------------------------------------------------------------------------------
files <- list.files(IN_DIR, pattern = "\\.rds$", full.names = TRUE)
if (length(files) == 0)
  stop("No model outputs found in ", IN_DIR,
       ".\n  Run 06_run_prevention_model.R first.", call. = FALSE)

message(sprintf("── 07_make_outputs.R : reading %d location file(s) from %s ──",
                length(files), IN_DIR))

wsd_list <- lapply(files, function(f) {
  obj <- readRDS(f)
  if (is.null(obj$wsd))
    stop("File ", basename(f), " has no $wsd table (re-run 06).", call. = FALSE)
  as.data.table(obj$wsd)
})
dt <- rbindlist(wsd_list, use.names = TRUE, fill = TRUE)

# ------------------------------------------------------------------------------
# 2. SELECT + ORDER the exact contract columns
# ------------------------------------------------------------------------------
missing <- setdiff(OUT_COLS, names(dt))
if (length(missing))
  stop("Model output is missing required column(s): ",
       paste(missing, collapse = ", "), call. = FALSE)

long <- dt[, ..OUT_COLS]
setorder(long, location, scenario, sex, year, age)

# ------------------------------------------------------------------------------
# 3. VALIDATION  (fail loudly before writing)
# ------------------------------------------------------------------------------
if (!identical(names(long), OUT_COLS))
  stop("Output column order does not match the contract.", call. = FALSE)
if (anyNA(long))
  stop("Compiled long table contains NA.", call. = FALSE)

# numeric stock/flow columns must be non-negative
for (cc in c("well", "sick", "newcases", "dead", "pop", "all.mx"))
  if (any(long[[cc]] < -1e-6)) stop("Column '", cc, "' has negative values.", call. = FALSE)

# effect multipliers on [0,1]
for (cc in c("eff_ir", "eff_cf"))
  if (any(long[[cc]] < 0 | long[[cc]] > 1))
    stop("Column '", cc, "' outside [0,1].", call. = FALSE)

# well-sick-dead accounting: well + sick <= pop ; all-cause deaths >= RHD deaths
if (long[, any(well + sick > pop + 1e-3)])
  stop("well + sick exceeds pop somewhere.", call. = FALSE)
if (long[, any(all.mx + 1e-6 < dead)])
  stop("all-cause deaths < RHD deaths somewhere.", call. = FALSE)

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
# 4. WRITE + report
# ------------------------------------------------------------------------------
fwrite(long,   paste0(OUT_DIR, "rhd_model_long.csv"))
saveRDS(long,  paste0(OUT_DIR, "rhd_model_long.rds"))

message("── Compiled standard long table ───────────────────────")
message(sprintf("  rows: %s | locations: %s | scenarios: %s | years %d-%d | ages %d-%d",
                formatC(nrow(long), format = "d", big.mark = ","),
                paste(unique(long$location), collapse = ", "),
                paste(unique(long$scenario), collapse = ", "),
                min(long$year), max(long$year), min(long$age), max(long$age)))
message("  columns: ", paste(names(long), collapse = ", "))

# concise headline: cumulative RHD deaths averted by scenario difference (all locs)
tot <- long[, .(dead = sum(dead), sick_end = sum(sick[year == max(year)])),
            by = scenario]
if (all(c("ref", "sap") %in% tot$scenario)) {
  averted <- tot[scenario == "ref", dead] - tot[scenario == "sap", dead]
  message(sprintf("  Cumulative RHD deaths averted (ref - sap): %s",
                  formatC(round(averted), format = "d", big.mark = ",")))
}
message("── 07_make_outputs.R complete ─────────────────────────")
message("  Wrote: ", OUT_DIR, "rhd_model_long.csv/.rds")
message("  Next: 08_economic_evaluation.R (benefit-cost & cost-effectiveness)")
