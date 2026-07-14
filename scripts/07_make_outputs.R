################################################################################
# INDONESIA INTEGRATED NCD MODEL — COMPILE OUTPUTS
# scripts/07_make_outputs.R
# ─────────────────────────────────────────────────────────────────────────────
# Reads scenario results from data/model/scenarios/ and produces clean summary
# tables and figures for reporting and slides.
#
# V1 outputs:
#   Tables — deaths averted, 40q30 trajectories, SDG 3.4 pace, slide audit tables
#   Figures — scenario comparisons, WPP/history bridge, deaths-averted by year and age
#
# PREREQUISITES: 06_run_scenarios.R must have been run first.
################################################################################

rm(list = ls())

if (!requireNamespace("here", quietly = TRUE))
  stop("Package 'here' is required. Install with: install.packages('here')", call. = FALSE)
source(here::here("R", "packages.R"))
library(here)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(scales)
})

source(here("R", "cause_registry.R"))   # SHARED_PARAMS

# ── PATHS ─────────────────────────────────────────────────────────────────────
SCEN_DIR    <- here("data", "model", "scenarios")
BASELINE_DIR <- here("data", "model", "baseline")
OUT_TABLES  <- here("outputs", "tables")
OUT_FIGURES <- here("outputs", "figures")
dir.create(OUT_TABLES,  recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIGURES, recursive = TRUE, showWarnings = FALSE)

PROJ_YEARS <- SHARED_PARAMS$run_years

################################################################################
# 1  LOAD SCENARIO RESULTS
################################################################################

message("\n── Checking prerequisite files ──────────────────────────────────────────")
required_files <- file.path(SCEN_DIR, c(
  "scenario_mx.csv", "scenario_deaths_averted.csv",
  "scenario_cum_averted.csv", "scenario_q4030.csv"
))
missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0)
  stop("Missing scenario output files — run scripts/06_run_scenarios.R first:\n",
       paste(" ", missing, collapse = "\n"))

message("\n── Loading scenario results ─────────────────────────────────────────────")

all_scenario_mx  <- read_csv(file.path(SCEN_DIR, "scenario_mx.csv"),
                             show_col_types = FALSE)
deaths_averted   <- read_csv(file.path(SCEN_DIR, "scenario_deaths_averted.csv"),
                             show_col_types = FALSE)
cum_averted      <- read_csv(file.path(SCEN_DIR, "scenario_cum_averted.csv"),
                             show_col_types = FALSE)
q4030_scen       <- read_csv(file.path(SCEN_DIR, "scenario_q4030.csv"),
                             show_col_types = FALSE)
q4030_baseline   <- read_csv(file.path(BASELINE_DIR, "baseline_q4030.csv"),
                             show_col_types = FALSE)

message("  scenario_mx rows          : ", comma(nrow(all_scenario_mx)))
message("  deaths_averted rows       : ", comma(nrow(deaths_averted)))
message("  q4030_scen rows           : ", comma(nrow(q4030_scen)))
message("  cum_averted periods       : ",
        paste(sort(unique(cum_averted$period)), collapse = ", "))
message("  Scenarios                 : ",
        paste(sort(unique(all_scenario_mx$scenario)), collapse = ", "))

################################################################################
# 2  SUMMARY TABLES
################################################################################

message("\n── Building summary tables ──────────────────────────────────────────────")

# ── T1: Deaths averted — cumulative 2025–2050 and 2025–2100 ──────────────────
# Both periods are cumulative from 2025, not sequential segments.
# This matches the README description and is what collaborators expect when
# comparing near-term (one generation) vs. full-century burden reduction.
PERIODS <- list(
  "2025-2050" = c(2025L, 2050L),
  "2025-2100" = c(2025L, 2100L)
)

t1_deaths_averted <- bind_rows(
  lapply(names(PERIODS), function(label) {
    yr_range <- PERIODS[[label]]
    deaths_averted |>
      filter(scenario != "baseline",
             year >= yr_range[1], year <= yr_range[2]) |>
      group_by(scenario, module) |>
      summarise(deaths_base      = sum(deaths_base,     na.rm = TRUE),
                deaths_averted   = sum(averted_display, na.rm = TRUE),
                deaths_averted_raw = sum(averted_raw,   na.rm = TRUE),  # signed; for QA
                .groups = "drop") |>
      mutate(period        = label,
             pct_reduction = round(100 * deaths_averted / pmax(deaths_base, 1), 2))
  })
) |>
  select(period, scenario, module, deaths_base,
         deaths_averted, deaths_averted_raw, pct_reduction) |>
  arrange(period, scenario, module)

write_csv(t1_deaths_averted, file.path(OUT_TABLES, "t1_deaths_averted_by_period.csv"))
message("  t1_deaths_averted_by_period.csv ✓  (periods: 2025-2050, 2025-2100 cumulative)")

# ── T2: 40q30 at key years by scenario and cause ─────────────────────────────
t2_q4030 <- q4030_scen |>
  filter(year %in% c(2025, 2030, 2040, 2050, 2075, 2100)) |>
  select(scenario, year, sex, cause, q4030_pct) |>
  pivot_wider(names_from = year, values_from = q4030_pct,
              names_prefix = "q_") |>
  arrange(scenario, sex, cause)

write_csv(t2_q4030, file.path(OUT_TABLES, "t2_q4030_by_scenario_year.csv"))
message("  t2_q4030_by_scenario_year.csv ✓")

# ── T3: SDG 3.4 pace ──────────────────────────────────────────────────────────
#
# Headline V1 SDG proxy: female all-cause 40q30.
# Secondary outputs: male all-cause 40q30 and cause-specific 40q30 by sex.
#
# SDG 3.4 benchmark context:
#   Classic benchmark : ((2/3)^(1/15) − 1) × 100 ≈ −2.7%/yr — the AROC
#     needed to achieve a one-third reduction over the full 2015–2030 window.
#     This is a mathematical constant independent of the data level.
#   Catch-up benchmark: depends on TWO model-/data-specific quantities:
#     (1) GBD Indonesia female all-cause 40q30 at 2015 — the SDG baseline
#     (2) Model baseline female all-cause 40q30 at 2025 — current level
#     Formula: ((2/3 × q4030_2015) / q4030_2025)^(1/5) − 1
#     So the catch-up rate is NOT a fixed constant: it changes with each
#     model run as q4030_2025 updates.

# Step 1: compute q4030_2015 from GBD all-cause mortality rates.
# Use the 8 five-year age groups spanning ages 30–69 (midpoints 32:5:67).
GBD_AC_FILE <- here("data", "gbd", "gbd_allcause_mx.csv")
if (!file.exists(GBD_AC_FILE))
  stop("gbd_allcause_mx.csv not found — run scripts/01_prepare_gbd_inputs.R first.")

gbd_q4030_2015 <- read_csv(GBD_AC_FILE, show_col_types = FALSE) |>
  filter(location_id == 11L, year == 2015L, sex == "Female",
         age_mid >= 30, age_mid <= 67) |>
  arrange(age_mid) |>
  mutate(
    mx_per_person = mx_all / 1e5,            # GBD rate per 100k → per person per year
    q5            = 1 - exp(-5 * mx_per_person)  # 5-year death probability
  ) |>
  summarise(q4030 = 1 - prod(1 - q5)) |>    # survival product over 8 groups (ages 30–69)
  pull(q4030)

# Step 2: model female all-cause 40q30 at 2025 (baseline scenario).
model_q4030_2025_f <- q4030_scen |>
  filter(year == 2025L, sex == "Female", cause == "All-cause",
         scenario == "baseline") |>
  pull(q4030_pct) / 100   # convert % → proportion

# Step 3: compute the proper catch-up AROC.
# SDG 3.4 target for Indonesia: one-third reduction from GBD 2015 level by 2030.
sdg_target_2030_f <- (2/3) * gbd_q4030_2015

SDG_CATCHUP_AROC <- round(100 * ((sdg_target_2030_f / model_q4030_2025_f)^(1/5) - 1), 1)

# Classic benchmark is a pure mathematical constant — does not depend on data level.
SDG_CLASSIC_AROC <- round(100 * ((2/3)^(1/15) - 1), 1)   # ≈ −2.7%/yr

message(sprintf("  GBD female 40q30 at 2015 : %.1f%%", gbd_q4030_2015 * 100))
message(sprintf("  Model female 40q30 at 2025: %.1f%%", model_q4030_2025_f * 100))
message(sprintf("  SDG 3.4 target by 2030    : %.1f%%", sdg_target_2030_f * 100))
message(sprintf("  SDG catch-up AROC (data-driven) : %.1f%%/yr", SDG_CATCHUP_AROC))
message(sprintf("  SDG classic AROC (mathematical) : %.1f%%/yr", SDG_CLASSIC_AROC))

t3_sdg <- q4030_scen |>
  filter(cause == "All-cause") |>
  left_join(
    q4030_scen |>
      filter(cause == "All-cause", year == 2025L) |>
      select(scenario, sex, q4030_2025 = q4030_pct),
    by = c("scenario", "sex")
  ) |>
  filter(year %in% c(2025L, 2030L, 2040L, 2050L)) |>
  mutate(
    # AROC undefined at t = 2025 (zero denominator in exponent) → NA_real_
    aroc_pct_per_year = if_else(
      year > 2025L,
      round(100 * ((q4030_pct / pmax(q4030_2025, 0.001))^(1 / (year - 2025)) - 1), 3),
      NA_real_
    ),
    pct_change_vs_2025         = round((q4030_pct - q4030_2025) / pmax(q4030_2025, 0.001) * 100, 2),
    sdg_catchup_aroc_benchmark = SDG_CATCHUP_AROC
  ) |>
  select(scenario, sex, year, q4030_pct, q4030_2025,
         aroc_pct_per_year, pct_change_vs_2025, sdg_catchup_aroc_benchmark) |>
  arrange(scenario, sex, year)

write_csv(t3_sdg, file.path(OUT_TABLES, "t3_sdg_3_4_pace.csv"))
message("  t3_sdg_3_4_pace.csv ✓")
message("  Primary metric : aroc_pct_per_year (AROC %/yr from 2025; NA at t=2025)")
message("  Benchmark      : catch-up ", SDG_CATCHUP_AROC, "%/yr (2025→2030, not classic 2015-base)")
message("  Headline proxy : female all-cause 40q30")

################################################################################
# 3  FIGURES
################################################################################

message("\n── Generating output figures ────────────────────────────────────────────")

theme_out <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        strip.background  = element_rect(fill = "grey92", colour = NA),
        strip.text        = element_text(face = "bold", size = 9),
        plot.title        = element_text(face = "bold", size = 12),
        plot.subtitle     = element_text(size = 9, colour = "grey40"),
        legend.position   = "bottom")

SCEN_COLS <- c(
  baseline    = "#888888",
  bp_fast     = "#C0392B", bp_slow      = "#E74C3C",
  statins_fast = "#2980B9", statins_slow = "#3498DB",
  sodium      = "#27AE60", tfa          = "#F39C12",
  diabetes_bp = "#8E44AD",
  all_fast    = "#E67E22", all_slow     = "#D35400",
  cancer_fast = "#16A085", cancer_slow  = "#1ABC9C"
)

CAUSE_LABS <- c(
  ihd             = "IHD",
  ischemic_stroke = "Ischaemic stroke",
  ich             = "ICH",
  hhd             = "HHD",
  cervical_ca     = "Cervical cancer"
)

# ── F1: All-cause 40q30 under all CVD scenarios — both sexes ─────────────────
cvd_scenarios <- c("baseline", "bp_fast", "bp_slow", "statins_fast",
                   "all_fast", "all_slow")

p_f1 <- q4030_scen |>
  filter(scenario %in% cvd_scenarios, cause == "All-cause") |>
  ggplot(aes(x = year, y = q4030_pct, colour = scenario,
             linetype = sex)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = SCEN_COLS[cvd_scenarios]) +
  scale_linetype_manual(values = c(Female = "solid", Male = "dashed")) +
  labs(title    = "F1: All-cause 40q30 under CVD intervention scenarios",
       subtitle = "Indonesia V1 | 2025–2100",
       x = NULL, y = "40q30 (%)", colour = NULL, linetype = "Sex") +
  theme_out
ggsave(file.path(OUT_FIGURES, "f1_allcause_q4030_cvd_scenarios.png"),
       p_f1, width = 10, height = 5, dpi = 150)

# ── F2: Cumulative deaths averted 2025–2050 — stacked bar by cause ────────────
p_f2 <- cum_averted |>
  filter(period == "2025-2050",                    # must filter: file now has both periods
         module %in% names(CVD_CAUSE_MAP), scenario != "baseline",
         scenario %in% cvd_scenarios) |>
  mutate(
    scenario_lab = factor(scenario, levels = rev(cvd_scenarios)),
    module_lab   = factor(recode(module, !!!CAUSE_LABS[names(CVD_CAUSE_MAP)]),
                          levels = c("IHD", "Ischaemic stroke", "ICH", "HHD"))
  ) |>
  ggplot(aes(x = scenario_lab, y = cum_averted / 1e3, fill = module_lab)) +
  geom_col(position = "stack") +
  coord_flip() +
  scale_fill_manual(values = c(IHD = "#C0392B", "Ischaemic stroke" = "#E74C3C",
                               ICH = "#F39C12", HHD = "#E67E22")) +
  labs(title    = "F2: Cumulative CVD deaths averted 2025–2050",
       subtitle = "Indonesia V1 | stacked by cause",
       x = NULL, y = "Deaths averted (thousands)", fill = "Cause") +
  theme_out
ggsave(file.path(OUT_FIGURES, "f2_cum_averted_2025_2050.png"),
       p_f2, width = 10, height = 5, dpi = 150)

# ── F3: SDG 3.4 pace — AROC from 2025 ────────────────────────────────────────
p_f3 <- t3_sdg |>
  filter(scenario %in% cvd_scenarios, sex == "Female", year > 2025) |>
  ggplot(aes(x = year, y = aroc_pct_per_year, colour = scenario)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_hline(yintercept = SDG_CATCHUP_AROC, colour = "grey40",
             linetype = "dashed", linewidth = 0.6) +
  annotate("text", x = 2031, y = SDG_CATCHUP_AROC + 0.4,
           label = paste0("SDG 3.4 catch-up: ", SDG_CATCHUP_AROC, "%/yr (2025→2030)"),
           size = 3, colour = "grey40", hjust = 0) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = SCEN_COLS[cvd_scenarios]) +
  scale_x_continuous(breaks = c(2030, 2040, 2050)) +
  labs(title    = "F3: SDG 3.4 pace — annualised rate of change in 40q30 from 2025 (Female)",
       subtitle = "Indonesia V1 | dashed = rate required to cut 40q30 by 1/3 between 2025 and 2030",
       x = NULL, y = "AROC (%/year from 2025)", colour = NULL) +
  theme_out
ggsave(file.path(OUT_FIGURES, "f3_sdg_3_4_pace_female.png"),
       p_f3, width = 10, height = 5, dpi = 150)

message("  f1_allcause_q4030_cvd_scenarios.png ✓")
message("  f2_cum_averted_2025_2050.png ✓")
message("  f3_sdg_3_4_pace_female.png ✓")

################################################################################
# 4  DECK FIGURES
# Saved to outputs/figures/ with the exact filenames the LaTeX deck references.
# docs/slides/indonesia_ncd_v1.tex uses \graphicspath{{../../outputs/figures/}}
# so these are always current without a manual copy step.
################################################################################

message("\n── Generating deck figures ──────────────────────────────────────────────")

# Shared theme for deck figures
theme_deck <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major  = element_line(colour = "grey92"),
        plot.title        = element_text(face = "bold", size = 13),
        plot.subtitle     = element_text(size = 10, colour = "grey45"),
        plot.caption      = element_text(size = 8,  colour = "grey55",
                                         hjust = 0),
        axis.text         = element_text(size = 10),
        axis.title        = element_text(size = 10),
        legend.position   = "right",
        legend.text       = element_text(size = 9))

CVD_CAUSE_COLS <- c(
  ihd             = "#8B2500",
  ischemic_stroke = "#E74C3C",
  ich             = "#F39C12",
  hhd             = "#E67E22"
)
CVD_CAUSE_LABS <- c(
  ihd             = "IHD",
  ischemic_stroke = "Isch. stroke",
  ich             = "ICH",
  hhd             = "HHD"
)

# ── D1: baseline_cvd_burden.png ──────────────────────────────────────────────
# Horizontal bar: cumulative 4-CVD baseline deaths 2025-2050 by cause.
# Use deaths_base from a non-baseline scenario — all scenarios share the same
# deaths_base (it is the baseline reference). "baseline" is never a row in
# scenario_cum_averted.csv because averted = 0 for the baseline itself.
d1_data <- deaths_averted |>
  filter(scenario == "all_fast",
         year >= 2025, year <= 2050,
         module %in% names(CVD_CAUSE_LABS)) |>
  group_by(module) |>
  summarise(cum_baseline_deaths = sum(deaths_base, na.rm = TRUE),
            .groups = "drop") |>
  mutate(
    cause_lab = factor(recode(module, !!!CVD_CAUSE_LABS),
                       levels = rev(c("IHD", "Isch. stroke", "ICH", "HHD"))),
    label     = paste0(round(cum_baseline_deaths / 1e6, 2), "M")
  )

p_d1 <- d1_data |>
  ggplot(aes(x = cum_baseline_deaths / 1e6, y = cause_lab,
             fill = module)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = label), hjust = -0.1, size = 3.8, fontface = "bold") +
  scale_fill_manual(values = CVD_CAUSE_COLS, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18)),
                     labels = scales::label_number(suffix = "M")) +
  labs(title    = "Projected four-CVD deaths without intervention, Indonesia 2025\u20132050",
       subtitle = "IHD dominates cumulative baseline deaths under WPP \u00d7 projected GBD cause fractions",
       x        = "Cumulative deaths 2025\u20132050",
       y        = NULL,
       caption  = "Source: Indonesia NCD V1 baseline | Projected GBD cause fractions | WPP 2024 demographic engine") +
  theme_deck +
  theme(legend.position = "none")

ggsave(file.path(OUT_FIGURES, "baseline_cvd_burden.png"),
       p_d1, width = 9, height = 5, dpi = 180)
message("  baseline_cvd_burden.png ✓")

# ── D2: deaths_averted.png ────────────────────────────────────────────────────
# Stacked horizontal bar: cumulative deaths averted 2025-2050 by scenario and cause.
deck_scenarios_averted <- c("all_fast", "all_slow", "bp_fast", "bp_slow",
                            "diabetes_bp", "sodium", "statins_fast", "tfa")
deck_scenario_labs <- c(
  all_fast    = "All CVD (fast)",  all_slow    = "All CVD (slow)",
  bp_fast     = "BP control (fast)", bp_slow     = "BP control (slow)",
  diabetes_bp = "Diabetes BP",    sodium      = "Sodium reduction",
  statins_fast = "Statins (fast)", tfa         = "TFA elimination"
)

d2_data <- cum_averted |>
  filter(period == "2025-2050",
         scenario %in% deck_scenarios_averted,
         module   %in% names(CVD_CAUSE_LABS)) |>
  mutate(
    scenario_lab = factor(recode(scenario, !!!deck_scenario_labs),
                          levels = rev(recode(deck_scenarios_averted,
                                              !!!deck_scenario_labs))),
    cause_lab    = factor(recode(module, !!!CVD_CAUSE_LABS),
                          levels = c("HHD", "IHD", "ICH", "Isch. stroke"))
  )

p_d2 <- d2_data |>
  ggplot(aes(x = cum_averted / 1e6, y = scenario_lab, fill = cause_lab)) +
  geom_col(width = 0.7, position = "stack") +
  scale_fill_manual(values = c("IHD" = "#8B2500", "ICH" = "#F39C12",
                               "HHD" = "#E67E22", "Isch. stroke" = "#E74C3C")) +
  scale_x_continuous(labels = scales::label_number(suffix = "M")) +
  labs(title    = "Cumulative 4-CVD deaths averted 2025\u20132050",
       subtitle = paste0("All CVD (fast) averts ",
                         round(sum(d2_data$cum_averted[
                           d2_data$scenario == "all_fast"]) / 1e6, 2),
                         "M \u2014 stroke and ICH account for most gains"),
       x        = "Deaths averted (millions)",
       y        = NULL,
       fill     = "Cause",
       caption  = "Source: Indonesia NCD V1 | 'All CVD (fast)' = BP + statins + sodium + TFA") +
  theme_deck

ggsave(file.path(OUT_FIGURES, "deaths_averted.png"),
       p_d2, width = 10, height = 6, dpi = 180)
message("  deaths_averted.png ✓")


# ── D2A: deaths_averted_by_year_area.png ─────────────────────────────────────
# Annual deaths averted by year, stacked by cause, for the all-fast CVD package.
# This shows the scale-up profile directly: fast scenarios reach full target by
# 2030, after which changes reflect WPP population/mortality dynamics and the
# projected cause-fraction baseline.
AVERTED_FOCUS_SCENARIO     <- "all_fast"
AVERTED_FOCUS_SCENARIO_LAB <- "All CVD (fast)"
AVERTED_YEAR_RANGE         <- c(2025L, 2050L)
AVERTED_AGE_RANGE          <- c(30L, 100L)

cause_fill_labs <- c(
  "IHD"          = "#8B2500",
  "Isch. stroke" = "#E74C3C",
  "ICH"          = "#F39C12",
  "HHD"          = "#E67E22"
)

cause_stack_levels <- c("HHD", "IHD", "ICH", "Isch. stroke")

d2_year_data <- deaths_averted |>
  filter(
    scenario == AVERTED_FOCUS_SCENARIO,
    year >= AVERTED_YEAR_RANGE[1], year <= AVERTED_YEAR_RANGE[2],
    module %in% names(CVD_CAUSE_LABS)
  ) |>
  group_by(year, module) |>
  summarise(
    deaths_averted     = sum(averted_display, na.rm = TRUE),
    deaths_averted_raw = sum(averted_raw,     na.rm = TRUE),
    deaths_base        = sum(deaths_base,     na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    cause_lab = factor(recode(module, !!!CVD_CAUSE_LABS),
                       levels = cause_stack_levels)
  ) |>
  arrange(year, cause_lab)

write_csv(
  d2_year_data,
  file.path(OUT_TABLES, "t1b_deaths_averted_annual_by_cause_all_fast.csv")
)

all_fast_averted_2025_2050 <- sum(d2_year_data$deaths_averted, na.rm = TRUE)

p_d2_year <- d2_year_data |>
  ggplot(aes(x = year, y = deaths_averted / 1e3, fill = cause_lab)) +
  geom_area(alpha = 0.95, colour = "white", linewidth = 0.12) +
  scale_fill_manual(values = cause_fill_labs, breaks = names(cause_fill_labs)) +
  scale_x_continuous(breaks = seq(2025, 2050, by = 5)) +
  scale_y_continuous(labels = scales::label_number(suffix = "k"),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = "Annual 4-CVD deaths averted by cause, 2025\u20132050",
    subtitle = paste0(
      AVERTED_FOCUS_SCENARIO_LAB, " | Stacked annual deaths averted | ",
      "cumulative total = ", round(all_fast_averted_2025_2050 / 1e6, 2), "M"
    ),
    x        = NULL,
    y        = "Annual deaths averted",
    fill     = "Cause",
    caption  = paste0(
      "Source: Indonesia NCD V1 | Deaths averted = baseline minus scenario | ",
      "negative annual cells clipped to zero for display"
    )
  ) +
  theme_deck

ggsave(file.path(OUT_FIGURES, "deaths_averted_by_year_area.png"),
       p_d2_year, width = 10, height = 6, dpi = 180)
message("  deaths_averted_by_year_area.png ✓")

# ── D2B: deaths_averted_by_age_area.png ──────────────────────────────────────
# Cumulative deaths averted by single-year age, stacked by cause, for the
# all-fast CVD package. This is calculated from scenario_mx.csv so age-specific
# results are preserved; scenario_deaths_averted.csv is already aggregated.
required_age_cols <- c("scenario", "module", "year", "sex", "age", "deaths_cause")
missing_age_cols <- setdiff(required_age_cols, names(all_scenario_mx))
if (length(missing_age_cols) > 0) {
  stop(
    "Cannot build deaths_averted_by_age_area.png because scenario_mx.csv is missing: ",
    paste(missing_age_cols, collapse = ", "),
    call. = FALSE
  )
}

d2_age_wide <- all_scenario_mx |>
  filter(
    scenario %in% c("baseline", AVERTED_FOCUS_SCENARIO),
    module %in% names(CVD_CAUSE_LABS),
    year >= AVERTED_YEAR_RANGE[1], year <= AVERTED_YEAR_RANGE[2],
    age >= AVERTED_AGE_RANGE[1], age <= AVERTED_AGE_RANGE[2]
  ) |>
  group_by(scenario, module, year, sex, age) |>
  summarise(deaths = sum(deaths_cause, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(names_from = scenario, values_from = deaths)

missing_baseline_age <- d2_age_wide |>
  filter(is.na(baseline) | is.na(.data[[AVERTED_FOCUS_SCENARIO]])) |>
  distinct(module, year, sex, age)

if (nrow(missing_baseline_age) > 0) {
  stop(
    "Cannot build age-specific deaths-averted plot: missing baseline or focus-scenario rows. First rows:\n",
    paste(utils::capture.output(print(utils::head(missing_baseline_age, 20))), collapse = "\n"),
    call. = FALSE
  )
}

d2_age_data <- d2_age_wide |>
  mutate(
    deaths_averted_raw = baseline - .data[[AVERTED_FOCUS_SCENARIO]],
    deaths_averted     = pmax(deaths_averted_raw, 0)
  ) |>
  group_by(age, module) |>
  summarise(
    deaths_averted     = sum(deaths_averted,     na.rm = TRUE),
    deaths_averted_raw = sum(deaths_averted_raw, na.rm = TRUE),
    deaths_base        = sum(baseline,           na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    cause_lab = factor(recode(module, !!!CVD_CAUSE_LABS),
                       levels = cause_stack_levels)
  ) |>
  arrange(age, cause_lab)

write_csv(
  d2_age_data,
  file.path(OUT_TABLES, "t1c_deaths_averted_by_age_all_fast_2025_2050.csv")
)

p_d2_age <- d2_age_data |>
  ggplot(aes(x = age, y = deaths_averted / 1e3, fill = cause_lab)) +
  geom_area(alpha = 0.95, colour = "white", linewidth = 0.10) +
  scale_fill_manual(values = cause_fill_labs, breaks = names(cause_fill_labs)) +
  scale_x_continuous(breaks = seq(30, 100, by = 10)) +
  scale_y_continuous(labels = scales::label_number(suffix = "k"),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = "Cumulative 4-CVD deaths averted by age and cause, 2025\u20132050",
    subtitle = paste0(
      AVERTED_FOCUS_SCENARIO_LAB, " | Both sexes combined | Single-year ages ",
      AVERTED_AGE_RANGE[1], "\u2013", AVERTED_AGE_RANGE[2]
    ),
    x        = "Age",
    y        = "Cumulative deaths averted",
    fill     = "Cause",
    caption  = paste0(
      "Source: Indonesia NCD V1 | Calculated from scenario_mx.csv age-specific baseline and scenario deaths | ",
      "negative cells clipped to zero for display"
    )
  ) +
  theme_deck

ggsave(file.path(OUT_FIGURES, "deaths_averted_by_age_area.png"),
       p_d2_age, width = 10, height = 6, dpi = 180)
message("  deaths_averted_by_age_area.png ✓")

# ── D3: q4030_trajectory.png ──────────────────────────────────────────────────
# Female all-cause 40q30 trajectory under CVD scenarios, 2025-2055.
deck_scens_q40 <- c("baseline", "all_fast", "all_slow", "bp_fast")
deck_scen_labs_q40 <- c(
  all_fast = "All CVD (fast)", all_slow = "All CVD (slow)",
  baseline = "Baseline",       bp_fast  = "BP control (fast)"
)
deck_line_types <- c("All CVD (fast)" = "solid", "All CVD (slow)" = "dashed",
                     "Baseline"       = "solid", "BP control (fast)" = "dashed")
deck_line_cols  <- c("All CVD (fast)" = "#1F6BAE", "All CVD (slow)" = "#74B9E7",
                     "Baseline"       = "#2C3E50",  "BP control (fast)" = "#C0392B")

p_d3 <- q4030_scen |>
  filter(scenario %in% deck_scens_q40, cause == "All-cause",
         sex == "Female", year <= 2055) |>
  mutate(scenario_lab = recode(scenario, !!!deck_scen_labs_q40)) |>
  ggplot(aes(x = year, y = q4030_pct / 100,
             colour = scenario_lab, linetype = scenario_lab)) +
  geom_vline(xintercept = 2030, colour = "grey60",
             linetype = "dotted", linewidth = 0.5) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = deck_line_cols) +
  scale_linetype_manual(values = deck_line_types) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 0.1)) +
  labs(title    = "Female all-cause 40q30 under the WPP envelope",
       subtitle = "Interventions accelerate the ongoing WPP-driven mortality decline",
       x        = NULL, y = "40q30 (%)",
       colour   = NULL, linetype = NULL,
       caption  = paste0("Source: Indonesia NCD V1 | Female | All-cause 40q30 under WPP 2024 envelope",
                         " | Vertical line = 2030 (SDG reference year)")) +
  theme_deck

ggsave(file.path(OUT_FIGURES, "q4030_trajectory.png"),
       p_d3, width = 9, height = 5.5, dpi = 180)
message("  q4030_trajectory.png ✓")

# ── D4: sdg_aroc_updated.png ──────────────────────────────────────────────────
# Horizontal bar: AROC 2025-2030 with both benchmarks.
# SDG_CATCHUP_AROC and SDG_CLASSIC_AROC are computed at lines 169/172 above.
# Using those values directly here — do not redefine.

d4_data <- t3_sdg |>
  filter(year == 2030L, sex == "Female",
         scenario %in% deck_scenarios_averted |
           scenario %in% c("baseline", "statins_fast")) |>
  mutate(
    scenario_lab = factor(
      recode(scenario,
             all_fast     = "All CVD (fast)",
             bp_fast      = "BP control (fast)",
             all_slow     = "All CVD (slow)",
             statins_fast = "Statins (fast)",
             baseline     = "Baseline"),
      levels = c("All CVD (fast)", "BP control (fast)", "All CVD (slow)",
                 "Statins (fast)", "Baseline")
    )
  ) |>
  filter(!is.na(scenario_lab))

p_d4 <- d4_data |>
  ggplot(aes(x = aroc_pct_per_year, y = scenario_lab)) +
  geom_vline(xintercept = SDG_CLASSIC_AROC,
             colour = "grey55", linetype = "dashed", linewidth = 0.7) +
  geom_vline(xintercept = SDG_CATCHUP_AROC,
             colour = "#C0392B", linetype = "dashed", linewidth = 0.7) +
  geom_col(aes(fill = aroc_pct_per_year <= SDG_CLASSIC_AROC),
           width = 0.7) +
  geom_text(aes(label = paste0(round(aroc_pct_per_year, 1), "%/yr"),
                hjust = ifelse(aroc_pct_per_year < -1, 1.15, -0.1)),
            size = 3.5, fontface = "bold", colour = "grey20") +
  scale_fill_manual(values = c("TRUE" = "#27AE60", "FALSE" = "#E74C3C"),
                    guide   = "none") +
  scale_x_continuous(
    limits = c(min(SDG_CATCHUP_AROC - 1, min(d4_data$aroc_pct_per_year, na.rm=T) - 0.5), 0.5),
    labels = scales::label_number(suffix = "%/yr")) +
  labs(title    = "Annual rate of change in female 40q30, 2025\u20132030",
       subtitle = paste0("Red dashed = 2025 catch-up (", SDG_CATCHUP_AROC,
                         "%/yr) | Grey dashed = classic 2015-base SDG pace (",
                         SDG_CLASSIC_AROC, "%/yr)"),
       x        = "Annual rate of change in female all-cause 40q30 (%/year)",
       y        = NULL,
       caption  = "Source: Indonesia NCD V1 | Female all-cause 40q30 | Directional result; point estimates only") +
  theme_deck

ggsave(file.path(OUT_FIGURES, "sdg_aroc_updated.png"),
       p_d4, width = 10, height = 5, dpi = 180)
message("  sdg_aroc_updated.png ✓")

# ── DA: cause_fractions_projected.png ─────────────────────────────────────────
# GBD historical anchor years (solid) + V1 logit-slope projection (dashed)
# Female, representative age group nearest 60-64.

CVD_GBD_NAMES <- c("Ischemic heart disease", "Ischemic stroke",
                   "Intracerebral hemorrhage", "Hypertensive heart disease")
CVD_HIST_LABS <- c(
  "Ischemic heart disease"    = "IHD",
  "Ischemic stroke"            = "Isch. stroke",
  "Intracerebral hemorrhage"   = "ICH",
  "Hypertensive heart disease" = "HHD"
)
# Local colour vector keyed by label (not module ID) for the two new plots.
CVD_LABEL_COLS <- c(
  "IHD"         = "#8B2500",
  "Isch. stroke" = "#E74C3C",
  "ICH"          = "#F39C12",
  "HHD"          = "#E67E22"
)

frac_file <- here("data", "gbd", "gbd_frac_annual.csv")
if (file.exists(frac_file)) {
  frac_df <- read_csv(frac_file, show_col_types = FALSE)
  ref_age  <- frac_df |>
    filter(sex == "Female", cause %in% CVD_GBD_NAMES) |>
    distinct(age_mid) |>
    mutate(d = abs(age_mid - 62)) |>
    arrange(d) |>
    slice(1) |> pull(age_mid)
  
  p_da <- frac_df |>
    filter(sex == "Female", age_mid == ref_age, cause %in% CVD_GBD_NAMES) |>
    mutate(
      cause_lab = factor(recode(cause, !!!CVD_HIST_LABS),
                         levels = c("IHD", "Isch. stroke", "ICH", "HHD")),
      period    = if_else(year <= 2023, "GBD historical", "V1 projected")
    ) |>
    ggplot(aes(x = year, y = frac * 100, colour = cause_lab, linetype = period)) +
    geom_vline(xintercept = 2023.5, colour = "grey50", linetype = "dotted",
               linewidth = 0.5) +
    annotate("rect", xmin = 2024, xmax = Inf, ymin = -Inf, ymax = Inf,
             alpha = 0.04, fill = "steelblue") +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = CVD_LABEL_COLS) +
    scale_linetype_manual(
      values = c("GBD historical" = "solid", "V1 projected" = "dashed")) +
    scale_x_continuous(breaks = c(2000, 2010, 2023, 2050, 2075, 2100)) +
    labs(
      title    = paste0("Four-CVD cause fractions: historical and V1 projected, ",
                        "Female age ", ref_age, "\u201364"),
      subtitle = paste0("Solid = GBD 2023 anchor interpolated | ",
                        "Dashed = V1 logit-slope projection | ",
                        "Shaded = projected period"),
      x = NULL, y = "Cause fraction (%)",
      colour = "Cause", linetype = NULL,
      caption = paste0("Source: GBD 2023 historical | Indonesia NCD V1 logit-slope projection ",
                       "re-anchored annually to WPP all-cause mortality")
    ) +
    theme_deck
  ggsave(file.path(OUT_FIGURES, "cause_fractions_projected.png"),
         p_da, width = 10, height = 5, dpi = 180)
  message("  cause_fractions_projected.png \u2713")
} else {
  message("  cause_fractions_projected.png  SKIPPED (gbd_frac_annual.csv not found)")
}


# ── DA2/DA3: incidence context plots ─────────────────────────────────────────
# Historical + projected incidence, using the same visual convention as the
# deaths history bridge:
#   dotted = raw GBD historical incident-case counts, if available
#   solid  = GBD historical incidence rates applied to WPP historical population
#   dashed = V1 WPP-rebased projection
#
# These overwrite the existing deck filenames, so the Rnw does not need new
# figure paths:
#   incidence_cases_60plus_by_cause.png
#   incidence_rates_60plus_by_cause_sex.png
#
# Interpretation note:
#   V1 projection is not only projected mortality cause fractions. It combines
#   WPP population + WPP all-cause mortality envelope + projected mortality cause
#   shares + module-specific incidence / CF / transition inputs. These incidence
#   plots show the state-module incidence inputs in their historical WPP context.

INCIDENCE_AGE_RANGE <- c(60L, 100L)
# GBD 2023 historical anchors run through 2023. Some intermediate files may
# carry a 2024 row, but WPP historical package data often stop at 2023 while
# the V1 projection population starts at 2025. We therefore omit 2024 from the
# historical bridge unless it is explicitly part of PROJ_YEARS. This avoids
# treating a denominator gap as a model failure.
HISTORICAL_LAST_YEAR <- min(2023L, min(PROJ_YEARS, na.rm = TRUE) - 1L)
incidence_file <- here("data", "gbd", "gbd_incidence_annual_1yr.csv")
pop_file       <- file.path(BASELINE_DIR, "pop_df.rds")

# Helper: WPP 2024 package stores age as character in some versions.
.wpp_age_to_int <- function(x) {
  x_chr <- as.character(x)
  x_chr[x_chr %in% c("100+", "100 plus", "100_plus")] <- "100"
  suppressWarnings(as.integer(x_chr))
}

# Helper: convert incidence rates to per-person-year probabilities/rates.
# Most pipeline inputs are already per person. If a file arrives as per 100,000,
# the heuristic below prevents multiplying per-100k rates by population.
.normalise_incidence_rate <- function(x) {
  x <- as.numeric(x)
  max_x <- suppressWarnings(max(x, na.rm = TRUE))
  med_x <- suppressWarnings(stats::median(x[x > 0], na.rm = TRUE))
  if (!is.finite(max_x)) return(x)
  if (max_x > 1 || (is.finite(med_x) && med_x > 1)) x / 1e5 else x
}

# Helper: find an optional raw incident-case count column. We only use columns
# with count-like names; generic `val` is intentionally avoided because in IHME
# extracts it may be either a rate or a count depending on selected metric.
.find_inc_count_col <- function(nms) {
  candidates <- c("incident_cases", "cases", "case_count", "incidence_count",
                  "incidence_n", "count", "n")
  hit <- intersect(candidates, nms)
  if (length(hit) == 0) NA_character_ else hit[1]
}

# Helper: ensure incidence table has V1 module IDs.
.standardise_incidence_modules <- function(df) {
  out <- df
  if (!"module" %in% names(out)) {
    if (!"cause" %in% names(out)) {
      stop("Incidence file must contain either `module` or `cause`.", call. = FALSE)
    }
    out <- out |>
      mutate(
        module = case_when(
          cause %in% c("ihd", "IHD", "Ischemic heart disease", "Ischaemic heart disease") ~ "ihd",
          cause %in% c("ischemic_stroke", "ischaemic_stroke", "istroke", "Ischemic stroke", "Ischaemic stroke") ~ "ischemic_stroke",
          cause %in% c("ich", "hstroke", "Intracerebral hemorrhage", "Intracerebral haemorrhage") ~ "ich",
          cause %in% c("hhd", "Hypertensive heart disease") ~ "hhd",
          TRUE ~ as.character(cause)
        )
      )
  }
  out |>
    mutate(
      module = recode(as.character(module),
                      "istroke" = "ischemic_stroke",
                      "ischaemic_stroke" = "ischemic_stroke",
                      "hstroke" = "ich",
                      .default = as.character(module)),
      year = as.integer(year),
      age  = as.integer(age),
      sex  = as.character(sex)
    )
}

if (file.exists(incidence_file) && file.exists(pop_file)) {
  inc_raw <- readr::read_csv(incidence_file, show_col_types = FALSE) |>
    .standardise_incidence_modules()
  
  inc_rate_col <- intersect(c("inc_rate", "incidence_rate", "rate", "IR", "ir"),
                            names(inc_raw))[1]
  if (is.na(inc_rate_col)) {
    stop(
      "Cannot build incidence plots: no incidence-rate column found in ",
      incidence_file, " (checked inc_rate, incidence_rate, rate, IR, ir).",
      call. = FALSE
    )
  }
  
  inc_count_col <- .find_inc_count_col(names(inc_raw))
  
  inc_1yr <- inc_raw |>
    mutate(
      inc_rate_raw  = as.numeric(.data[[inc_rate_col]]),
      inc_rate_prob = .normalise_incidence_rate(inc_rate_raw),
      raw_cases     = if (!is.na(inc_count_col)) as.numeric(.data[[inc_count_col]]) else NA_real_
    ) |>
    filter(
      module %in% setdiff(names(CVD_CAUSE_LABS), "hhd"),
      !is.na(inc_rate_prob),
      age >= INCIDENCE_AGE_RANGE[1],
      age <= INCIDENCE_AGE_RANGE[2],
      year <= HISTORICAL_LAST_YEAR | year %in% PROJ_YEARS
    ) |>
    mutate(
      cause_lab = factor(recode(module, !!!CVD_CAUSE_LABS),
                         levels = c("IHD", "Isch. stroke", "ICH"))
    )
  
  # Projection population: baseline pop_df.rds, normally 2025+.
  pop_proj <- readRDS(pop_file) |>
    mutate(
      year = as.integer(year),
      age  = as.integer(age),
      sex  = as.character(sex),
      pop  = as.numeric(pop)
    ) |>
    select(year, sex, age, pop)
  
  # Historical population: WPP historical mx/pop package when available.
  hist_years_inc <- inc_1yr |>
    filter(year <= HISTORICAL_LAST_YEAR) |>
    distinct(year) |>
    arrange(year) |>
    pull(year) |>
    as.integer()
  
  if (length(hist_years_inc) > 0 && requireNamespace("wpp2024", quietly = TRUE)) {
    wpp_env_inc <- new.env(parent = baseenv())
    data("popAge1dt", package = "wpp2024", envir = wpp_env_inc)
    
    pop_hist <- wpp_env_inc$popAge1dt |>
      filter(name == "Indonesia", year %in% hist_years_inc, !is.na(age)) |>
      transmute(
        year = as.integer(year),
        age  = .wpp_age_to_int(age),
        Female = as.numeric(popF) * 1e3,
        Male   = as.numeric(popM) * 1e3
      ) |>
      filter(!is.na(age),
             age >= INCIDENCE_AGE_RANGE[1],
             age <= INCIDENCE_AGE_RANGE[2]) |>
      pivot_longer(c(Female, Male), names_to = "sex", values_to = "pop")
  } else {
    if (length(hist_years_inc) > 0) {
      warning(
        "wpp2024 is not installed, so WPP-adjusted historical incidence anchors ",
        "are skipped. Install wpp2024 to show historical incidence context.",
        call. = FALSE
      )
    }
    pop_hist <- tibble(year = integer(), sex = character(), age = integer(), pop = numeric())
  }
  
  # Combine historical and projected WPP population.  The incidence input can
  # include bridge years (especially 2024) that are not always present in either
  # WPP historical tables or the V1 projection population file.  Rather than
  # failing at this harmless join gap, build the exact population keys needed by
  # the incidence table and linearly interpolate / nearest-fill WPP population
  # within each sex-age stratum.  This preserves the intended convention:
  # historical incidence rates are put on the WPP population denominator, while
  # projected incidence rates use the V1 WPP-rebased population.
  pop_available_inc <- bind_rows(pop_hist, pop_proj) |>
    mutate(
      year = as.integer(year),
      age  = as.integer(age),
      sex  = as.character(sex),
      pop  = as.numeric(pop)
    ) |>
    filter(!is.na(year), !is.na(age), !is.na(sex), is.finite(pop), pop >= 0) |>
    distinct(year, sex, age, .keep_all = TRUE)
  
  # Downstream incidence diagnostics use this all-years WPP population object.
  # It includes historical WPP population when available and V1 projection
  # population for projection years, with one row per year-sex-age.
  pop_all_inc <- pop_available_inc
  
  # If historical WPP population was unavailable, skip historical incidence rows
  # rather than extrapolating 2025+ population backward across decades.
  if (length(hist_years_inc) > 0 && nrow(pop_hist) == 0) {
    warning(
      "Historical WPP population was unavailable; incidence history will be skipped, ",
      "but projected incidence figures will still be written.",
      call. = FALSE
    )
    inc_1yr <- inc_1yr |>
      filter(year %in% PROJ_YEARS)
    hist_years_inc <- integer()
  }
  
  needed_pop_keys_inc <- inc_1yr |>
    filter(year %in% c(hist_years_inc, PROJ_YEARS)) |>
    distinct(year, sex, age) |>
    arrange(sex, age, year)
  
  interpolate_pop_one_stratum <- function(keys, key) {
    avail <- pop_available_inc |>
      filter(sex == key$sex[[1]], age == key$age[[1]]) |>
      arrange(year)
    
    if (nrow(avail) == 0) {
      keys$pop <- NA_real_
      return(keys)
    }
    
    avail <- avail |>
      filter(is.finite(pop)) |>
      distinct(year, .keep_all = TRUE)
    
    if (nrow(avail) == 1) {
      keys$pop <- avail$pop[[1]]
      return(keys)
    }
    
    keys$pop <- as.numeric(stats::approx(
      x = avail$year,
      y = avail$pop,
      xout = keys$year,
      rule = 2,
      ties = "ordered"
    )$y)
    keys
  }
  
  pop_needed_inc <- needed_pop_keys_inc |>
    group_by(sex, age) |>
    group_modify(interpolate_pop_one_stratum) |>
    ungroup()
  
  inc_joined <- inc_1yr |>
    left_join(pop_needed_inc, by = c("year", "sex", "age"))
  
  missing_pop_inc <- inc_joined |>
    filter(is.na(pop), year %in% PROJ_YEARS | year %in% hist_years_inc) |>
    distinct(module, year, sex, age)
  
  missing_pop_proj_inc <- missing_pop_inc |>
    filter(year %in% PROJ_YEARS)
  missing_pop_hist_inc <- missing_pop_inc |>
    filter(year %in% hist_years_inc)
  
  if (nrow(missing_pop_proj_inc) > 0) {
    stop(
      "Cannot build projected incidence plots: missing V1/WPP population rows. First rows:\n",
      paste(utils::capture.output(print(utils::head(missing_pop_proj_inc, 20))), collapse = "\n"),
      call. = FALSE
    )
  }
  
  if (nrow(missing_pop_hist_inc) > 0) {
    warning(
      "Dropping historical incidence rows with no WPP population match. First rows:\n",
      paste(utils::capture.output(print(utils::head(missing_pop_hist_inc, 20))), collapse = "\n"),
      call. = FALSE
    )
    inc_joined <- inc_joined |>
      filter(!(is.na(pop) & year %in% hist_years_inc))
  }
  
  inc_joined <- inc_joined |>
    mutate(
      incident_cases_wpp = inc_rate_prob * pop,
      series = case_when(
        year <= HISTORICAL_LAST_YEAR       ~ "WPP-adjusted historical",
        year %in% PROJ_YEARS                ~ "V1 WPP projection",
        TRUE                                ~ NA_character_
      )
    ) |>
    filter(!is.na(series))
  
  # Raw historical cases are optional and included only if the input file
  # carries a clearly count-like column.
  incidence_raw_cases_60plus <- inc_1yr |>
    filter(year <= HISTORICAL_LAST_YEAR, !is.na(raw_cases)) |>
    group_by(module, cause_lab, year) |>
    summarise(incident_cases = sum(raw_cases, na.rm = TRUE), .groups = "drop") |>
    mutate(series = "Raw GBD counts")
  
  incidence_cases_60plus <- bind_rows(
    inc_joined |>
      group_by(module, cause_lab, year, series) |>
      summarise(
        population_60plus = sum(pop, na.rm = TRUE),
        incident_cases    = sum(incident_cases_wpp, na.rm = TRUE),
        inc_rate_60plus   = incident_cases / pmax(population_60plus, 1),
        inc_rate_per100k  = inc_rate_60plus * 1e5,
        .groups = "drop"
      ),
    incidence_raw_cases_60plus |>
      mutate(
        population_60plus = NA_real_,
        inc_rate_60plus   = NA_real_,
        inc_rate_per100k  = NA_real_
      ) |>
      select(module, cause_lab, year, series, population_60plus,
             incident_cases, inc_rate_60plus, inc_rate_per100k)
  ) |>
    mutate(
      series = factor(series,
                      levels = c("Raw GBD counts", "WPP-adjusted historical", "V1 WPP projection")),
      cause_lab = factor(as.character(cause_lab), levels = c("IHD", "Isch. stroke", "ICH"))
    ) |>
    arrange(cause_lab, series, year)
  
  readr::write_csv(
    incidence_cases_60plus,
    file.path(OUT_TABLES, "t0e_incidence_cases_60plus_by_cause.csv")
  )
  
  p_inc_cases <- incidence_cases_60plus |>
    ggplot2::ggplot(ggplot2::aes(x = year, y = incident_cases / 1e3,
                                 colour = cause_lab, linetype = series)) +
    ggplot2::geom_vline(xintercept = min(PROJ_YEARS, na.rm = TRUE) - 0.5,
                        colour = "grey50", linetype = "dotted", linewidth = 0.5) +
    ggplot2::geom_line(linewidth = 0.9, na.rm = TRUE) +
    ggplot2::geom_point(
      data = incidence_cases_60plus |>
        filter(series != "V1 WPP projection" | year == min(PROJ_YEARS, na.rm = TRUE)),
      ggplot2::aes(shape = series),
      size = 2.0, stroke = 0.8, fill = "white", na.rm = TRUE
    ) +
    ggplot2::scale_colour_manual(values = CVD_LABEL_COLS, drop = FALSE) +
    ggplot2::scale_linetype_manual(values = c(
      "Raw GBD counts"          = "dotted",
      "WPP-adjusted historical" = "solid",
      "V1 WPP projection"       = "dashed"
    ), drop = FALSE) +
    ggplot2::scale_shape_manual(values = c(
      "Raw GBD counts"          = 16,
      "WPP-adjusted historical" = 21,
      "V1 WPP projection"       = 21
    ), drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = c(2000, 2010, 2023, 2025, 2050, 2075, 2100)) +
    ggplot2::scale_y_continuous(labels = scales::label_number(suffix = "k"),
                                expand = expansion(mult = c(0, 0.05))) +
    ggplot2::labs(
      title    = "Annual incident cases ages 60+: raw GBD, WPP-adjusted history, and V1 baseline",
      subtitle = paste0(
        "Dotted = raw GBD counts | Solid = GBD incidence rates on WPP historical population | ",
        "Dashed = V1 WPP-rebased projection"
      ),
      x = NULL,
      y = "Annual incident cases",
      colour = "Cause",
      linetype = "Series",
      shape = "Series",
      caption = paste0(
        "Source: GBD 2023 historical incidence; WPP 2024 population; ",
        "Indonesia NCD V1 baseline incidence projection | HHD excluded because it is a direct-mortality module in V1"
      )
    ) +
    theme_deck
  
  ggplot2::ggsave(
    file.path(OUT_FIGURES, "incidence_cases_60plus_by_cause.png"),
    p_inc_cases, width = 10, height = 5.3, dpi = 180
  )
  message("  incidence_cases_60plus_by_cause.png ✓")
  
  incidence_rates_60plus <- inc_joined |>
    group_by(module, cause_lab, year, sex, series) |>
    summarise(
      population_60plus = sum(pop, na.rm = TRUE),
      incident_cases    = sum(incident_cases_wpp, na.rm = TRUE),
      inc_rate_60plus   = incident_cases / pmax(population_60plus, 1),
      inc_rate_per100k  = inc_rate_60plus * 1e5,
      .groups = "drop"
    ) |>
    mutate(
      series = factor(series,
                      levels = c("WPP-adjusted historical", "V1 WPP projection")),
      cause_lab = factor(as.character(cause_lab), levels = c("IHD", "Isch. stroke", "ICH"))
    ) |>
    arrange(cause_lab, sex, series, year)
  
  readr::write_csv(
    incidence_rates_60plus,
    file.path(OUT_TABLES, "t0d_incidence_rates_60plus_by_cause_sex.csv")
  )
  
  p_inc_rate <- incidence_rates_60plus |>
    ggplot2::ggplot(ggplot2::aes(x = year, y = inc_rate_per100k,
                                 colour = sex, linetype = series)) +
    ggplot2::geom_vline(xintercept = min(PROJ_YEARS, na.rm = TRUE) - 0.5,
                        colour = "grey50", linetype = "dotted", linewidth = 0.5) +
    ggplot2::geom_line(linewidth = 0.85) +
    ggplot2::facet_wrap(~ cause_lab, scales = "free_y", ncol = 2) +
    ggplot2::scale_x_continuous(breaks = c(2000, 2010, 2023, 2025, 2050, 2075, 2100)) +
    ggplot2::scale_y_continuous(labels = scales::label_number()) +
    ggplot2::scale_colour_manual(values = c("Female" = "#C0392B", "Male" = "#1F6BAE")) +
    ggplot2::scale_linetype_manual(values = c("WPP-adjusted historical" = "solid",
                                              "V1 WPP projection" = "dashed")) +
    ggplot2::labs(
      title    = "Incidence rates ages 60+: WPP-adjusted history and V1 baseline",
      subtitle = paste0(
        "Population-weighted across single-year ages ", INCIDENCE_AGE_RANGE[1],
        "–", INCIDENCE_AGE_RANGE[2],
        " | Solid = historical GBD rates on WPP population | Dashed = V1 projection"
      ),
      x = NULL,
      y = "Incident cases per 100,000 population ages 60+",
      colour = "Sex",
      linetype = "Series",
      caption = paste0(
        "Source: GBD 2023 incidence rates; WPP 2024 population; ",
        "Indonesia NCD V1 baseline incidence projection | HHD excluded because it is a direct-mortality module in V1"
      )
    ) +
    theme_deck
  
  ggplot2::ggsave(
    file.path(OUT_FIGURES, "incidence_rates_60plus_by_cause_sex.png"),
    p_inc_rate, width = 10, height = 5.8, dpi = 180
  )
  message("  incidence_rates_60plus_by_cause_sex.png ✓")
  
  # Optional diagnostic: age-standardised rates using 2025 WPP 60+ age weights.
  # This helps distinguish true rate-projection behaviour from age-composition
  # shifts within the open 60+ population.
  #
  # `pop_needed_inc` is the incidence-specific WPP population bridge built above
  # from historical WPP population plus the V1 projection population.  Earlier
  # versions of this diagnostic used `pop_all_inc`; that object is not created
  # in the current incidence block, so use the bridge object directly.
  std_weights <- pop_needed_inc |>
    filter(year == min(PROJ_YEARS, na.rm = TRUE),
           age >= INCIDENCE_AGE_RANGE[1], age <= INCIDENCE_AGE_RANGE[2]) |>
    group_by(sex) |>
    mutate(weight = pop / pmax(sum(pop, na.rm = TRUE), 1)) |>
    ungroup() |>
    select(sex, age, weight)
  
  if (nrow(std_weights) > 0) {
    incidence_rates_60plus_std <- inc_joined |>
      left_join(std_weights, by = c("sex", "age")) |>
      filter(!is.na(weight)) |>
      group_by(module, cause_lab, year, sex, series) |>
      summarise(inc_rate_per100k = sum(inc_rate_prob * weight, na.rm = TRUE) * 1e5,
                .groups = "drop") |>
      mutate(
        series = factor(series,
                        levels = c("WPP-adjusted historical", "V1 WPP projection")),
        cause_lab = factor(as.character(cause_lab), levels = c("IHD", "Isch. stroke", "ICH"))
      )
    
    readr::write_csv(
      incidence_rates_60plus_std,
      file.path(OUT_TABLES, "t0f_incidence_rates_60plus_age_std_by_cause_sex.csv")
    )
    
    p_inc_rate_std <- incidence_rates_60plus_std |>
      ggplot2::ggplot(ggplot2::aes(x = year, y = inc_rate_per100k,
                                   colour = sex, linetype = series)) +
      ggplot2::geom_vline(xintercept = min(PROJ_YEARS, na.rm = TRUE) - 0.5,
                          colour = "grey50", linetype = "dotted", linewidth = 0.5) +
      ggplot2::geom_line(linewidth = 0.85) +
      ggplot2::facet_wrap(~ cause_lab, scales = "free_y", ncol = 2) +
      ggplot2::scale_x_continuous(breaks = c(2000, 2010, 2023, 2025, 2050, 2075, 2100)) +
      ggplot2::scale_y_continuous(labels = scales::label_number()) +
      ggplot2::scale_colour_manual(values = c("Female" = "#C0392B", "Male" = "#1F6BAE")) +
      ggplot2::scale_linetype_manual(values = c("WPP-adjusted historical" = "solid",
                                                "V1 WPP projection" = "dashed")) +
      ggplot2::labs(
        title    = "Age-standardised incidence rates ages 60+: historical context and V1 baseline",
        subtitle = "Standardised to the 2025 WPP 60+ age distribution within sex",
        x = NULL,
        y = "Age-standardised incident cases per 100,000 ages 60+",
        colour = "Sex",
        linetype = "Series",
        caption = "Source: GBD incidence rates; WPP 2024 population; Indonesia NCD V1 baseline incidence projection"
      ) +
      theme_deck
    
    ggplot2::ggsave(
      file.path(OUT_FIGURES, "incidence_rates_60plus_age_std_by_cause_sex.png"),
      p_inc_rate_std, width = 10, height = 5.8, dpi = 180
    )
    message("  incidence_rates_60plus_age_std_by_cause_sex.png ✓")
  }
  
  # Optional diagnostic: ICH male single-age rates, grouped for readability.
  ich_male_diag <- inc_joined |>
    filter(module == "ich", sex == "Male") |>
    mutate(age_group = cut(age,
                           breaks = c(60, 65, 70, 75, 80, 85, 90, 95, 101),
                           right = FALSE,
                           labels = c("60-64", "65-69", "70-74", "75-79",
                                      "80-84", "85-89", "90-94", "95+"))) |>
    group_by(year, age_group, series) |>
    summarise(inc_rate_per100k = mean(inc_rate_prob, na.rm = TRUE) * 1e5,
              .groups = "drop") |>
    filter(!is.na(age_group))
  
  if (nrow(ich_male_diag) > 0) {
    p_ich_diag <- ich_male_diag |>
      ggplot2::ggplot(ggplot2::aes(x = year, y = inc_rate_per100k,
                                   colour = age_group, linetype = series)) +
      ggplot2::geom_vline(xintercept = min(PROJ_YEARS, na.rm = TRUE) - 0.5,
                          colour = "grey50", linetype = "dotted", linewidth = 0.5) +
      ggplot2::geom_line(linewidth = 0.75) +
      ggplot2::scale_x_continuous(breaks = c(2000, 2010, 2023, 2025, 2050, 2075, 2100)) +
      ggplot2::labs(
        title    = "ICH male incidence diagnostic: age-specific rates",
        subtitle = "Used to review whether the aggregate male ICH decline is age-composition or input-rate driven",
        x = NULL,
        y = "Incident cases per 100,000",
        colour = "Age",
        linetype = "Series",
        caption = "Source: GBD incidence rates; WPP 2024 population; Indonesia NCD V1 baseline incidence projection"
      ) +
      theme_deck
    
    ggplot2::ggsave(
      file.path(OUT_FIGURES, "ich_male_age_specific_incidence_diagnostic.png"),
      p_ich_diag, width = 10, height = 5.8, dpi = 180
    )
    message("  ich_male_age_specific_incidence_diagnostic.png ✓")
  }
} else {
  message("  incidence plots SKIPPED (requires gbd_incidence_annual_1yr.csv and baseline/pop_df.rds)")
}

# ── DB: deaths_historical_projected.png ──────────────────────────────────────
# Raw GBD historical counts + WPP-adjusted historical anchors + V1 WPP baseline.
#
# Why both historical series are shown:
#   Raw GBD counts are the observed historical anchors.
#   WPP-adjusted historical anchors apply each GBD cause fraction to WPP mx/pop
#   for the same year, putting the historical cause composition onto the same
#   denominator/envelope as the V1 projection. This makes the 2025 WPP-projected
#   baseline visually interpretable instead of appearing as an unexplained jump.

hist_file     <- here("data", "gbd", "gbd_cause_deaths_n.csv")
frac_1yr_file <- here("data", "gbd", "gbd_frac_annual_1yr.csv")

if (file.exists(hist_file) && file.exists(frac_1yr_file)) {
  hist_raw <- read_csv(hist_file, show_col_types = FALSE)
  
  # Column name for death counts varies by pipeline version; detect defensively.
  deaths_col <- intersect(c("count", "deaths_n", "deaths", "val", "n"),
                          names(hist_raw))[1]
  if (is.na(deaths_col)) {
    stop(
      "deaths_historical_projected: cannot find death count column in ",
      "gbd_cause_deaths_n.csv (checked: count, deaths_n, deaths, val, n). ",
      "Available columns: ", paste(names(hist_raw), collapse = ", "),
      call. = FALSE
    )
  }
  
  hist_years <- hist_raw |>
    filter(cause %in% CVD_GBD_NAMES) |>
    distinct(year) |>
    arrange(year) |>
    pull(year) |>
    as.integer()
  
  cause_levels <- c("IHD", "Isch. stroke", "ICH", "HHD")
  
  # 1) Raw GBD cause death counts.
  hist_deaths_raw <- hist_raw |>
    filter(cause %in% CVD_GBD_NAMES) |>
    mutate(.deaths = .data[[deaths_col]]) |>
    group_by(cause, year) |>
    summarise(deaths = sum(.deaths, na.rm = TRUE), .groups = "drop") |>
    mutate(
      cause_lab = factor(recode(cause, !!!CVD_HIST_LABS), levels = cause_levels),
      series    = "Raw GBD counts"
    ) |>
    select(year, deaths, cause_lab, series)
  
  # 2) WPP-adjusted historical anchors.
  #    Requires WPP historical mx/pop from the wpp2024 package. The active model
  #    uses WPP 2025+ outputs, but for this historical bridge figure we need
  #    WPP historical years too.
  if (!requireNamespace("wpp2024", quietly = TRUE)) {
    warning(
      "wpp2024 is not installed, so WPP-adjusted historical anchors are skipped ",
      "in deaths_historical_projected.png.",
      call. = FALSE
    )
    
    hist_deaths_wpp <- tibble(
      year = integer(), deaths = numeric(),
      cause_lab = factor(character(), levels = cause_levels),
      series = character()
    )
  } else {
    wpp_env <- new.env(parent = baseenv())
    data("popAge1dt", package = "wpp2024", envir = wpp_env)
    data("mx1dt",     package = "wpp2024", envir = wpp_env)
    
    .wpp_age_to_int <- function(x) {
      x_chr <- as.character(x)
      x_chr[x_chr %in% c("100+", "100 plus", "100_plus")] <- "100"
      suppressWarnings(as.integer(x_chr))
    }
    
    wpp_pop_hist <- wpp_env$popAge1dt |>
      filter(name == "Indonesia", year %in% hist_years, !is.na(age)) |>
      transmute(
        year = as.integer(year),
        age  = .wpp_age_to_int(age),
        Female = as.numeric(popF) * 1e3,
        Male   = as.numeric(popM) * 1e3
      ) |>
      filter(!is.na(age), age >= 0, age <= 100) |>
      pivot_longer(c(Female, Male), names_to = "sex", values_to = "pop")
    
    wpp_mx_hist <- wpp_env$mx1dt |>
      filter(name == "Indonesia", year %in% hist_years, !is.na(age)) |>
      transmute(
        year = as.integer(year),
        age  = .wpp_age_to_int(age),
        Female = as.numeric(mxF),
        Male   = as.numeric(mxM)
      ) |>
      filter(!is.na(age), age >= 0, age <= 100) |>
      pivot_longer(c(Female, Male), names_to = "sex", values_to = "mx_wpp")
    
    wpp_hist <- wpp_mx_hist |>
      left_join(wpp_pop_hist, by = c("year", "sex", "age"))
    
    missing_wpp_hist <- wpp_hist |>
      filter(is.na(mx_wpp) | is.na(pop)) |>
      distinct(year, sex, age)
    
    if (nrow(missing_wpp_hist) > 0) {
      stop(
        "deaths_historical_projected: missing WPP historical mx/pop rows. First rows:\n",
        paste(utils::capture.output(print(utils::head(missing_wpp_hist, 20))),
              collapse = "\n"),
        call. = FALSE
      )
    }
    
    frac_hist_1yr <- read_csv(frac_1yr_file, show_col_types = FALSE) |>
      filter(cause %in% CVD_GBD_NAMES, year %in% hist_years) |>
      mutate(
        year = as.integer(year),
        age  = as.integer(age),
        cause_lab = recode(cause, !!!CVD_HIST_LABS)
      ) |>
      select(cause_lab, year, sex, age, frac)
    
    # Expand to the full WPP age-sex spine. Missing rows outside GBD support
    # become structural zeroes for this bridge plot.
    hist_deaths_wpp <- tidyr::expand_grid(
      cause_lab = cause_levels,
      year      = hist_years,
      sex       = c("Female", "Male"),
      age       = 0:100
    ) |>
      left_join(frac_hist_1yr, by = c("cause_lab", "year", "sex", "age")) |>
      mutate(frac = tidyr::replace_na(frac, 0)) |>
      left_join(wpp_hist, by = c("year", "sex", "age")) |>
      group_by(cause_lab, year) |>
      summarise(deaths = sum(mx_wpp * frac * pop, na.rm = TRUE), .groups = "drop") |>
      mutate(
        cause_lab = factor(cause_lab, levels = cause_levels),
        series    = "WPP-adjusted historical"
      ) |>
      select(year, deaths, cause_lab, series)
  }
  
  # 3) V1 projected baseline on the WPP denominator/envelope.
  proj_deaths <- all_scenario_mx |>
    filter(scenario == "baseline", module %in% names(CVD_CAUSE_LABS)) |>
    group_by(module, year) |>
    summarise(deaths = sum(deaths_cause, na.rm = TRUE), .groups = "drop") |>
    mutate(
      cause_lab = factor(recode(module, !!!CVD_CAUSE_LABS), levels = cause_levels),
      series    = "V1 WPP projection"
    ) |>
    select(year, deaths, cause_lab, series)
  
  deaths_comb <- bind_rows(hist_deaths_raw, hist_deaths_wpp, proj_deaths) |>
    mutate(
      series = factor(
        series,
        levels = c("Raw GBD counts", "WPP-adjusted historical", "V1 WPP projection")
      )
    )
  
  # Output a companion table so the raw-vs-WPP rebasing is auditable.
  deaths_hist_compare <- deaths_comb |>
    filter(year <= 2025) |>
    mutate(cause_lab = as.character(cause_lab)) |>
    pivot_wider(names_from = series, values_from = deaths) |>
    arrange(factor(cause_lab, levels = cause_levels), year)
  
  write_csv(
    deaths_hist_compare,
    file.path(OUT_TABLES, "t0_cvd_raw_gbd_vs_wpp_adjusted_history.csv")
  )
  message("  t0_cvd_raw_gbd_vs_wpp_adjusted_history.csv ✓")
  
  p_db <- deaths_comb |>
    ggplot(aes(x = year, y = deaths / 1e3,
               colour = cause_lab, linetype = series)) +
    geom_vline(xintercept = 2024.5, colour = "grey50", linetype = "dotted",
               linewidth = 0.5) +
    geom_line(linewidth = 0.85) +
    geom_point(
      data = deaths_comb |>
        filter(series != "V1 WPP projection" | year == min(PROJ_YEARS)),
      aes(shape = series),
      size = 2.2,
      stroke = 0.8,
      fill = "white"
    ) +
    scale_colour_manual(values = CVD_LABEL_COLS) +
    scale_linetype_manual(values = c(
      "Raw GBD counts"          = "dotted",
      "WPP-adjusted historical" = "solid",
      "V1 WPP projection"       = "dashed"
    )) +
    scale_shape_manual(values = c(
      "Raw GBD counts"          = 16,
      "WPP-adjusted historical" = 21,
      "V1 WPP projection"       = 21
    )) +
    scale_x_continuous(breaks = c(2000, 2010, 2023, 2050, 2075, 2100)) +
    labs(
      title    = "Four-CVD annual deaths: raw GBD, WPP-adjusted history, and V1 baseline",
      subtitle = paste0(
        "Dotted = raw GBD counts | Solid = GBD cause fractions on WPP historical mx/pop | ",
        "Dashed = V1 WPP-rebased projection"
      ),
      x = NULL,
      y = "Annual deaths (thousands)",
      colour = "Cause",
      linetype = "Series",
      shape = "Series",
      caption = paste0(
        "Source: GBD 2023 historical counts and fractions; WPP 2024 mx/pop; ",
        "Indonesia NCD V1 baseline projection"
      )
    ) +
    theme_deck
  
  ggsave(file.path(OUT_FIGURES, "deaths_historical_projected.png"),
         p_db, width = 10, height = 5, dpi = 180)
  message("  deaths_historical_projected.png ✓")
} else {
  message(
    "  deaths_historical_projected.png  SKIPPED (requires gbd_cause_deaths_n.csv ",
    "and gbd_frac_annual_1yr.csv)"
  )
}

# ── DC: deaths_60plus_comparison.png ─────────────────────────────────────────
# Annual deaths ages 60+, baseline vs all_fast, faceted by CVD cause.
p_dc <- all_scenario_mx |>
  filter(scenario %in% c("baseline", "all_fast"),
         module %in% names(CVD_CAUSE_LABS),
         age >= 60, year <= 2060) |>
  group_by(scenario, module, year) |>
  summarise(deaths = sum(deaths_cause, na.rm = TRUE), .groups = "drop") |>
  mutate(
    cause_lab    = factor(recode(module, !!!CVD_CAUSE_LABS),
                          levels = c("IHD", "Isch. stroke", "ICH", "HHD")),
    scenario_lab = recode(scenario,
                          "baseline" = "Baseline",
                          "all_fast" = "All CVD (fast)")
  ) |>
  ggplot(aes(x = year, y = deaths / 1e3,
             colour = scenario_lab, linetype = scenario_lab)) +
  geom_line(linewidth = 0.85) +
  facet_wrap(~ cause_lab, scales = "free_y", ncol = 2) +
  scale_colour_manual(
    values = c("Baseline" = "#2C3E50", "All CVD (fast)" = "#1F6BAE")) +
  scale_linetype_manual(
    values = c("Baseline" = "solid", "All CVD (fast)" = "dashed")) +
  labs(
    title    = "Annual CVD deaths ages 60+: baseline vs all-fast CVD package, 2025\u20132060",
    subtitle = "Faceted by cause | Both sexes combined | Ages 60\u2013100",
    x = NULL, y = "Annual deaths (thousands)",
    colour = NULL, linetype = NULL,
    caption = "Source: Indonesia NCD V1 | Ages 60+ | Both sexes | Point estimates only"
  ) +
  theme_deck +
  theme(strip.text = element_text(face = "bold", size = 9))
ggsave(file.path(OUT_FIGURES, "deaths_60plus_comparison.png"),
       p_dc, width = 10, height = 6, dpi = 180)
message("  deaths_60plus_comparison.png \u2713")


# ── DD: deck_numbers.tex for the Rnw deck ────────────────────────────────────
# The collaborator Rnw uses LaTeX macros for headline values and scenario
# assumptions. Write them here so the deck is regenerated from the same output
# tables and figures. Assumption values are read from an optional scenario
# assumptions file when available; otherwise the deck renders placeholders rather
# than hard-coding stale values.

DECK_NUMBERS_FILE <- here("deck_numbers.tex")

fmt_num <- function(x, digits = 1) {
  if (length(x) == 0 || is.na(x) || !is.finite(x)) return("??")
  format(round(as.numeric(x), digits), nsmall = digits, big.mark = ",", trim = TRUE)
}
fmt_millions <- function(x, digits = 2) {
  if (length(x) == 0 || is.na(x) || !is.finite(x)) return("??")
  paste0(format(round(as.numeric(x) / 1e6, digits), nsmall = digits, big.mark = ",", trim = TRUE), "M")
}
fmt_pct <- function(x, digits = 1) {
  if (length(x) == 0 || is.na(x) || !is.finite(x)) return("??")
  paste0(format(round(as.numeric(x), digits), nsmall = digits, trim = TRUE), "\\%")
}
tex_text <- function(x) {
  if (length(x) == 0 || is.na(x)) return("??")
  x <- as.character(x)
  x <- gsub("_", "\\\\_", x, fixed = TRUE)
  x <- gsub("&", "\\\\&", x, fixed = TRUE)
  x <- gsub("(?<!\\\\)%", "\\\\%", x, perl = TRUE)
  x
}
tex_cmd <- function(name, value) paste0("\\newcommand{\\", name, "}{", as.character(value), "}")
tex_cmd_text <- function(name, value) tex_cmd(name, tex_text(value))

# Headline result numbers from deaths_averted/cum_averted.
headline_window <- c(2025L, 2050L)
headline_scenario <- "all_fast"
headline_modules <- names(CVD_CAUSE_LABS)

headline_rows <- deaths_averted |>
  filter(scenario == headline_scenario,
         year >= headline_window[1], year <= headline_window[2],
         module %in% headline_modules)

headline_total_cvd <- sum(headline_rows$deaths_base, na.rm = TRUE)
headline_averted   <- sum(headline_rows$averted_display, na.rm = TRUE)
headline_pct       <- if (headline_total_cvd > 0) 100 * headline_averted / headline_total_cvd else NA_real_
headline_stroke_ich <- headline_rows |>
  filter(module %in% c("ischemic_stroke", "ich")) |>
  summarise(x = sum(deaths_base, na.rm = TRUE)) |>
  pull(x)

headline_aroc <- t3_sdg |>
  filter(scenario == headline_scenario, sex == "Female", year == 2030L) |>
  summarise(x = dplyr::first(aroc_pct_per_year)) |>
  pull(x)

q_base_2030 <- q4030_scen |>
  filter(scenario == "baseline", sex == "Female", cause == "All-cause", year == 2030L) |>
  summarise(x = dplyr::first(q4030_pct)) |>
  pull(x)
q_all_2030 <- q4030_scen |>
  filter(scenario == headline_scenario, sex == "Female", cause == "All-cause", year == 2030L) |>
  summarise(x = dplyr::first(q4030_pct)) |>
  pull(x)
headline_q_reduction <- if (length(q_base_2030) > 0 && length(q_all_2030) > 0 &&
                            isTRUE(is.finite(q_base_2030)) && isTRUE(is.finite(q_all_2030)) &&
                            isTRUE(q_base_2030 > 0)) {
  100 * (q_base_2030 - q_all_2030) / q_base_2030
} else NA_real_

# Optional scenario assumptions file. This intentionally supports several
# simple schemas so the plotting script does not need to know how the scenario
# builder stores metadata.
assumption_files <- c(
  file.path(SCEN_DIR, "scenario_assumptions.csv"),
  file.path(SCEN_DIR, "scenario_metadata.csv"),
  file.path(OUT_TABLES, "scenario_assumptions.csv"),
  here("data", "model", "scenario_assumptions.csv")
)
assumption_file <- assumption_files[file.exists(assumption_files)][1]
scenario_assumptions <- NULL
if (!is.na(assumption_file)) {
  scenario_assumptions <- suppressMessages(readr::read_csv(assumption_file, show_col_types = FALSE))
}

get_assumption <- function(keys, scenario = headline_scenario, default = "??") {
  if (is.null(scenario_assumptions) || nrow(scenario_assumptions) == 0) return(default)
  keys <- as.character(keys)
  dat <- scenario_assumptions
  
  # Schema 1: scenario, key/name/parameter, value.
  key_col <- intersect(c("key", "name", "parameter", "assumption", "field"), names(dat))[1]
  val_col <- intersect(c("value", "val", "setting"), names(dat))[1]
  if (!is.na(key_col) && !is.na(val_col)) {
    dat2 <- dat
    if ("scenario" %in% names(dat2)) {
      dat2 <- dat2 |> filter(.data$scenario %in% c(scenario, "all_fast", "all", "global", "default"))
    }
    hit_df <- dat2 |>
      filter(.data[[key_col]] %in% keys) |>
      slice(1)
    hit <- if (nrow(hit_df) > 0) hit_df[[val_col]][1] else NA
    if (length(hit) > 0 && !is.na(hit)) return(as.character(hit))
  }
  
  # Schema 2: one row per scenario with wide columns.
  if ("scenario" %in% names(dat)) {
    dat <- dat |> filter(.data$scenario %in% c(scenario, "all_fast", "all", "global", "default"))
  }
  for (k in keys) {
    if (k %in% names(dat) && nrow(dat) > 0) {
      val <- dat[[k]][1]
      if (!is.na(val)) return(as.character(val))
    }
  }
  default
}

macro_lines <- c(
  "% Auto-generated by scripts/07_make_outputs.R. Do not edit by hand.",
  paste0("% Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  tex_cmd("numTotalCVD", fmt_millions(headline_total_cvd, 2)),
  tex_cmd("numAllFastAverted", fmt_millions(headline_averted, 2)),
  tex_cmd("numAllFastPct", fmt_pct(headline_pct, 1)),
  tex_cmd("numStrokeICH", fmt_millions(headline_stroke_ich, 1)),
  tex_cmd("numQReduction", fmt_pct(headline_q_reduction, 1)),
  tex_cmd("numAllFastAROC", paste0(fmt_num(headline_aroc, 1), "\\%/yr")),
  tex_cmd("numSDGCatchup", paste0(fmt_num(SDG_CATCHUP_AROC, 1), "\\%/yr")),
  tex_cmd_text("numBPBaselineControl", get_assumption(c("bp_baseline_control", "baseline_bp_control", "htn_baseline_control", "baseline_control"))),
  tex_cmd_text("numBPTargetControl", get_assumption(c("bp_target_control", "target_bp_control", "htn_target_control", "target_control"))),
  tex_cmd_text("numBPStartYear", get_assumption(c("bp_start_year", "control_start_year", "htn_start_year"), default = "??")),
  tex_cmd_text("numBPTargetYear", get_assumption(c("bp_target_year", "control_target_year", "htn_target_year"), default = "??")),
  tex_cmd_text("numStatinsTargetCoverage", get_assumption(c("statins_target_coverage", "statin_target_coverage", "statins_target"))),
  tex_cmd_text("numSodiumReduction", get_assumption(c("sodium_reduction", "salt_reduction", "salteff", "sodium_target"))),
  tex_cmd_text("numTFATarget", get_assumption(c("tfa_target", "tfa_target_tfa", "target_tfa"))),
  tex_cmd("numScenarioStartYear", as.character(headline_window[1])),
  tex_cmd("numScenarioEndYear", as.character(headline_window[2])),
  tex_cmd_text("scenarioCVDPackage", "all-fast CVD package"),
  tex_cmd_text("numCervicalEnvelopeFlag", get_assumption(c("cervical_envelope_flag", "cervical_validation_flag"), default = "requires review")),
  tex_cmd_text("numCervicalMismatchAges", get_assumption(c("cervical_mismatch_ages", "cervical_calibration_ages"), default = "15--34"))
)

writeLines(macro_lines, DECK_NUMBERS_FILE)
message("  deck_numbers.tex ✓  (", normalizePath(DECK_NUMBERS_FILE), ")")

################################################################################
# 5  COMPLETE
################################################################################

message("\n── 07_make_outputs.R complete ───────────────────────────────────────────")
message("  Tables  : ", normalizePath(OUT_TABLES))
message("  Figures : ", normalizePath(OUT_FIGURES))
message("  Slides  : Rnw deck reads figures from outputs/figures and dynamic macros from deck_numbers.tex")
