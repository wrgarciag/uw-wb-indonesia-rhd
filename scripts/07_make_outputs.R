# ==============================================================================
# 8. RUN incremental and benefit analysis
# ==============================================================================

inc <- data.table(
  year             = years,
  rhd_deaths_avert = ref$rhd_deaths - sap$rhd_deaths,
  severe_avert     = ref$severe     - sap$severe,
  inc_cost         = sap$cost - ref$cost
)
inc[, deaths_avert := rhd_deaths_avert]

df      <- 1 / (1 + disc_rate)^(years - years[1])          # discount factors, base 2021
gdp_pc  <- gdp_pc_2019 * (1 + gdp_growth)^(years - 2019)
vsl     <- vsl_mult * gdp_pc

benefit_vsl   <- inc$deaths_avert * vsl        # value of health (dominant term)
benefit_gdp   <- inc$deaths_avert * gdp_pc     # productivity proxy
benefit_total <- benefit_vsl + benefit_gdp

tot_cost_disc    <- sum(inc$inc_cost  * df)
tot_benefit_disc <- sum(benefit_total * df)
tot_deaths_avert <- sum(inc$deaths_avert)
bcr            <- tot_benefit_disc / tot_cost_disc
net_ben        <- tot_benefit_disc - tot_cost_disc
cost_per_death <- tot_cost_disc / tot_deaths_avert
cost_per_daly  <- tot_cost_disc / (sum(inc$deaths_avert * df) * dalys_per_death)

# ==============================================================================
# 9. REPORT
# ==============================================================================
cat("\n============================================================\n")
cat(sprintf(" RHD SECONDARY PREVENTION (screening + SAP)\n Horizon %d-%d | scale-up ramp %d-%d\n",
            min(years), max(years), ramp_start, ramp_end))
cat("============================================================\n\n")

cat("--- BUDGET IMPACT (US$ millions) ---\n")
bi <- sap[, .(year, screening = c_screen/1e6, SAP = c_sap/1e6,
              HF = c_hf/1e6, surgery = c_surg/1e6, total = cost/1e6)]
bi <- merge(bi, ref[, .(year, total_ref = cost/1e6)], by = "year")
bi[, incremental := total - total_ref]
print(bi[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 1) else x)])

cat("\n--- HEALTH IMPACT (cumulative) ---\n")
cat(sprintf("  RHD deaths averted                   : %s\n",
            format(round(tot_deaths_avert), big.mark = ",")))
cat(sprintf("  Severe RHD cases averted (2030 stock): %s\n",
            format(round(inc$severe_avert[n_years]), big.mark = ",")))

cat("\n--- ECONOMIC RESULTS (3% discounting; 2019 US$) ---\n")
cat(sprintf("  Total incremental cost    : $%.2f billion\n", tot_cost_disc/1e9))
cat(sprintf("  Total monetised benefit   : $%.2f billion\n", tot_benefit_disc/1e9))
cat(sprintf("  Net benefit               : $%.2f billion\n", net_ben/1e9))
cat(sprintf("  Benefit-cost ratio        : %.2f\n", bcr))
cat(sprintf("  Cost per death averted    : $%s\n", format(round(cost_per_death), big.mark = ",")))
cat(sprintf("  Cost per DALY averted      : $%s\n\n", format(round(cost_per_daly), big.mark = ",")))

results <- list(budget_impact = bi, incremental = inc,
                summary = data.table(bcr, net_benefit = net_ben,
                                     total_cost = tot_cost_disc, total_benefit = tot_benefit_disc,
                                     deaths_averted = tot_deaths_avert,
                                     cost_per_death, cost_per_daly))
saveRDS(results, "rhd_results.rds")

fwrite(ref, "rhd_reference.csv")
fwrite(bi,  "rhd_budget_impact.csv")
fwrite(inc, "rhd_incremental.csv")
cat("Saved: rhd_budget_impact.csv, rhd_incremental.csv, rhd_results.rds\n")

