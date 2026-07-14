
# Run GBD data processing script to get baseline_rates
source("021_get_base_rates_indonesia.R")

# Run TPS transition probabilities script to get dt_tps
source("022_get_tps_indonesia.R")

# Run script to make tps and bgmx trend forecast

source("023_get_tps_bgmx_indonesia.R")

