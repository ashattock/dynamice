###########################################################
# LAUNCH
#
# Main launch function for running DyanMICE measles model for
# WHO EPI50 analysis.
#
###########################################################

# Set working directory to sourced file
if (interactive()) setwd(getSrcDirectory(function() {}))

# Load all required packages and functions
source("dependencies.R")

message("Running DynaMICE for EPI50 analysis")

# Set options (see options.R)
o = set_options(do_step = 1)

# Step 1) Run vaccine scenarios
run_scenarios()  # See scenarios.R

# Step 2) Produce outputs and figures
run_results()  # See results.R

# Finish up
message("* Finished!")

