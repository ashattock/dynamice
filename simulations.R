############################################################
# SIMULATIONS
#
# Create, simulate, and post-process all simulations.
#
############################################################

# ---------------------------------------------------------
# Parent function for creating and running all simulations
# ---------------------------------------------------------
run_simulations = function() {
  
  # Only continue if specified by do_step
  if (!is.element(1, o$do_step)) return()
  
  message("* Running vaccine simulations")
  
  # Prepare coverage input data
  if (o$reload_coverage)
    create_coverage()
  
  # Generate full set of simulations to run
  sims = get_simulations()
  
  # ---- Submit simulations to cluster ----
  
  # Number of jobs to be run
  n_jobs = nrow(sims) # sum(sims$run)
  
  browser()
  
  # Skip if nothing to run
  if (n_jobs > 0) {
    
    # Submit all jobs to the cluster (see myRfunctions.R)
    submit_cluster_jobs(n_jobs, "submit.sh", "run_sim")
    
    # Throw an error if any cluster jobs failed (see myRfunctions.R)
    stop_if_errors(o$pth$log, o$err_file, err_tol = 1)
    
    # Remove all log files if desired
    if (o$rm_cluster_log) 
      unlink(paste0(o$pth$log, "*"), force = TRUE)
  }
  
  browser()
  
  # ---- Concatenate output ----
  
  # Aggregate results for each branch and indicator
  # run_aggregate(sims, timeframe)  # See aggregate.R
}

# ---------------------------------------------------------
# Generate full set of simulations to run
# ---------------------------------------------------------
get_simulations = function() {
  
  # ---- Sanity check on scenario selection ----
  
  # Load all scenarios from config file
  scenarios = fread(paste0(o$pth$config, "scenarios.csv"))
  
  # Throw an error if any unrecognised scenarios
  unknown = setdiff(o$scenarios, scenarios$scenario)
  if (length(unknown) > 0)
    stop("Unrecognised scenarios defined: ", paste(unknown, collapse = ", "))
  
  # ---- Full set of simulations ----
  
  # Grid of countries and scenarios to run
  sims = expand_grid(country  = o$countries, 
                     scenario = o$scenarios, 
                     r0 = c(NA, o$vary_r0)) %>%
    # Append scenario ID...
    mutate(id = get_simulation_id(.), 
           .before = 1) %>%
    arrange(scenario, country) %>%
    as.data.table()
  
  # ---- Skip existing sims ----
  
  # message(" > Identifying previously completed simulations")
  # 
  # # Extract IDs of sims that have already been run for each country
  # exist_pth = o$pth[[paste1("post", timeframe)]]
  # exist_id  = str_remove(list.files(exist_pth), ".rds$")
  # 
  # # Logical whether simulation should be run / rerun
  # run_sim = !(sims$id %in% exist_id)
  # if (o$overwrite) run_sim[] = TRUE
  # 
  # # Skip any existing sims (unless overwriting)
  # sims %<>%
  #   cbind(run = run_sim) %>%
  #   mutate(job_num = cumsum(run * 1), 
  #          job_num = ifelse(run, job_num, NA))
  
  # Save scenario dataframe to file
  saveRDS(sims, file = paste0(o$pth$sims, "all_simulations.rds"))
  
  # ---- Display number of sims ----
  
  # # Number of sims
  # n_total = nrow(sims)
  # n_run   = sum(sims$run)
  # 
  # # Report total number of sims
  # message(" > Total number of simulations: ", thou_sep(n_total))
  # 
  # # Report number of sims we'll run now
  # message("  - Skipping: ", thou_sep(n_total - n_run))
  # message("  - Simulating: ", thou_sep(n_run))
  
  return(sims)
}

# ---------------------------------------------------------
# Create simulation ID convention 
# ---------------------------------------------------------
get_simulation_id = function(sim) {
  
  # Combine scenario details to create sim ID
  ids = sim %>%
    mutate(r0_str = ifelse(is.na(r0), "def", r0),
           r0_str = str_pad(r0_str, 2, pad = "0")) %>%
    mutate(id = paste1(country, r0_str, scenario)) %>%
    pull(id)
  
  return(ids)
}
  
# ---------------------------------------------------------
# All steps to actually simulate the model
# ---------------------------------------------------------
run_sim = function(job_id) {
  
  # ---- Details of this simulation ----
  
  # Load full scenario vaccination details
  scenarios_dt = fread(paste0(o$pth$config, "scenarios.csv"))
  
  # Load full set of simulations
  sims_dt = readRDS(paste0(o$pth$sims, "all_simulations.rds"))
  
  # Select simulation assocaited with this job_id
  sim = sims_dt %>%
    slice(job_id) %>% # filter(job_num == job_id)
    left_join(y  = scenarios_dt, 
              by = "scenario") %>%
    select(-scenario_name)
  
  message(" > Running ", sim$id)
  
  # ---- Load data for this simulation ----
  
  # Prepare model input data
  data = prepare_data(sim)  # See prepare.R
  
  # ---- Simulate model ----
  
  # Run DynaMICE model
  run_model(sim, data)  # See model.R
  
  # Estimate disease burden from model outcomes
  run_burden(sim, data)  # See model.R
}

