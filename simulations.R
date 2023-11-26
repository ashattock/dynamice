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
  
  # Generate full set of simulations to run
  sims = get_simulations()
  
  # ---- Run simulations ----
  
  # Number of jobs to be run
  n_jobs = sum(sims$run)
  
  browser()
  
  # Skip if nothing to run
  if (n_jobs > 0) {
    
    # Specify wrapper function for running this set of simulations
    sim_fn = paste0("run_sim_", timeframe)
    
    # Submit all jobs to the cluster (see myRfunctions.R)
    submit_cluster_jobs(n_jobs, "bash_submit.sh", sim_fn)
    
    # Throw an error if any cluster jobs failed (see myRfunctions.R)
    stop_if_errors(o$pth$log, o$err_file, err_tol = 1)
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
  
  browser()
  
  # Grid of countries, EIR, and seeds
  sims = expand_grid(country  = o$countries, 
                     eir      = o$sim$eir, 
                     seed     = 1 : o$sim$n_seeds[[timeframe]], 
                     scenario = run_scenarios) %>%
    # Append scenario ID...
    mutate(id = get_simulation_id(.)) %>%
    arrange(scenario, country, eir, seed) %>%
    as.data.table()
  
  browser()
  
  # prepare coverage inputs
  vac_strategies <- c(
    "nomcv",                 # (1) no vaccination
    "mcv1",                  # (2) MCV1 only
    "mcv1-mcv2",             # (3) MCV1 + MCV2
    "mcv1-mcv2-sia",         # (4) MCV1 + MCV2 + SIA
    "mcv1-sia",              # (5) MCV1 + SIA
    "mcv1-mcv2alt1",         # (6) MCV1 + MCV2(early intro, fast rollout)
    "mcv1-mcv2alt1-sia",     # (7) MCV1 + MCV2(early intro, fast rollout) + SIA
    "mcv1-mcv2-siaalt1",     # (8) MCV1 + MCV2 + SIA(zero dose first)
    "mcv1-mcv2-siaalt2",     # (9) MCV1 + MCV2 + SIA(already vaccinated first)
    "mcv1-siaalt1",          # (10) MCV1 + SIA(zero dose first)
    "mcv1-siaalt2",          # (11) MCV1 + SIA(already vaccinated first)
    "mcv1-mcv2alt1-siaalt1", # (12) MCV1 + MCV2(early intro) + SIA(zero dose first)
    "mcv1-mcv2alt2")         # (13) MCV1 + MCV2(early intro, gradual rollout)
  
  # set SIAs implementation method for each scenario
  # additional modifications needed in Rcpp for extent-specific approaches
  # ex. random reach for subnational campaign & 7.7% at national level for national and rollover-nat campaigns
  # 0: no SIA
  # 1: random reach (baseline assumption)
  # 2: 7.7% less likely to be reached at national level
  # 3: zero-dose first
  # 4: already-vaccinated first
  # 5: 7.7% less likely to be reached at subnational level
  set_sia         <- c (0, 0, 0, 2, 2, 0, 2, 3, 4, 3, 4, 3, 0)
  
  # set routine vaccination parameters to distinguish between scenarios
  # 0: no routine MCV
  # 1: MCV1 only
  # 2: MCV1 + MCV2
  set_routine <- c (0, 1, 2, 2, 1, 2, 2, 2, 2, 1, 1, 2, 2)
  
  # ---- Generate vaccine coverage scenarios ----
  
  # Prepare coverage input data - update when the data are changed
  if (o$reload_coverage)
    lapply(vac_strategies, create_vaccine_coverage_routine_sia)
  
  
  
  
  
  # ---- Skip existing sims ----
  
  message(" > Identifying previously completed simulations")
  
  # Extract IDs of sims that have already been run for each country
  exist_pth = o$pth[[paste1("post", timeframe)]]
  exist_id  = str_remove(list.files(exist_pth), ".rds$")
  
  # Logical whether simulation should be run / rerun
  run_sim = !(sims$id %in% exist_id)
  if (o$overwrite) run_sim[] = TRUE
  
  # Skip any existing sims (unless overwriting)
  sims %<>%
    cbind(run = run_sim) %>%
    mutate(job_num = cumsum(run * 1), 
           job_num = ifelse(run, job_num, NA))
  
  # Save scenario dataframe to file
  saveRDS(sims, file = paste0(o$pth$sim_input, "simulations_", timeframe, ".rds"))
  
  # ---- Display number of sims ----
  
  # Number of sims
  n_total = nrow(sims)
  n_run   = sum(sims$run)
  
  # Report total number of sims
  message(" > Total number of simulations: ", thou_sep(n_total))
  
  # Report number of sims we'll run now
  message("  - Skipping: ", thou_sep(n_total - n_run))
  message("  - Simulating: ", thou_sep(n_run))
  
  return(sims)
}
  
# ---------------------------------------------------------
# Create simulation ID convention 
# ---------------------------------------------------------
get_simulation_id = function(sim) {
  
  browser()
  
  # Format IDs from details, padding EIR and seed values
  ids = paste(sim$country, 
              sim$eir  %>% str_pad(4, pad = "0"), 
              sim$seed %>% str_pad(3, pad = "0"), 
              sim$scenario, 
              sep = "_")
  
  return(ids)
}
  
  
# ---------------------------------------------------------
# All steps to actually simulate the model
# ---------------------------------------------------------
run_sim = function(job_id) {  
  
  browser()
  
  # ---- Run model for SIA activities ----
  
  # run model by different SIA methods
  for (isia in c(1,2,5)){
    
    if (isia == 2){
      sel_scns <- 1:length(vac_strategies) # main assumption
      set_sia  <- c (0, 0, 0, 2, 2, 0, 2, 3, 4, 3, 4, 3, 0)
    } else {
      sel_scns <- c(2,3,4,5,7) # evaluate different SIA assumptions
      set_sia [c(4,5,7,14)] <- isia
    }
    
    for (index in sel_scns){
      
      scenario_name  <- vac_strategies [index]
      message(" - ", scenario_name)
      
      scenario_number <- sprintf("scenario%02d", index)
      
      # run model and estimate cases
      burden_estimate_file <- runScenario_rcpp (
        scenario_name          = scenario_name,
        routine            = set_routine[index],
        using_sia              = set_sia[index])
      
      browser()
      
      # separately estimate dalys
      burden_estimate_file <- paste0 ("central_burden_estimate_",
                                      scenario_name, ".csv")
      
      # merge outputs into csv files
      get_burden_estimate(
        # vaccine_coverage_subfolder = var$vaccine_coverage_subfolder,
        scenario_name              = vac_strategies [index],
        save_scenario              = scenario_number,
        routine                = set_routine [index],
        using_sia                  = set_sia         [index],
        folder_date                = "20230401")
    }
    
    browser()
    
    # move files to a specified folder
    res_files <- list.files(o$pth$central)
    dir.create (paste0 ("previous_res/20230401/siareach_", isia, "/"))
    file.rename (from = paste0 (o$pth$central, res_files),
                 to = paste0 ("previous_res/20230401/siareach_", isia, "/", res_files))
  }
  
  # analyse burden estimates under different R0
  set_sia <- c (0, 0, 0, 2, 2, 0, 2, 3, 4, 3, 4, 3, 0) # set back to the oringinal assumptions
  
  for (ir0 in c(seq(6,26,2))) {
    # vary R0 values
    data_r0 <- copy (data_r0) [ , r0 := ir0]
    source ("R/functions_rcpp.R")
    
    if (ir0 %in% c(6,16,26)){
      sel_scns <- c(2,3,4,5)
    } else {
      sel_scns <- 2
    }
    
    for (index in sel_scns){
      
      scenario_name   <- vac_strategies [index]
      scenario_number <- sprintf ("scenario%02d", index)
      
      # run model and estimate cases
      burden_estimate_file <- runScenario_rcpp (
        scenario_name = scenario_name,
        routine   = set_routine[index],
        using_sia     = set_sia[index])
      
      # separately estimate deaths
      burden_estimate_file <- paste0 ("central_burden_estimate_",
                                      scenario_name, ".csv")
      
      # merge outputs into csv files
      get_burden_estimate(
        scenario_name              = vac_strategies [index],
        save_scenario              = sprintf ("scenario%02d", index),
        routine                = set_routine [index],
        using_sia                  = set_sia         [index],
        folder_date                = "20230401")
      
      # rename file
      file.rename (from = paste0(o$pth$central, burden_estimate_file),
                   to = paste0(o$pth$central,
                               "central_burden_estimate_",
                               scenario_name,
                               "_r0-", ir0, ".csv"))
    }
  }
  
  browser()
  
  # move files to a specified folder
  res_files <- list.files(o$pth$central)
  dir.create (paste0 ("previous_res/20230401/siareach_2/senanl_r0"))
  file.rename (from = paste0(o$pth$central, res_files),
               to = paste0("previous_res/20230401/siareach_2/senanl_r0/", res_files))
}

