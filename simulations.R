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
  
  # ---- Run simulations ----
  
  # Number of jobs to be run
  n_jobs = nrow(sims) # sum(sims$run)
  
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
  
  # Load all scenarios from config file
  scenarios = fread(paste0(o$pth$config, "scenarios.csv"))
  
  # Interpretatiion of set_routine
  #  0: no routine MCV
  #  1: MCV1 only
  #  2: MCV1 + MCV2
  
  # Interpretatiion of set_sia
  #  0: no SIA
  #  1: random reach (baseline assumption)
  #  2: 7.7% less likely to be reached at national level
  #  3: zero-dose first
  #  4: already-vaccinated first
  #  5: 7.7% less likely to be reached at subnational level
  
  # Grid of countries and scenarios to run
  sims = expand_grid(country  = o$countries, 
                     scenario = o$scenarios) %>%
    # Append scenario vaccination details...
    left_join(y  = scenarios, 
              by = "scenario") %>%
    select(-scenario_name) %>%
    # Append scenario ID...
    mutate(id = get_simulation_id(.)) %>%
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
  # 
  # # Save scenario dataframe to file
  # saveRDS(sims, file = paste0(o$pth$sim_input, "simulations_", timeframe, ".rds"))
  
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
  
  # Combine details to create unique simulation ID
  ids = paste1(
    sim$country, 
    sim$set_routine, 
    sim$set_sia, 
    sim$scenario)
  
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
      sel_scns <- 1:length(scenarios) # main assumption
      set_sia  <- c (0, 0, 0, 2, 2, 0, 2, 3, 4, 3, 4, 3, 0)
    } else {
      sel_scns <- c(2,3,4,5,7) # evaluate different SIA assumptions
      set_sia [c(4,5,7,14)] <- isia
    }
    
    for (index in sel_scns){
      
      scenario_name  <- scenarios [index]
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
        scenario_name              = scenarios [index],
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
      
      scenario_name   <- scenarios [index]
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
        scenario_name              = scenarios [index],
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

