###########################################################
# RESULTS
#
# Call plotting functions as defined by o$plot_xxx flags (see
# options.R). All plotting functions themselves live in
# plotting.R.
#
###########################################################

# ---------------------------------------------------------
# Call plots as defined by o$plot_xxx flags
# ---------------------------------------------------------
run_results = function() {
  
  # Only continue if specified by do_step
  if (!is.element(3, o$do_step)) return()
  
  message("* Producing results")
  
  # Create main output files
  if (o$do_output_files) {
    
    # Output files for VIMC
    vimc_output()
    
    # Output files for EPI50
    epi50_output()
  }
  
  # Plot coverage
  # if (o$plot_coverage)
  #   plot_coverage_all()
  
  # Plot everything
  # if (o$plot_everything)
  #   plot_everything()
  
  # Plot deaths and DALYs averted
  if (o$plot_burden_averted)
    plot_burden_averted()
}

# ---------------------------------------------------------
# Main output csv files for VIMC
# ---------------------------------------------------------
vimc_output = function() {
  
  message(" > Producing VIMC output files")
  
  # ---- Format results template ----
  
  # Load central estimates template
  template_dt = readRDS(paste0(o$pth$input, "template.rds"))
  
  # Trivial columns for which we are to insert values
  keep_cols = names(which(colSums(is.na(template_dt)) != nrow(template_dt)))
  
  # Redefine year scope
  vimc_dt = template_dt %>%
    filter(country %in% o$countries) %>%
    select(-year) %>%
    unique() %>%
    expand_grid(year = o$years) %>%
    select(all_of(keep_cols)) %>%
    as.data.table()
  
  # Loop through scenarios
  for (scenario in o$scenarios) {
    
    # ---- Load all results for this scenario
    
    # We'll load results for all countrues for this scenario
    load_names = paste1(o$countries_all, scenario)
    load_files = paste0(o$pth$compiled, load_names, ".rds")
    
    # Load these results
    results_dt = rbindlist(lapply(load_files, readRDS))
    
    # ---- Central estimates ----
    
    # Central results for this scenario 
    model_dt = results_dt %>%
      filter(scenario == !!scenario, 
             is.na(r0)) %>%  # This identifies central estimate
      select(-scenario, -r0) %>%
      pivot_wider(names_from  = "metric", 
                  values_from = "value") %>%
      as.data.table()
    
    # Construct file name path to save to
    save_name = paste1("total", scenario)
    save_file = paste0(o$pth$results, save_name, ".csv")
    
    # Save all country results to file
    fwrite(model_dt, file = save_file)
    
    # Append model results to VIMC template
    output_dt = vimc_dt %>%
      left_join(y  = model_dt, 
                by = c("year", "age", "country"))
    # select(all_of(names(template_dt)))
    
    # Construct file name path to save to
    save_name = paste1("central", scenario)
    save_file = paste0(o$pth$results, save_name, ".csv")
    
    # Save results to file
    fwrite(output_dt, file = save_file)
  }
  
  # ---- Uncertainty simulations ----
  
  # TODO: I'm not sure how LSHTM handle uncertainty
}

# ---------------------------------------------------------
# Output for EPI50 analysis
# ---------------------------------------------------------
epi50_output = function() {
  
  message(" > Producing EPI50 output file")
  
  # Convert scenarios into named list
  scen_list = o$scenarios %>%
    setNames(c("no_vaccine", "vaccine")) %>%
    as.list()
  
  # Function for loading results of given scenario
  load_fn = function(ref) {
    
    # File name and path for loading
    load_name = paste1("total", scen_list[[ref]])
    load_file = paste0(o$pth$results, load_name, ".csv")
    
    # Load file and set to long format
    result_dt = fread(load_file) %>%
      select(country, year, age, deaths, dalys) %>%
      mutate(scenario = ref, .before = 1)  %>%
      pivot_longer(cols = c("deaths", "dalys"), 
                   names_to  = "metric", 
                   values_to = "value") %>%
      arrange(country, year, age, metric) %>%
      as.data.table()
    
    return(result_dt)
  }
  
  # Apply loading function and format into single datatable
  epi50_dt = names(scen_list) %>%
    lapply(load_fn) %>%
    rbindlist()
  
  # Construct file name path to save to
  save_file = paste0(o$pth$output, "epi50_dynamice_results.rds")
  
  # Save results to file
  saveRDS(epi50_dt, file = save_file)
}

# ---------------------------------------------------------
# Easy load of all central results
# ---------------------------------------------------------
load_central_results = function() {
  
  # TODO: Throw helpful warning if vimc_output() has not yet been called
  
  # Function for loading a single file
  load_fn = function(scenario) {
    
    # File name and path
    load_name = paste1("central", scenario)
    load_file = paste0(o$pth$results, load_name, ".csv")
    
    # Load results and append scenario details
    result_dt = fread(load_file) %>%
      mutate(scenario = scenario, 
             .before = 1) %>%
      select(-disease)
   
    return(result_dt) 
  }
  
  # Call loading function for each scenario
  central_dt = lapply(o$scenarios, load_fn) %>%
    rbindlist() %>%
    # Convert to tidy format...
    pivot_longer(cols = o$metrics, 
                 names_to = "metric") %>%
    # Tidy up...
    select(country, country_name, scenario, 
           year, age, metric, value) %>%
    arrange(metric, scenario, country, year, age) %>%
    as.data.table()
  
  return(central_dt)
}

