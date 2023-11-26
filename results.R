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
  if (!is.element(2, o$do_step)) return()
  
  message("* Producing results")
  
  # Create output files for VIMC
  if (o$do_vimc_files)
    vimc_output()
  
  browser()
  
  # Plot coverage
  if (o$plot_coverage)
    plot_coverage_all()
  
  # Plot everything
  #
  # TODO: Seperate this out to individual plots
  if (o$plot_everything)
    plot_everything()
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
    load_names = paste1(o$countries, scenario)
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

