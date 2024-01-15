###########################################################
# OPTIONS
#
# Set key options for all things model related. The output
# of this function, o (a list), lives in the global environment,
# so can be referenced throughout the pipeline.
#
###########################################################

# ---------------------------------------------------------
# Set model options and assumptions
# ---------------------------------------------------------
set_options = function(do_step = NA, quiet = FALSE) {

  if (!quiet) message("* Setting options")

  # Several global R settings to make life easier
  default_R_options()  # See auxiliary.R

  # Initiate options list
  o = list(do_step = do_step)

  # Prepare output directory system
  o = prepare_dirs(o)  # See directories.R

  # ---- Simulation settings ----

  # Countries to simulate (ISO3 or "all")
  o$countries = "all" # qc(AGO, CHN, COD, ETH, IDN, IND, MDG, MWI, NGA, PAK, PHL, SOM, UGA, UKR)
  
  # Scenarios to simulate (scenario ID or "all")
  o$scenarios = c("nomcv", "mcv1_mcv2_sia")
  
  # Years to analyse
  o$years = 1974 : 2024
  
  # Age range 
  o$ages = 0 : 100
  
  # Assume a fixed R0 for the central run
  o$fix_r0 = 15.9  # Set to NA to turn off
  
  # Vary R0 for uncertainty simulations
  o$vary_r0 = seq(6, 26, 2)
  
  # Disability for non-fatal case
  #
  # TODO: Ask what this represents, assuming the product of disability
  #       weight and average duration of infection until recovery
  o$disability_weight = 0.002
  
  # Model metrics to report
  o$metrics = c("deaths", "dalys", "yll")
  
  # ---- Cluster settings ----
  
  # Detect user - required for checking cluster jobs
  o$user = Sys.info()[["user"]]
  
  # Flag for overwriting any existing simulations
  o$overwrite = FALSE
  
  # Time allocated for each cluster job
  o$job_time  = "00:05:00"  # Use "HH:MM:SS" format
  o$job_queue = "30min"  # Queue to use (check out scicore wiki if unfamiliar)
  
  # Memory to allocate for each cluster job
  o$job_memory = "500MB"
  
  # Set an upper limit for jobs that can be run at any one time
  o$job_limit = 2000
  
  # Define names for cluster log and error files
  o$log_file = "scicore_log.txt"
  o$err_file = "scicore_error.txt"
  
  # Action to take if user is already running cluster jobs
  o$cluster_conflict_action = "error"  # Set to 'none' to turn off
  
  # ---- Results and plotting flags ----
  
  # Produce main output files
  o$do_output_files = TRUE

  # Turn figures on or off
  o$plot_coverage       = FALSE  # TODO: Not yet integrated
  o$plot_everything     = FALSE  # TODO: Not yet integrated
  o$plot_burden_averted = TRUE

  # ---- Plotting settings ----

  # Saved figure size
  o$save_width  = 14
  o$save_height = 10

  # Units of figures sizes
  o$save_units = "in"

  # Plotting resolution (in dpi)
  o$save_resolution = 300

  # Image format for saving figure
  #
  # NOTE: Use a character vector to save with multiple formats at once
  o$figure_format = "png" # Classic options: "png", "pdf", or "svg"
  
  # ---- Parse key input options ----
  
  # Parse countries and scenarios
  o = parse_countries(o, quiet)
  o = parse_scenarios(o, quiet)
  
  return(o)
}

# ---------------------------------------------------------
# Parse countries to simulate
# ---------------------------------------------------------
parse_countries = function(o, quiet) {
  
  # All countries - defined in csv file
  o$countries_all = fread(paste0(o$pth$config, "countries.csv"))$country
  
  # Parse the 'all' option
  if (length(o$countries) == 1 && o$countries == "all") 
    o$countries = o$countries_all
  
  # Throw an error if any unknown countries defined
  unknown = setdiff(o$countries, o$countries_all)
  if (length(unknown) > 0)
    stop("Unknown countries selected: ", paste(unknown, collapse = ", "))
  
  # Display which country(ies) we are running
  if (quiet == FALSE) {
    
    # Number of countries we are running
    k = length(o$countries)
    n = length(o$countries_all)
    
    # All countries
    if (k == n)
      message(" > Running all ", n, " countries")
    
    # One country
    if (k == 1)
      message(" > Running 1 country: ", o$countries)
    
    # Multiple countries, but not all
    if (k > 1 && k < n)
      message(" > Running ", k, " countries: ", 
              paste0(o$countries, collapse = ", "))
  }
  
  return(o)
}

# ---------------------------------------------------------
# Parse scenarios to simulate
# ---------------------------------------------------------
parse_scenarios = function(o, quiet) {
  
  # All feasible scenarios - defined in directories.R
  all_scenarios = fread(paste0(o$pth$config, "scenarios.csv"))$scenario
  
  # Parse the 'all' option
  if (length(o$scenarios) == 1 && o$scenarios == "all") 
    o$scenarios = all_scenarios
  
  # Throw an error if any unknown scenarios defined
  unknown = setdiff(o$scenarios, all_scenarios)
  if (length(unknown) > 0)
    stop("Unknown scenarios selected: ", paste(unknown, collapse = ", "))
  
  # Display which country(ies) we are running
  if (quiet == FALSE) {
    
    # Number of countries we are running
    k = length(o$scenarios)
    n = length(all_scenarios)
    
    # All countries
    if (k == n)
      message(" > Running all ", n, " vaccine scenarios")
    
    # One country
    if (k == 1)
      message(" > Running 1 vaccine scenario: ", o$scenarios)
    
    # Multiple countries, but not all
    if (k > 1 && k < n)
      message(" > Running ", k, " vaccine scenarios: ", 
              paste0(o$scenarios, collapse = ", "))
  }
  
  return(o)
}

