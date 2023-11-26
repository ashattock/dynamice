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

  # Countries to simulate
  o$countries = qc(AGO, CHN, COD, ETH, IDN, IND, MDG,
                   MWI, NGA, PAK, PHL, SOM, UGA, UKR)
  
  # Scenarios to simulate
  o$scenarios = qc(nomcv, mcv1)
  
  # Years to analyse
  o$years = 1980 : 2024
  
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

  # ---- Data settings ----
  
  # Force re-construct coverage data
  o$reload_coverage = FALSE
  
  # ---- Cluster settings ----
  
  # Detect user - required for checking cluster jobs
  o$user = Sys.info()[["user"]]
  
  # Flag for overwriting any existing simulations
  o$overwrite = FALSE
  
  # Time allocated for each cluster job
  o$job_time  = "00:15:00"  # Use "HH:MM:SS" or "D-HH:MM:SS" format
  o$job_queue = "30min"  # Queue to use (check out scicore wiki if unfamiliar)
  
  # Memory to allocate for each cluster job
  o$job_memory = "1GB"
  
  # Set an upper limit for jobs that can be run at any one time
  o$job_limit = 2000
  
  # Define names for cluster log and error files
  o$log_file = "scicore_log.txt"
  o$err_file = "scicore_error.txt"
  
  # Action to take if user is already running cluster jobs
  o$cluster_conflict_action = "error"  # Set to 'none' to turn off
  
  # ---- Plotting flags ----

  # Turn figures on or off
  o$plot_coverage   = TRUE
  o$plot_everything = TRUE

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

  return(o)
}

