###########################################################
# SET DIRECTORIES
#
# Set and get directories in one place in the name of consistency
# and ease. Creates any directories that do not currently exist.
#
# OUTPUTS:
#	- A list of relevant directories (within o$pth) which can
#   be referenced throughout in the pipeline
#
###########################################################

# ---------------------------------------------------------
# Define paths for project inputs and outputs
# ---------------------------------------------------------
prepare_dirs = function(o) {

  # Initiate list with reference to main code directory
  #
  # NOTE: We should already be in this code directory
  pth = list(code = getwd())
  
  # Functino for setting path relative to cwd
  set_pth = function(...)
    dir = file.path(pth$code, ...)
  
  # ---- Input and configuration directories ----

  # Directories for all input files
  pth$config   = set_pth("config")
  pth$input    = set_pth("input")
  pth$coverage = set_pth("input", "coverage")

  # ---- Output directories ----
  
  # Parent path for all output files
  pth$output   = set_pth("output")

  # Path to model output, results, and figures
  pth$sims     = set_pth("output", "1_simulations")
  pth$burden   = set_pth("output", "2_burden")
  pth$compiled = set_pth("output", "3_compiled")
  pth$results  = set_pth("output", "4_results")
  pth$figures  = set_pth("output", "5_figures")
  
  # Directory for cluster logs
  pth$log = set_pth("log")

  # Append paths to o list
  o = make_dirs(o, pth)

  return(o)
}

# ---------------------------------------------------------
# Make all output directories and append to o list
# ---------------------------------------------------------
make_dirs = function(o, pth) {

  # Iterate through dirs
  for (pth_name in names(pth)) {
    this_pth = pth[[pth_name]]

    # If it does not already exist, create it
    if (!dir.exists(this_pth))
      dir.create(this_pth, recursive = TRUE)

    # Add a file separator to end of dir path
    pth[[pth_name]] = paste0(this_pth, file_sep())
  }

  # Append to o list
  o$pth = pth

  return(o)
}

