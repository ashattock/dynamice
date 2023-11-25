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
set_options = function(do_step = NA) {

  message("* Setting options")

  # Several global R settings to make life easier
  default_R_options()  # See auxiliary.R

  # Initiate options list
  o = list(do_step = do_step)

  # Prepare output directory system
  # o = prepare_dirs(o)  # See directories.R

  # ---- Country scope ----

  # list countries for analysis
  o$countries = qc(AGO, CHN, COD, ETH, IDN, IND, MDG,
                   MWI, NGA, PAK, PHL, SOM, UGA, UKR)

  # ---- Time settings ----

  # Years to analyse
  o$analysis_years = 1980 : 2024

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
