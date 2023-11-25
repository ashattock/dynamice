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

  # Plot coverage
  if (o$plot_coverage)
    plot_coverage_all()

  # Plot everything
  #
  # TODO: Seperate this out to individual plots
  if (o$plot_everything)
    plot_everything()
}

