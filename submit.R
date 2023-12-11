############################################################
# SUBMIT
#
# Submit cluster tasks from a cluster array.
#
############################################################

source("dependencies.R")

# Extract input from bash file
args = commandArgs(trailingOnly = TRUE)

# Name arguments when provided from bash
if (length(args) > 0) {
  job_fn = as.character(args[1])
  job_id = as.numeric(args[2])
  
} else {  # Otherwise set defaults - use for testing and debugging
  job_fn = "run_sim"
  job_id = 49
}

# Reload simulation options into global environment
o = set_options(quiet = TRUE) # See options.R

# Run job_fn function for job_id and catch any error
tryCatch(
  get(job_fn)(job_id),
  
  # Error handler
  error = function(e) {
    
    # Concatenate (ideally useful) error message
    err_message = paste0(" ! ", e$message, " (array ID: ", job_id, ")")
    
    message(err_message)
    
    # Append this error message to an error log file
    write(err_message, file = paste0(o$pth$log, o$err_file), append = TRUE)
  },
  
  # Close up
  finally = message("Closing cluster job")
)

