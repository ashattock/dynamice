###########################################################
# PREPARE
#
# Prepare for model simulation. Load model parameters and
# model input data.
#
###########################################################

# ---------------------------------------------------------
# Prepare model parameters
# ---------------------------------------------------------
prepare_params = function(data) {
  
  # MODEL PARAMETERS:
  #   gamma := recovery rate per timestep
  #   tstep := number of timesteps per year
  #     amp := amplification scale for seasonality
  #     ve1 := vaccine efficacy for first dose by each age group
  # ve2plus := vaccine efficacy for two and more doses
  #   vage2 := Age at second dose (country-specific)
  
  # Time and duration
  tstep	= 1000	# Number of timesteps in a year
  dinf  = 14		# Duration of infection (days)
  
  # age-dependent vaccine efficacy for first dose, based on a linear model (Hughes et al. 2020)
  ve1_intcp = 0.64598  # intercept of the linear model
  ve1_slope = 0.01485  # slope of the linear model, per month of age
  ve2plus   = 0.98     # vaccine efficacy for two and more doses
  
  age_ve1 = ve1_intcp + ve1_slope * 12 * c(1:(3*52)/52, 4:101)  # based on age in months
  age_ve1 = ifelse(age_ve1 >= ve2plus, ve2plus, age_ve1)
  
  # country-specific age at vaccination for MCV2
  vage2 = data$vax_age$mcv2
  
  # Parameters for Rcpp functions
  p = list(
    gamma   = 1 / (dinf * tstep/365),
    tstep   = tstep,
    amp     = 0.05,
    ve1     = age_ve1,
    ve2plus = ve2plus, 
    vage2   = vage2)
  
  return(p)
}

# ---------------------------------------------------------
# Load data for this country and this scenario
# ---------------------------------------------------------
prepare_data = function(country, scenario) {
  
  # ---- Coverage data ----
  
  # Coverage data for this scenario and country
  d = list(
    coverage_routine = load_coverage(country, scenario, "routine"), 
    coverage_sia     = load_coverage(country, scenario, "sia"))
  
  # ---- Load other input data for this country ----
  
  # All data rds files in data directory
  data_files = list.files(o$pth$data, pattern = "^data_.+\\.rds$")
  
  # Loop through files
  for (file in data_files) {
    
    # Shorthand reference for this data
    ref = str_remove_all(file, "(^data_|.rds$)")
    
    # Load data file
    data = readRDS(paste0(o$pth$data, file))
    
    # Filter datatable for this country
    if (is.data.frame(data))
      d[[ref]] = filter(data, country == !!country)
    
    # ... or select list element for this country
    if (!is.data.frame(data))
      d[[ref]] = data[[country]]
  }
  
  # ---- Possible data overwrites ----

  # Update timeliness if MCV1 not given at 39 weeks (9 months)
  if (d$vax_age$mcv1 != 39)
    d$timeliness[!is.na(age), timeliness := ifelse(age < d$vax_age$mcv1, 0, 1)]
  
  # Assume a fixed R0 if requested
  if (!is.na(o$fix_r0))
    d$rnought[, r0 := o$fix_r0]
  
  return(d)
}

# ---------------------------------------------------------
# Load coverage data
# ---------------------------------------------------------
load_coverage = function(country, scenario, type) {
  
  # Construct file name and path
  file_name = paste1(type, scenario)
  file_path = paste0(o$pth$coverage, file_name, ".csv")
  
  # Load coverage data and filter for years of interest
  coverage_data = fread(file_path) %>%
    filter(country_code %in% !!country, 
           year %in% o$analysis_years) %>%
    select(country = country_code, year, coverage)
  
  # Remove trivial values for SIA
  if (type == "sia") {
    coverage_data %<>%
      filter(coverage > 0)
  }
  
  return(coverage_data)
}

