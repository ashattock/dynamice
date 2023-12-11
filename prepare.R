###########################################################
# PREPARE
#
# Prepare for model simulation. Load model parameters and
# model input data.
#
###########################################################

# ---------------------------------------------------------
# Prepare coverage data by splitting into routine and SIA
# ---------------------------------------------------------
run_prepare = function() {
  
  # Only continue if specified by do_step
  if (!is.element(1, o$do_step)) return()
  
  message("* Preparing model resources")
  
  # Function for constructing file paths
  get_path = function(x)
    paste0(o$pth$coverage, paste1(x, scenario), ".csv")
  
  # Function for converting age in years to week reference
  #
  # For age < 3 years, weekly age (0 year: 1-52, 1 year: 53-104, 2 year: 105-156)
  # For age >= 3 years, yearly age (3 year: 157, 100 year: 254)
  year2week = function(year) {
    
    # Convert to week value
    week = ifelse(
      test = year < 3, 
      yes  = round(year * 52), 
      no   = round((year - 2) + (3 * 52)))
    
    return(week)
  }
  
  # Repeat for each scenario of interest
  for (scenario in o$scenarios) {
    
    # Routine vaccine coverage
    routine_dt = fread(get_path("coverage")) %>%
      filter(country %in% o$countries, 
             vaccine != "SIA") %>%
      select(vaccine, country, year, coverage)
    
    # SIA vaccine coverage
    sia_dt = fread(get_path("coverage")) %>%
      filter(country %in% o$countries, 
             vaccine == "SIA", 
             coverage > 0) %>%
      # Convert age in years to weeks...
      mutate(a0 = year2week(age_first), 
             a1 = year2week(age_last)) %>%
      # Some fine adaptions...
      mutate(a0 = ifelse(age_first < 3, a0 + 1, a0), # Starting next weekly age
             a0 = pmax(a0, 1),      # Bound below by 1
             a1 = pmax(a1, 1)) %>%  # Bound below by 1
      select(-age_first, -age_last)
    
    # Save coverage to file for routine vaccination and SIAs
    fwrite(routine_dt, file = get_path("routine"))
    fwrite(sia_dt,     file = get_path("sia"))
  }
}

# ---------------------------------------------------------
# Prepare model parameters
# ---------------------------------------------------------
prepare_params = function(data) {
  
  # MODEL PARAMETERS:
  #   gamma := Recovery rate per timestep
  #   tstep := Number of timesteps per year
  #     amp := Amplification scale for seasonality
  #     ve1 := Vaccine efficacy for first dose by each age group
  # ve2plus := Vaccine efficacy for two and more doses
  #   vage2 := Age at second dose (country-specific)
  
  # Time and duration
  tstep	= 1000	# Number of timesteps in a year
  dinf  = 14		# Duration of infection (days)
  
  # Age-dependent vaccine efficacy for first dose, based on a linear model (Hughes et al. 2020)
  ve1_intcp = 0.64598  # Intercept of the linear model
  ve1_slope = 0.01485  # Slope of the linear model, per month of age
  ve2plus   = 0.98     # Vaccine efficacy for two and more doses
  
  # First dose efficacy by age
  age_ve1 = ve1_intcp + ve1_slope * 12 * c(1:(3*52)/52, 4:101)  # Based on age in months
  age_ve1 = ifelse(age_ve1 >= ve2plus, ve2plus, age_ve1)
  
  # Country-specific age at vaccination for MCV2
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
# Prepare data for this country and this scenario
# ---------------------------------------------------------
prepare_data = function(sim) {
  
  # ---- Coverage data ----
  
  # Coverage data for this scenario and country
  data = list(
    coverage_routine = load_coverage(sim, "routine"), 
    coverage_sia     = load_coverage(sim, "sia"))
  
  # ---- Load other input data for this country ----
  
  # All data rds files in inupt directory
  data_files = list.files(o$pth$input, pattern = "^data_.+\\.rds$")
  
  # Loop through files
  for (file in data_files) {
    
    # Shorthand reference for this data
    ref = str_remove_all(file, "(^data_|.rds$)")
    
    # Load data file
    this_data = readRDS(paste0(o$pth$input, file))
    
    # Filter datatable for this country
    if (is.data.frame(this_data))
      data[[ref]] = this_data[country == sim$country]
    
    # ... or select list element for this country
    if (!is.data.frame(this_data))
      data[[ref]] = this_data[[sim$country]]
  }
  
  # ---- Set basic reproduction number ----
  
  # If provided, use scenario-specific R0
  if (!is.na(sim$r0))
    data$r0 = sim$r0
  
  # If not provided but fixed, use user-defined R0
  if (is.na(sim$r0) && !is.na(o$fix_r0))
    data$r0 = o$fix_r0
  
  # If not provided nor fixed, use country-specific R0
  if (is.na(sim$r0) && is.na(o$fix_r0))
    data$r0 = data$rnought$r0
  
  # Remove now-redundant data
  data$rnought = NULL
  
  # ---- Other data considerations ----
  
  # Update timeliness if MCV1 not given at 39 weeks (9 months)
  if (data$vax_age$mcv1 != 39)
    data$timeliness[!is.na(age), timeliness := ifelse(age < data$vax_age$mcv1, 0, 1)]
  
  return(data)
}

# ---------------------------------------------------------
# Load coverage data
# ---------------------------------------------------------
load_coverage = function(sim, type) {
  
  # Construct file name and path
  file_name = paste1(type, sim$scenario)
  file_path = paste0(o$pth$coverage, file_name, ".csv")
  
  # Load coverage data and filter for years of interest
  coverage_data = fread(file_path) %>%
    filter(country %in% sim$country, 
           year %in% o$years)
  
  # Remove trivial values for SIA
  if (type == "sia") {
    coverage_data %<>%
      filter(coverage > 0)
  }
  
  return(coverage_data)
}

