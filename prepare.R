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
# Load data for this country and this scenario
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
    select(-country) %>%
    rename(country = country_code) %>%
    filter(country %in% sim$country, 
           year %in% o$years)
  
  # Remove trivial values for SIA
  if (type == "sia") {
    coverage_data %<>%
      filter(coverage > 0)
  }
  
  return(coverage_data)
}

# ---------------------------------------------------------
# Splits VIMC coverage files into routine (MCV1, MCV2) and SIAs
# ---------------------------------------------------------
create_coverage = function() {
  
  browser()
  
  # Repeat for each scenario of interest
  for (scenario in o$scenarios) {
    
    # Coverage file for this scenario - to be split between routine and SIA
    file = list(
      coverage = paste0(o$pth$coverage, "coverage_", scenario, ".csv"), 
      routine  = paste0(o$pth$coverage, "routine_",  scenario, ".csv"), 
      sia      = paste0(o$pth$coverage, "sia_",      scenario, ".csv"))
    
    # read vaccine coverage data file
    vaccov <- fread(file = file$coverage, na.strings = "<NA>")
    
    # select routine vaccination coverage
    keep_cols_routine <- c("vaccine", "country_code", "country", "year", "coverage")
    routine <- vaccov [activity_type != "campaign", ..keep_cols_routine]
    
    # select campaigns with coverage > 0 and with information of target population size
    keep_cols_sia <- c("vaccine", "country_code", "country", "year", "extent", "mid_day",
                       "age_first", "age_last", "age_range_verbatim", "target", "coverage", "coverage_subnat")
    sia <- vaccov [activity_type == "campaign" & ( !is.na(target) & !is.na(coverage) & coverage != 0),
                   ..keep_cols_sia]
    
    # check if national or subnational coverage >100%
    if (nrow(sia[coverage > 1 | coverage_subnat > 1]) > 0)
      stop ("national or subnational coverage > 100%")
    
    # -----------------------------------------------------------------------------------
    # calculate values for a0 and a1 to match the age groups
    sia [, `:=` (a0 = 0, a1 = 0)] #reached = round (as.numeric(target) * as.numeric(coverage))
    
    # for age < 3 years, weekly age (0 year: 1-52, 1 year: 53-104, 2 year: 105-156)
    # for age >= 3 years, yearly age (3 year: 157, 100 year: 254)
    find.a <- function (x) {
      t0 <- ifelse (x < 3, round (x * 52), round ((x - 2) + (3 * 52)) )
      return (t0) }
    
    sia [, `:=` (a0 = find.a(age_first), a1 = find.a(age_last))]
    sia [age_first < 3, a0 := a0 + 1] # starting the next weekly age
    
    # set age == 0 year as 1 week
    sia [a0 == 0, a0 := 1]
    sia [a1 == 0, a1 := 1]
    
    # write vaccine coverage data for routine vaccination and SIAs
    fwrite(x = routine, file = file$routine)
    fwrite(x = sia,     file = file$sia)
  }
}

