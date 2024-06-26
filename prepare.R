###########################################################
# PREPARE
#
# Prepare for model simulation. Load model parameters and
# model input data.
#
###########################################################

# ---------------------------------------------------------
# Parent function for all preparation processes
# ---------------------------------------------------------
run_prepare = function() {
  
  # Only continue if specified by do_step
  if (!is.element(1, o$do_step)) return()
  
  message("* Preparing model resources")
  
  # Prepare coverage files
  prepare_coverage()
  
  # Prepare raw vaccine schedule data
  prepare_schedule()
}

# ---------------------------------------------------------
# Prepare coverage data by splitting into routine and SIA
# ---------------------------------------------------------
prepare_coverage = function() {
  
  # Function for constructing file paths
  get_path = function(x)
    paste0(o$pth$coverage, paste1(x, scenario), ".csv")
  
  # Function for converting age in years to week reference
  year2week = function(year) {
    
    # For age < 3 years, weekly age (0 year: 1-52, 1 year: 53-104, 2 year: 105-156)
    # For age >= 3 years, yearly age (3 year: 157, 100 year: 254)
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
# Prepare routine vaccine schedule data by country
# ---------------------------------------------------------
prepare_schedule = function() {
  
  # ---- Load and clean raw data ----
  
  # Raw data file from WHO
  #
  # SOURCE: immunizationdata.who.int/pages/schedule-by-disease/measles.html
  schedule_file = "vaccination_schedule.csv"
  schedule_path = paste0(o$pth$input, schedule_file)
  
  # Map unit references to scalers required for week format
  unit_dict = data.table(
    unit = c("Y", "M", "W"), 
    mult = c(52, 52/12, 1))
  
  # Regular expression for extracting year, month, week units
  unit_exp = paste(unit_dict$unit, collapse = ",")
  unit_exp = paste0("[", unit_exp, "]")
  
  # Load and clean data
  schedule_dt = fread(schedule_path) %>%
    # Retain only necessary columns...
    select(country  = ISO_3_CODE, 
           year     = YEAR, 
           schedule = SCHEDULEROUNDS, 
           target   = TARGETPOP, 
           age      = AGEADMINISTERED) %>%
    # Retain general routine doses for countries of interest...
    filter(country  %in% o$countries,
           schedule %in% c(1, 2),
           target == "") %>%
    # Parse age string into value and unit...
    mutate(unit  = str_extract(age, unit_exp), 
           value = str_extract(age, "[0-9]+")) %>%
    replace_na(list(unit = mode(unit))) %>%
    # Append multiplier to convert all ages to weeks...
    left_join(y  = unit_dict, 
              by = "unit") %>%
    mutate(weeks = mult * as.numeric(value)) %>%
    # Take the mean where we have multiple entires...
    group_by(country, schedule) %>%
    summarise(age = round(mean(weeks))) %>%
    ungroup() %>%
    # If over 3 years old, convert to annual increments...
    mutate(age = ifelse(
      test = age > 52 * 3, 
      yes  = 52 * 3 + floor(age / 52 - 2), 
      no   = age)) %>%
    # Convert to wide format...
    mutate(schedule = paste0("mcv", schedule)) %>%
    pivot_wider(names_from  = schedule, 
                values_from = age) %>%
    # Tidy up...
    arrange(country) %>%
    as.data.table()
  
  # ---- Impute any missing values with regional average ----
  
  # Only needed to any data missing
  if (any(is.na(schedule_dt))) {
    
    # Load country-region details
    region_dt = fread(paste0(o$pth$config, "regions.csv"))
    
    # Take the mean from each region
    regional_mean = schedule_dt %>%
      left_join(y  = region_dt, 
                by = "country") %>%
      group_by(region) %>%
      summarise(mcv1_mean = round(mean(mcv1, na.rm = TRUE)), 
                mcv2_mean = round(mean(mcv2, na.rm = TRUE))) %>%
      ungroup() %>%
      as.data.table()
    
    # Impute missing values with regional mean
    schedule_dt %<>%    
      left_join(y  = region_dt, 
                by = "country") %>%  
      left_join(y  = regional_mean, 
                by = "region") %>%
      mutate(mcv1 = ifelse(is.na(mcv1), mcv1_mean, mcv1), 
             mcv2 = ifelse(is.na(mcv2), mcv2_mean, mcv2)) %>%
      select(country, mcv1, mcv2)
  }
  
  # ---- Write (or overwrite data file) ----
  
  # File path and name to save to
  save_file = paste0(o$pth$input, "data_vax_age.rds")
  
  # Save the file in input dir
  saveRDS(schedule_dt, file = save_file)
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
  ve1_intcp = 0.70  # Intercept of the linear model
  ve1_slope = 0.02  # Slope of the linear model, per month of age
  ve2plus   = 0.98  # Vaccine efficacy for two and more doses
  
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
    all_data = readRDS(paste0(o$pth$input, file))
    
    # Filter datatable for this country
    if (is.data.frame(all_data))
      this_data = all_data[country == sim$country]
    
    # ... or select list element for this country
    if (!is.data.frame(all_data))
      this_data = all_data[[sim$country]]
    
    # Impute values from regional countries if missing
    if (is.null(this_data) || nrow(this_data) == 0)
      this_data = impute_missing_data(ref, all_data, sim$country)
    
    # Store data
    data[[ref]] = this_data
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

# ---------------------------------------------------------
# Impute values from regional countries if missing
# ---------------------------------------------------------
impute_missing_data = function(ref, all_data, country) {
  
  # Load country-region details
  region_dt = fread(paste0(o$pth$config, "regions.csv"))
  
  # Region of country in question
  region = region_dt %>%
    filter(country == !!country) %>% 
    pull(region)
  
  # ---- CFR ----
  
  # Check data reference
  if (ref == "cfr_portnoy_21") {
    
    # Mean of countries in region
    this_data = all_data %>%
      left_join(y  = region_dt, 
                by = "country") %>%
      filter(region == !!region) %>%
      group_by(year, age) %>%
      summarise(cfr = mean(cfr)) %>%
      ungroup() %>%
      mutate(country = country) %>%
      select(all_of(names(all_data))) %>%
      as.data.table()
  }
  
  # ---- R0 ----
  
  # Check data reference
  if (ref == "rnought") {
    
    # Mean of countries in region
    this_data = all_data %>%
      left_join(y  = region_dt, 
                by = "country") %>%
      group_by(region) %>%
      summarise(r0 = mean(r0)) %>%
      ungroup() %>%
      filter(region == !!region) %>%
      mutate(country = country) %>%
      select(all_of(names(all_data))) %>%
      as.data.table()
  }
  
  # ---- Timeliness ----
  
  # Check data reference
  if (ref == "timeliness") {
    
    # Assume prompt timeliness
    this_data = all_data %>%
      select(age, timeliness) %>%
      unique() %>%
      mutate(country = country, 
             .before = 1) %>%
      mutate(prop_final_cov = 
               ifelse(is.na(age), 1, NA))
  }
  
  # ---- Life expectancy ----
  
  # Check data reference
  if (ref == "life_exp") {
    
    # Mean of countries in region
    this_data = all_data %>%
      left_join(y  = region_dt, 
                by = "country") %>%
      filter(region == !!region) %>%
      group_by(year) %>%
      summarise(value = mean(value)) %>%
      ungroup() %>%
      mutate(country = country) %>%
      select(all_of(names(all_data))) %>%
      as.data.table()
  }
  
  # ---- Population ----
  
  # Check data reference
  if (ref == "population") {
    
    # Assume trivial
    this_data = all_data %>%
      mutate(country = !!country, 
             value   = 0) %>%
      unique()
  }
  
  # ---- Contact rate ----
  
  # Check data reference
  if (ref == "contact") {
    
    # Regional neighbours to take mean over
    neighbours = region_dt %>%
      filter(region == !!region) %>% 
      pull(country) %>%
      intersect(names(all_data))
    
    # Extract data from regional neighbours
    data_list = all_data[neighbours]
    data_dims = dim(data_list[[1]])
    
    # Combine into 3D array
    #
    # See: stackoverflow.com/questions/26018216/calculating-mean-of-multiple-matrices-in-r
    data_array = array(
      data = do.call(cbind, data_list), 
      dim  = c(data_dims, length(neighbours)))
    
    # Take the mean to derive country values
    this_data = apply(data_array, c(1, 2), mean)
  }
  
  return(this_data)
}

