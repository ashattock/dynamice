# run_scenario
# execute simulations for MR-MAPs scenarios
# update: 2023/04/07

run_scenarios = function() {
  
  # Only continue if specified by do_step
  if (!is.element(1, o$do_step)) return()
  
  message("* Running vaccine scenarios")
  
  # ---- Load data ----
  
  message (" - Loading data")
  
  # Load all data
  data = load_data()
  
  # add data for country-specific age at vaccination
  # https://immunizationdata.who.int/pages/schedule-by-disease/measles.html
  # input monthly age and then covert to weekly age for dynaMICE structure
  data_vage <- data.table(
    country = o$countries,
    mcv1 = ceiling (c(10.5,  9, 9, 9, 8,
                      9,  9, 9, 9, 9,
                      9, 12, 9, 9)/12*52),
    mcv2 = ceiling (c(20,   15,   18, 15, 18,
                      13.5, 16.5, 16.5, 15, 15,
                      16.5,   72,   15, 16.5)/12*52))
  # adjust to yearly age for those >= 3 years old
  data_vage[mcv2 > 52*3, mcv2 := 52*3 + floor(mcv2/52-2)]
  
  # update timeliness data for countries not giving MCV1 to 39 week old (9 months)
  for (ictry in data_vage[mcv1 != 39, country]) {
    c_age <- data_vage[country == ictry, mcv1]
    data$timeliness[country == ictry & !is.na(age), timeliness := ifelse (age < c_age, 0, 1)]
  }
  
  # assume a fixed R0 for the central run
  # median R0 for least developed countries, vaccine era, from Guerra et al. (2017)
  if (!is.na(o$fix_r0))
    data$r0[country %in% o$countries, r0 := o$fix_r0]
  
  # ---- Configure scenarios ----
  
  message (" - Configuring scenarios")
  
  # set up variables
  var <- list(
    log_name                   = "test_log") 
  
  # prepare coverage inputs
  vac_strategies <- c(
    "nomcv",                 # (1) no vaccination
    "mcv1",                  # (2) MCV1 only
    "mcv1-mcv2",             # (3) MCV1 + MCV2
    "mcv1-mcv2-sia",         # (4) MCV1 + MCV2 + SIA
    "mcv1-sia",              # (5) MCV1 + SIA
    "mcv1-mcv2alt1",         # (6) MCV1 + MCV2(early intro, fast rollout)
    "mcv1-mcv2alt1-sia",     # (7) MCV1 + MCV2(early intro, fast rollout) + SIA
    "mcv1-mcv2-siaalt1",     # (8) MCV1 + MCV2 + SIA(zero dose first)
    "mcv1-mcv2-siaalt2",     # (9) MCV1 + MCV2 + SIA(already vaccinated first)
    "mcv1-siaalt1",          # (10) MCV1 + SIA(zero dose first)
    "mcv1-siaalt2",          # (11) MCV1 + SIA(already vaccinated first)
    "mcv1-mcv2alt1-siaalt1", # (12) MCV1 + MCV2(early intro) + SIA(zero dose first)
    "mcv1-mcv2alt2")         # (13) MCV1 + MCV2(early intro, gradual rollout)
  
  # set SIAs implementation method for each scenario
  # additional modifications needed in Rcpp for extent-specific approaches
  # ex. random reach for subnational campaign & 7.7% at national level for national and rollover-nat campaigns
  # 0: no SIA
  # 1: random reach (baseline assumption)
  # 2: 7.7% less likely to be reached at national level
  # 3: zero-dose first
  # 4: already-vaccinated first
  # 5: 7.7% less likely to be reached at subnational level
  set_sia         <- c (0, 0, 0, 2, 2, 0, 2, 3, 4, 3, 4, 3, 0)
  
  # set routine vaccination parameters to distinguish between scenarios
  # 0: no routine MCV
  # 1: MCV1 only
  # 2: MCV1 + MCV2
  set_vaccination <- c (0, 1, 2, 2, 1, 2, 2, 2, 2, 1, 1, 2, 2)
  
  # ---- Generate vaccine coverage scenarios ----
  
  # Prepare coverage input data - update when the data are changed
  if (o$reload_coverage)
    lapply(vac_strategies, create_vaccine_coverage_routine_sia)
  
  # ---- Run model for SIA activities ----
  
  # run model by different SIA methods
  for (isia in c(1,2,5)){
    
    if (isia == 2){
      sel_scns <- 1:length(vac_strategies) # main assumption
      set_sia  <- c (0, 0, 0, 2, 2, 0, 2, 3, 4, 3, 4, 3, 0)
    } else {
      sel_scns <- c(2,3,4,5,7) # evaluate different SIA assumptions
      set_sia [c(4,5,7,14)] <- isia
    }
    
    for (index in sel_scns){
      
      scenario_name  <- vac_strategies [index]
      message("  > ", scenario_name)
      
      scenario_number <- sprintf("scenario%02d", index)
      
      # run model and estimate cases
      burden_estimate_file <- runScenario_rcpp (
        scenario_name          = scenario_name,
        save_scenario          = scenario_number,
        log_name               = var$log_name,
        vaccination            = set_vaccination [index],
        using_sia              = set_sia         [index],
        sim_years              = 1980:2020, 
        data                   = data)
      
      browser()
      
      # separately estimate dalys
      burden_estimate_file <- paste0 ("central_burden_estimate_",
                                      scenario_name, ".csv")
      
      # merge outputs into csv files
      get_burden_estimate(
        # vaccine_coverage_subfolder = var$vaccine_coverage_subfolder,
        scenario_name              = vac_strategies [index],
        save_scenario              = scenario_number,
        log_name                   = var$log_name,
        vaccination                = set_vaccination [index],
        using_sia                  = set_sia         [index],
        folder_date                = "20230401",
        sim_years                  = 1980:2020, 
        data                       = data)
    }
    
    browser()
    
    # move files to a specified folder
    res_files <- list.files(o$pth$central)
    dir.create (paste0 ("previous_res/20230401/siareach_", isia, "/"))
    file.rename (from = paste0 (o$pth$central, res_files),
                 to = paste0 ("previous_res/20230401/siareach_", isia, "/", res_files))
  }
  
  # analyse burden estimates under different R0
  set_sia <- c (0, 0, 0, 2, 2, 0, 2, 3, 4, 3, 4, 3, 0) # set back to the oringinal assumptions
  
  for (ir0 in c(seq(6,26,2))) {
    # vary R0 values
    data_r0 <- copy (data_r0) [ , r0 := ir0]
    source ("R/functions_rcpp.R")
    
    if (ir0 %in% c(6,16,26)){
      sel_scns <- c(2,3,4,5)
    } else {
      sel_scns <- 2
    }
    
    for (index in sel_scns){
      
      scenario_name   <- vac_strategies [index]
      scenario_number <- sprintf ("scenario%02d", index)
      
      # run model and estimate cases
      burden_estimate_file <- runScenario_rcpp (
        scenario_name              = scenario_name,
        save_scenario              = scenario_number,
        log_name                   = var$log_name,
        vaccination                = set_vaccination [index],
        using_sia                  = set_sia [index],
        sim_years                  = 1980:2020,
        data                       = data)
      
      # separately estimate deaths
      burden_estimate_file <- paste0 ("central_burden_estimate_",
                                      scenario_name, ".csv")
      
      # merge outputs into csv files
      get_burden_estimate(
        scenario_name              = vac_strategies [index],
        save_scenario              = sprintf ("scenario%02d", index),
        log_name                   = var$log_name,
        vaccination                = set_vaccination [index],
        using_sia                  = set_sia         [index],
        folder_date                = "20230401",
        sim_years                  = 1980:2020, 
        data                       = data)
      
      # rename file
      file.rename (from = paste0(o$pth$central, burden_estimate_file),
                   to = paste0(o$pth$central,
                               "central_burden_estimate_",
                               scenario_name,
                               "_r0-", ir0, ".csv"))
    }
  }
  
  browser()
  
  # move files to a specified folder
  res_files <- list.files(o$pth$central)
  dir.create (paste0 ("previous_res/20230401/siareach_2/senanl_r0"))
  file.rename (from = paste0(o$pth$central, res_files),
               to = paste0("previous_res/20230401/siareach_2/senanl_r0/", res_files))
}

# ---------------------------------------------------------
# Load all data
# ---------------------------------------------------------
load_data = function() {
  
  # TODO: Better to do this for one country at a time?
  
  # Initiate data list
  data = list()
  
  # All data rds files in data directory
  all_files = list.files(o$pth$data, pattern = ".+\\.rds$")
  
  # Loop through files
  for (file in all_files) {
    
    # Shorthand reference for this data
    ref = str_remove_all(file, "(^data_|.rds$)")
    
    # Load data and store in list
    data[[ref]] = readRDS(paste0(o$pth$data, file))
  }
  
  return(data)
}

