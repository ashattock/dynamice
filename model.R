###########################################################
# MODEL
#
# Main transmission model functionality. 
#
###########################################################

runScenario_rcpp <- function (
    scenario_name,
    vaccination,            # Whether children are vaccinated. 0: No vaccination; 1: Only MCV1; 2: MCV1 and MCV2
    using_sia) {              # Whether supplementary immunization campaigns are used. 0: no SIA; 1: with SIA  (Portnoy), 2: with SIA (7.7%)
  
  # ----------------------------------------------------------------------------
  # Run model
  # ----------------------------------------------------------------------------
  for (iso3 in o$countries) {
    out_run <- model(
      iso3          = iso3,
      scenario_name = scenario_name, 
      vaccination   = vaccination,
      using_sia     = using_sia)
    
    browser()
  }
}

# ------------------------------------------------------------------------------
# Run measles model for a given country and scenario
# ------------------------------------------------------------------------------
model = function(iso3, scenario_name, vaccination, using_sia) {
  
  # NOTES ON INPUTS: 
  #
  # vaccination := A numeric indicator that determines vaccination programmes
  # for children: 0 - No vaccination, 1 - Only MCV1,  2 - MCV1 and MCV2.
  #
  # using_sia := A numeric indicator that determines whether supplementary
  # immunisation activities (SIAs) are implemented and how SIAs are distributed
  # between zero-dose and already-vaccinated populations: 0 - no SIA, 1 - SIAs
  # based on a weighted logistic function fitted with Portnoy's data, and 2 -
  # SIAs based on an assumption that 7.7% of the population are never reached by
  # vaccination.
  
  # ---- Load data and parameters ----
  
  # Model input data
  d = prepare_data(iso3, scenario_name)
  
  # Global model parameters
  p = prepare_params(d)
  
  # ---- Model set up ----
  
  # country-specific timeliness curve
  country_timeliness <- d$timeliness[!is.na(age), timeliness]
  timeliness_ages    <- d$timeliness[!is.na(age), age]
  
  # expand 0-2 years old to weekly age strata
  s         <- 52 # number of finer stages within an age band (weekly ages, so 52)
  jt        <- 3  # how many ages to expand to s (or to weeks)
  beta_full <- matrix (0, ncol = 254, nrow = 254)
  
  beta_full[(1:(s*jt)), (1:(s*jt))] <- expandMatrix (
    A = d$contact [1:jt, 1:jt]/s,  # needs to be divided by 52 so that the mean total number of contacts stays the same
    expand_rows =  s, expand_cols =  s,
    rescale_rows = FALSE, rescale_cols = FALSE)
  
  beta_full[1:(s*jt),((s*jt)+1):(ncol(beta_full))] <- expandMatrix(
    A = d$contact [1:jt,(jt+1):ncol(d$contact)],
    expand_rows = s, expand_cols = 1,
    rescale_rows = F, rescale_cols = F)
  
  beta_full[((s*jt)+1):(nrow(beta_full)), 1:(s*jt)] <- expandMatrix(
    A = d$contact [(jt+1):nrow(d$contact),1:jt]/s,  # adjust to ensure the mean total number of contacts stays the same
    expand_rows = 1, expand_cols = s,
    rescale_rows = F, rescale_cols = F)
  
  beta_full[((s*jt)+1):(nrow(beta_full)), ((s*jt)+1):(ncol(beta_full))] <-
    d$contact [(jt+1):nrow(d$contact),(jt+1):ncol(d$contact)]
  
  beta_full_unadj <- beta_full
  
  # infection rate under target R0
  beta_tstep  <- d$rnought$r0 * p$gamma
  
  # setup inputs for DynaMICE Rcpp model
  n_years   <- length(o$years)
  
  y_out       <- array(0, c(254, 14, n_years))   # numbers of age groups, compartments, years
  case0d_out  <- array(0, c(254, n_years))
  case1d_out  <- array(0, c(254, n_years))
  case2d_out  <- array(0, c(254, n_years))
  pop_out     <- array(0, c(254, n_years))
  pop0d_out   <- array(0, c(254, n_years))
  popSus_out  <- array(0, c(254, n_years))
  dose_out    <- array(0, c(254, n_years))       # number of doses administrated
  reach0_out  <- array(0, c(254, n_years))       # number of zero-dose population reached
  fvp_out     <- array(0, c(254, n_years))       # number of fully vaccinated population reached
  
  init_Comp     <- matrix(0, 254, 14)
  init_Comp[,2] <- 0.95
  init_Comp[,3] <- 0.05
  
  t_spinup <- 1:1e5 # assume a fixed period for equilibrium period
  
  # ---- Main model loop ----
  
  message("  > ", iso3, ": running model")
  
  # run model by yearly input
  for (y in o$years) {
    
    message("   ~ ", y)
    
    # ---- Spin up ----
    
    #old script groups those aged 70-80, but division is by actual population size
    pop.vector <- d$population[year == y, value]
    
    # first expand polymod matrix (contact_tstep) and population vector and
    # then divide by population sizes, otherwise it doesn't work.
    pop.vector_full <- c(rep(pop.vector[1:jt]/s, each = s), pop.vector[(jt +1):length(pop.vector)])
    
    # change zero values to 1 to avoid division by 0
    pop.vector_full[pop.vector_full==0] <- 1
    
    # adjust contact reciprocity by population structure of each year
    beta_full <- matrix(0, 254, 254)
    beta_full_R0 <- matrix(0, 254, 254)
    for (i in 1:254){
      for (j in 1:254) {
        beta_full [i, j] <- (beta_full_unadj[i, j] * pop.vector_full[i] +
                               beta_full_unadj[j, i] * pop.vector_full[j])/(2*pop.vector_full[i])
        
        # calculate infection rate of the contact matrix based on the R0 definition
        # transform from "contactees per contactor" to "contactors per contactee"
        # This step can be skipped, since the largest eigenvalue remains unchanged.
        beta_full_R0 [i, j] <-  beta_full [i, j] * (pop.vector_full[i] / pop.vector_full[j])
      }
    }
    # make sure the contact matrix to represent target R0
    beta_full <- (beta_tstep / Re(eigen(beta_full_R0, only.values=T)$values[1])) * beta_full
    
    # run spin-up period
    if (y == o$years[1])
      out_Comp <- rcpp_spinup(init_Comp, p, beta_full, pop.vector_full, length(t_spinup))
    
    # ---- Coverage and timeliness ----
    
    if (vaccination >= 1) {
      
      # Maximum coverage can (obviously) only be 100%
      # To estimate proportion that is vaccinated at each week, we first calculate the number of individuals remaining susceptible
      # Then we calculate the number of individuals that should be vaccinated each week, in order to remain 1 - coverage susceptibles at the end of the timeliness data
      # In essence, this becomes the inverse of the cumulative timeliness curve
      cycov <- d$coverage_routine[year == y & vaccine == "MCV1", coverage]
      
      # Not use the following adjustment for coverage, as it leads to reduced the ...
      # final size of vaccinated population under the assumption of perfect timeliness
      # cycov <- d$coverage_routine [year == y & vaccine == "MCV1", coverage] /
      #   d$timeliness [is.na(age), prop_final_cov]
      
      if (length(cycov) == 0) {   # check if vaccine is not yet introduced and thereby, coverage value missing for this year
        cycov <- 0
      } else if (is.na(cycov)) {
        cycov <- 0
      }
      
      country_year_timeliness_mcv1 <- 1 - min(cycov,1) * country_timeliness
      
      country_year_timeliness_mcv1 <- -diff(country_year_timeliness_mcv1) /
        (country_year_timeliness_mcv1 [1:(length(country_year_timeliness_mcv1)-1)])
      
      # Timeliness is reported by week in the first year, and by month in the second year. Assume there is no vaccination in between
      country_year_timeliness_mcv1 [is.na(country_year_timeliness_mcv1)] <- 0
      country_year_timeliness_mcv1 [is.nan(country_year_timeliness_mcv1)] <- 0
      country_year_timeliness_mcv1_allages <- rep(0, 254)
      country_year_timeliness_mcv1_allages [round(timeliness_ages)] <- country_year_timeliness_mcv1
      
    } else {
      country_year_timeliness_mcv1_allages <- rep(0, 254)
    }
    
    if(vaccination == 2){
      country_year_mcv2 <- d$coverage_routine [year == y & vaccine == "MCV2", coverage]
    } else {
      country_year_mcv2 <- 0
    }
    
    # if ( length (country_year_mcv2) == 0 || is.na(country_year_mcv2) ) {
    #   country_year_mcv2 <- 0
    # }
    
    # set up SIA inputs
    cy_coverage_sia <- d$coverage_sia [year == y]
    if ((using_sia >= 1) && (dim(cy_coverage_sia)[1] > 0)) {
      # set up timesteps based on day of the year
      setorder(cy_coverage_sia, mid_day)
      sia_days <- cy_coverage_sia$mid_day
      t_sia_days <- c()
      for (iday in unique(sia_days)){
        t_sia_days <- c(t_sia_days,
                        round(p$tstep*(iday/365))-1 + (1:sum(sia_days == iday)))
      }
      
      sia_input <- list (
        sia_implement = using_sia, # choose SIA implementation approach
        a0 = cy_coverage_sia$a0,
        a1 = cy_coverage_sia$a1,
        siacov = as.double(cy_coverage_sia$coverage),
        siacov_subnat = as.double(cy_coverage_sia$ coverage_subnat),
        sia_subnat = as.integer(ifelse (cy_coverage_sia$extent %in% c("sub-national", "unknown"), 1, 0)), # distinguish subnational campaigns
        sia_tstep = as.integer(c(t_sia_days, 2000))  # add a larger-than-max-timestep number to stop loops in Rcpp
      )
    } else {
      sia_input <- list (
        sia_implement = as.integer(0),
        a0 = as.integer(0),
        a1 = as.integer(0),
        siacov = as.double(0),
        siacov_subnat = as.double(0),
        sia_subnat = as.integer(0),
        sia_tstep = as.integer(0)
      )
    }
    
    # ---- Simulate model ----
    
    t_start <- length(t_spinup) + (y-o$years[1])*p$tstep + 1
    
    outp <- rcpp_vaccine_oney(
      out_Comp,
      p,
      sia_input,
      beta_full,
      pop.vector_full,
      country_year_timeliness_mcv1_allages,
      country_year_mcv2,
      t_start)
    
    # ---- Store timestep output ----
    
    out_Comp <- outp$out_Comp
    
    # Store model output from this timestep
    y_out    [, , (y-o$years[1])+1] <- out_Comp
    case0d_out [, (y-o$years[1])+1] <- outp$cases_0d*pop.vector_full     # new cases among 0-dose
    case1d_out [, (y-o$years[1])+1] <- outp$cases_1d*pop.vector_full     # new cases among 1-dose
    case2d_out [, (y-o$years[1])+1] <- outp$cases_2d*pop.vector_full     # new cases among >=2-dose
    pop_out    [, (y-o$years[1])+1] <- rowSums(out_Comp[, 1:13])*pop.vector_full         # all compartments
    pop0d_out  [, (y-o$years[1])+1] <- rowSums(out_Comp[, 1:4])*pop.vector_full          # sum of M, S, I, R
    popSus_out [, (y-o$years[1])+1] <- rowSums(out_Comp[, c(2,5,8,11)])*pop.vector_full  # sum of S, V1S, V2S, V3S
    dose_out   [, (y-o$years[1])+1] <- outp$doses*pop.vector_full
    reach0_out [, (y-o$years[1])+1] <- outp$reach_d0*pop.vector_full
    fvp_out    [, (y-o$years[1])+1] <- outp$fvps*pop.vector_full
  }
  
  # Summarise complete model output
  output = list(
    cases0d  = rbind(colSums(case0d_out[1:52,]), colSums(case0d_out[53:104,]), colSums(case0d_out[105:156,]), case0d_out[157:254,]),
    cases1d  = rbind(colSums(case1d_out[1:52,]), colSums(case1d_out[53:104,]), colSums(case1d_out[105:156,]), case1d_out[157:254,]),
    cases2d  = rbind(colSums(case2d_out[1:52,]), colSums(case2d_out[53:104,]), colSums(case2d_out[105:156,]), case2d_out[157:254,]),
    pops     = rbind(colSums(   pop_out[1:52,]), colSums(   pop_out[53:104,]), colSums(   pop_out[105:156,]),    pop_out[157:254,]),
    pops0d   = rbind(colSums( pop0d_out[1:52,]), colSums( pop0d_out[53:104,]), colSums( pop0d_out[105:156,]),  pop0d_out[157:254,]),
    popsSus  = rbind(colSums(popSus_out[1:52,]), colSums(popSus_out[53:104,]), colSums(popSus_out[105:156,]), popSus_out[157:254,]),
    doses    = rbind(colSums(  dose_out[1:52,]), colSums(  dose_out[53:104,]), colSums(  dose_out[105:156,]),   dose_out[157:254,]),
    reachs0d = rbind(colSums(reach0_out[1:52,]), colSums(reach0_out[53:104,]), colSums(reach0_out[105:156,]), reach0_out[157:254,]),
    fvps     = rbind(colSums(   fvp_out[1:52,]), colSums(   fvp_out[53:104,]), colSums(   fvp_out[105:156,]),    fvp_out[157:254,]))
  
  save_name = paste1(iso3, scenario_name, vaccination, using_sia)
  save_file = paste0(o$pth$output, save_name, ".rds")
  
  # Save model output
  saveRDS(output, file = save_file)
}

# ------------------------------------------------------------------------------
# Create burden estimate csv files
# ------------------------------------------------------------------------------
get_burden_estimate <- function (
    scenario_name,
    save_scenario,
    log_name,
    vaccination,
    using_sia,
    folder_date,
    data
) {
  
  browser()
  
  # ---- Merge and process results ----
  
  # define folder name
  foldername <- paste0 (
    folder_date,
    "_v", vaccination,
    "_s", using_sia,
    "_deter")
  
  # merge results
  output_files <- list.files (path = paste0 ("outcome/", save_scenario, "/",foldername,"/"),
                              recursive = T, full.names = T)
  
  # Load VIMC results template 
  data_template = readRDS(paste0(o$pth$data, "data_template.rds"))
  
  # set up format of output file
  years          <- o$years
  ages           <- 0 : 100
  template    	 <- copy (setDT (data_template)) [year %in% years]   # drop rows if simulation period is shorter
  report_years   <- sort (unique (template$year))
  country_names  <- unique (subset (template, select = c("country", "country_name")))
  c_names        <- country_names$country_name
  names(c_names) <- country_names$country
  
  # file name for burden estimates
  burden_estimate_file <- paste0 ("central_burden_estimate_", scenario_name)
  
  browser() # We may want to load coverage for all countries here...
  
  # coverage file
  coverage_routine <- copy(fread(paste0(o$pth$coverage,
                                        "routine_",
                                        scenario_name,
                                        ".csv")))
  
  # read RDS files
  all_runs <- rbindlist (lapply (output_files, function (filename, ...) {
    res <- withCallingHandlers (
      readRDS (filename),
      warning = function(w) {warning(w, filename);}
    )
    filefinal <- stringr::str_extract (filename, "[^/]+$")      # remove path but keep filename
    res2 <- data.table (country     = rep (stringr::str_sub(filefinal, 1, 3), length(ages)*length(years)), # 101 ages * 121 years
                        year        = rep (years, each = length(ages)),
                        age         = rep (ages, length(years)),
                        #cases       = as.vector(res$cases),
                        cases0d     = as.vector(res$cases0d),
                        cases1d     = as.vector(res$cases1d),
                        cases2d     = as.vector(res$cases2d),
                        pops        = as.vector(res$pops),
                        pops0d      = as.vector(res$pops0d),
                        popsSus     = as.vector(res$popsSus),
                        doses       = as.vector(res$doses),
                        reachs0d    = as.vector(res$reachs0d),
                        fvps        = as.vector(res$fvps)
    )
    return (res2)
  }))
  
  # add country names and disease (Measles) to match template file
  all_runs[, c("country_name", "disease") := list (c_names[country], "Measles")]
  
  # select output years
  all_runs <- subset (all_runs, year %in% report_years)
  
  # ---- add columns for remaining life expectancy & MCV1 ----
  
  sel_countries <- unique (all_runs$country)
  
  # remaining life expectancy
  lexp_remain <- tailor_data_lexp_remain (sel_countries)
  
  
  all_runs <- lexp_remain [all_runs,
                           .(i.country, year, age, cases0d, cases1d, cases2d,
                             pops, pops0d, popsSus,  doses, reachs0d, fvps,
                             country_name, disease, value),
                           on = .(country_code = country,
                                  age          = age,
                                  year         = year)]
  
  # rename column names for output
  setnames (all_runs,
            old = c("i.country", "value"      ),
            new = c("country"  , "remain_lexp"))
  
  
  # MCV1 coverage
  coverage_routine_MCV1 <- coverage_routine [(vaccine == "MCV1") & (country_code %in% sel_countries)]
  
  all_runs <- coverage_routine_MCV1 [all_runs,
                                     .(i.country, i.year, age,
                                       cases0d, cases1d, cases2d,
                                       pops, pops0d, popsSus,
                                       doses, reachs0d, fvps, country_name,
                                       disease, coverage, remain_lexp),
                                     on = .(country_code = country,
                                            year         = year) ]
  
  
  # rename column "coverage" to "MCV1"
  setnames (x = all_runs,
            old = c("i.country", "i.year", "coverage"),
            new = c("country"  , "year"  , "MCV1"    ))
  
  browser() # Use data$cfr here
  
  # ---- calculate deaths and DALYs ----
  data_cfr_21 <- setDT (copy (data_cfr_portnoy_21))
  
  # The data contains CFR estimates between 1981 and 2020, for age between 0 to 99.
  # Extrapolate CFRs for year 1980 and 2021-2100 and for age 100.
  
  min_year = min (all_runs [, year])
  max_year = max (all_runs [, year])
  data_cfr_21 <- data_cfr_21 [year %in% min_year:max_year & country %in% sel_countries]
  
  if (min_year < 1981) {
    data_cfr_add <- rbindlist (lapply (min_year:1981, function(i) copy (data_cfr_21 [year == 1981, ])[, year := i]))
    data_cfr_21  <- rbind     (data_cfr_add, data_cfr_21, use.names = TRUE)
  }
  
  if (max_year > 2020) {
    data_cfr_add <- rbindlist (lapply (2021:max_year, function(i) copy (data_cfr_21 [year == 2020, ])[, year := i]))
    data_cfr_21  <- rbind     (data_cfr_21, data_cfr_add, use.names = TRUE)
  }
  
  data_cfr_21  <- rbind (data_cfr_21,
                         copy (data_cfr_21 [age == 99, ])[, age := 100],
                         use.names = TRUE)
  setorder (data_cfr_21, country, year, age)
  
  all_runs <- data_cfr_21 [all_runs,
                           .(disease, year, age, country, country_name,
                             pops, pops0d, popsSus, cases0d, cases1d, cases2d, doses,
                             reachs0d, fvps, cfr, remain_lexp),
                           on = .(country = country,
                                  year    = year,
                                  age     = age)]
  
  # estimate deaths
  all_runs [, `:=` (deaths0d = cases0d * cfr,
                    deaths1d = cases1d * cfr,
                    deaths2d = cases2d * cfr)]
  
  # calculate DALYs = (YLDs) + (YLLs)
  all_runs [, dalys := (((cases0d+cases1d+cases2d) - (deaths0d+deaths1d+deaths2d)) * 0.002) +
              ((deaths0d+deaths1d+deaths2d) * remain_lexp)]
  
  # adjust columns for output
  select.cols <- c("disease", "country", "country_name", "year", "age",
                   "pops", "pops0d", "popsSus",
                   "cases0d", "cases1d", "cases2d", "deaths0d", "deaths1d", "deaths2d",
                   "dalys", "doses", "reachs0d", "fvps")
  all_runs <- subset (all_runs, select = select.cols)
  
  # save burden estimates to file
  fwrite(x    = all_runs[order(country, year, age)],
         file = paste0(o$pth$central, burden_estimate_file, ".csv"))
  
  
  # release memory for the next round
  remove (list = c("all_runs", "sel_countries"))
  
  message("burden estimate csv file generated")
}
