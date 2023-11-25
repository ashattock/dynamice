# functions_rcpp.R
# main functions for running the DynaMICE model based on Rcpp files
# update: 2022/01/20

# ------------------------------------------------------------------------------
#' Run the Rcpp measles model for a selected vaccination strategy
#'
#' A function that executes the Rcpp measles model under a selected vaccination
#' scenario, including a pre-specified set of countries.
# ------------------------------------------------------------------------------
#' @param scenario_name Name of the vaccination scenario selected or being
#' analysed.
#' @param save_scenario A folder name for saving results from a selected
#' scenario, denoted by a two-digit number. e.g. "scenario08".
#' @param log_name A file name for keeping a log.
#' @param vaccination A numeric indicator that determines vaccination programmes
#'  for children: 0 - No vaccination, 1 - Only MCV1, and 2 - MCV1 and MCV2.
#' @param using_sia A numeric indicator that determines whether supplementary
#' immunisation activities (SIAs) are implemented and how SIAs are distributed
#' between zero-dose and already-vaccinated populations: 0 - no SIA, 1 - SIAs
#' based on a weighted logistic function fitted with Portnoy's data, and 2 -
#' SIAs based on an assumption that 7.7% of the population are never reached by
#' vaccination.
#' @param sim_years A numeric vector containing calendar years included for model
#'  simulation.
#'
#' @importFrom foreach %dopar% %:% foreach
#' @import data.table
#'
#' @examples
#' \dontrun{
#' runScenario_rcpp (
#'   scenario_name              = "campaign-only-default",
#'   save_scenario               = scenario_number,
#'   log_name                   = "test_log",
#'   vaccination                = 0,
#'   using_sia                  = 1
#'   sim_years                  = 1980:2100
#'   )
#'   }
runScenario_rcpp <- function (
    scenario_name,
    save_scenario,
    log_name,
    vaccination,            # Whether children are vaccinated. 0: No vaccination; 1: Only MCV1; 2: MCV1 and MCV2
    using_sia,              # Whether supplementary immunization campaigns are used. 0: no SIA; 1: with SIA  (Portnoy), 2: with SIA (7.7%)
    sim_years,              # calendar years for simulation
    data
) {
  
  # --------------------------------------------------------------------------
  # define global parameters
  # --------------------------------------------------------------------------
  ages 		    <- c (0:100)	              # Numeric vector with age-strata that are reported (Model ALWAYS models 101 age-groups; 0-2yo in weekly age-strata, and 3-100yo in annual age-strata)
  dinf		    <- 14				                # duration of infection (days)
  tstep			  <- 1000				              # number of timesteps in a year
  
  # age-dependent vaccine efficacy for first dose, based on a linear model (Hughes et al. 2020)
  ve1_intcp   <- 0.64598                  # intercept of the linear model
  ve1_slope   <- 0.01485                  # slope of the linear model, per month of age
  ve2plus     <- 0.98                     # vaccine efficacy for two and more doses
  
  age_ve1 <- ve1_intcp + ve1_slope * 12 * c(1:(3*52)/52, 4:101)  # based on age in months
  age_ve1 <- ifelse (age_ve1 >= ve2plus, ve2plus, age_ve1)
  
  # parameters for Rcpp functions
  parms_rcpp <- list (gamma         = 1 / (dinf * tstep/365),    # rate of losing infectivity
                      tstep         = tstep,
                      amp           = 0.05,		                   # amplitude for seasonality
                      ve1           = age_ve1,
                      ve2plus       = ve2plus)
  
  browser()
  
  # --------------------------------------------------------------------------
  # prepare input data
  # --------------------------------------------------------------------------
  # define filename of coverage data
  data_coverage_routine <- paste0(o$pth$coverage,
                                  "routine_",
                                  scenario_name,
                                  ".csv")
  
  data_coverage_sia <- paste0 (o$pth$coverage,
                               "sia_",
                               scenario_name,
                               ".csv")
  
  browser()
  
  # import/read data
  coverage_routine	<- copy(fread(data_coverage_routine))[year %in% sim_years]
  coverage_sia		  <- copy(fread(data_coverage_sia))[year %in% sim_years]
  timeliness  		  <- setDT(data_timeliness)
  rnought	    		  <- setDT(data_r0)
  population  		  <- setDT(data_pop)
  
  # use synthetic contact matrices
  contact_list <- sapply (o$countries,
                          function(cty){data_contact_syn[[cty]]},
                          simplify = FALSE, USE.NAMES = TRUE)
  
  browser()
  
  # ----------------------------------------------------------------------------
  # Run model
  # ----------------------------------------------------------------------------
  for (iso3 in o$countries) {
    out_run <- runCountry_rcpp (iso3               = iso3,
                                years              = as.numeric (sim_years),
                                vaccination        = vaccination,
                                using_sia          = using_sia,
                                parms_rcpp         = parms_rcpp,
                                c_coverage_routine = coverage_routine[country_code == iso3,],
                                c_coverage_sia     = coverage_sia[country_code == iso3 & coverage != 0,],
                                c_timeliness       = timeliness[country_code == iso3,],
                                c_contact          = contact_list[[iso3]],
                                c_rnought          = rnought[country_code == iso3, r0],
                                c_population       = population[country_code == iso3,],
                                save_scenario      = save_scenario,
                                log_name           = log_name)
  }
}

# ------------------------------------------------------------------------------
#' Execute the Rcpp measles model for a single country run
#'
#' A function nested under \code{\link{runScenario_rcpp}} to run Rcpp codes
#'  for measles vaccination, given a particular country and a variable set of
#'  probabilistic sensitivity analysis (PSA).
# ------------------------------------------------------------------------------
#' @param iso3 ISO-3 code of the selected country.
#' @param years A vector containing continuous calender years for simulation.
#' @param vaccination A numeric indicator that determines vaccination programmes
#'  for children: 0 - No vaccination, 1 - Only MCV1,  2 - MCV1 and MCV2.
#' @param using_sia A numeric indicator that determines whether supplementary
#' immunisation activities (SIAs) are implemented and how SIAs are distributed
#' between zero-dose and already-vaccinated populations: 0 - no SIA, 1 - SIAs
#' based on a weighted logistic function fitted with Portnoy's data, and 2 -
#' SIAs based on an assumption that 7.7% of the population are never reached by
#' vaccination.
#' @param parms_rcpp A list includes time-invariant parameters for Rcpp functions:
#' recovery rate per timestep (\code{gamma}), number of timesteps per year (\code{tstep}),
#'  amplification scale for seasonality (\code{amp}), vaccine efficacy for first dose
#' by each age group (\code{ve1}), vaccine efficacy for two and more doses (\code{ve2plus}).
#'  assuming 'take' ('all-or-not') protection for each age group (weekly age
#' groups for age 0, 1, and 2; yearly age groups for age between 3 and 100).
#' @param c_coverage_routine A data frame for routine vaccination coverage under
#' a selected scenario for a specific country.
#' @param c_coverage_sia A data frame for SIA coverage under a selected scenario
#'  for a specific country.
#' @param c_timeliness A data frame for timeliness estimates by age for a
#' specific country.
#' @param c_contact A data frame for the contact matrix by age (0-100 years old)
#'  for a specific country.
#' @param c_rnought A numeric variable of R0 value for the selected country
#' \code{iso3}.
#' @param c_population A data frame for population size by age for a specific
#' country.
#' @param save_scenario A folder name for saving results from a selected
#' scenario, denoted by a two-digit number. e.g. "scenario08".
#' @param log_name A file name for keeping a log.
#'
#' @importFrom Rcpp sourceCpp
#' @import data.table
#'
#' @examples
#' \dontrun{
#' runCountry_rcpp (
#'   iso3               = "BGD",
#'   years              = 1980:2020,
#'   vaccination        = 1,
#'   using_sia          = 1,
#'   parms_rcpp         = list (gamma = 1/(14*1000/365),
#'                              tstep = 1000,
#'                              amp   = 0.05,
#'                              ve1   = c(rep(0,40), rep(0.7, 214)),
#'                              ve2plus = 0.98),
#'   c_coverage_routine = coverage_routine[country_code == "BGD"],
#'   c_coverage_sia     = coverage_sia[country_code == "BGD" & coverage!=0],
#'   c_timeliness       = timeliness[country_code == "BGD"],
#'   c_contact          = contact[["BGD"]],
#'   c_rnought          = 10,
#'   c_population       = population[country_code == "BGD"],
#'   save_scenario      = "scenario01",
#'   foldername         = NULL,
#'   log_name           = "test_log"
#'   )
#'   }
runCountry_rcpp <- function (
    #variables specific for loop
  iso3,
  years,
  
  #infection dynamic variables
  vaccination,
  using_sia,
  parms_rcpp,
  
  #input data
  c_coverage_routine,
  c_coverage_sia,
  c_timeliness,
  c_contact,
  c_rnought,
  c_population,
  
  # dynaMICE model options
  save_scenario,
  foldername,
  log_name
) {
  
  
  
  
  
  
  
  
  # country-specific timeliness curve
  country_timeliness <- c_timeliness [!is.na(age), timeliness]
  timeliness_ages    <- c_timeliness [!is.na(age), age]
  
  # country-specific age at vaccination for MCV2
  parms_rcpp$vage2 <- as.integer (data_vage [country_code == iso3, "mcv2"])
  
  # expand 0-2 years old to weekly age strata
  s         <- 52 # number of finer stages within an age band (weekly ages, so 52)
  jt        <- 3  # how many ages to expand to s (or to weeks)
  beta_full <- matrix (0, ncol = 254, nrow = 254)
  
  beta_full[(1:(s*jt)), (1:(s*jt))] <- expandMatrix (
    A = c_contact [1:jt, 1:jt]/s,  # needs to be divided by 52 so that the mean total number of contacts stays the same
    expand_rows =  s, expand_cols =  s,
    rescale_rows = FALSE, rescale_cols = FALSE)
  
  beta_full[1:(s*jt),((s*jt)+1):(ncol(beta_full))] <- expandMatrix(
    A = c_contact [1:jt,(jt+1):ncol(c_contact)],
    expand_rows = s, expand_cols = 1,
    rescale_rows = F, rescale_cols = F)
  
  beta_full[((s*jt)+1):(nrow(beta_full)), 1:(s*jt)] <- expandMatrix(
    A = c_contact [(jt+1):nrow(c_contact),1:jt]/s,  # adjust to ensure the mean total number of contacts stays the same
    expand_rows = 1, expand_cols = s,
    rescale_rows = F, rescale_cols = F)
  
  beta_full[((s*jt)+1):(nrow(beta_full)), ((s*jt)+1):(ncol(beta_full))] <-
    c_contact [(jt+1):nrow(c_contact),(jt+1):ncol(c_contact)]
  
  beta_full_unadj <- beta_full
  
  # infection rate under target R0
  beta_tstep  <- c_rnought*parms_rcpp$gamma
  
  # setup inputs for DynaMICE Rcpp model
  n_years   <- length(years)
  
  y_out       <- array(0, c(254, 14, n_years))   # numbers of age groups, compartments, years
  #case_out    <- array(0, c(254, n_years))
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
  
  writelog (log_name, paste0(iso3, "; Generated input data & started model run"))
  
  # run model by yearly input
  for (y in years) {
    
    #old script groups those aged 70-80, but division is by actual population size
    pop.vector <- c_population[year == y, value]
    
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
    if (y == years[1]){
      out_Comp <- rcpp_spinup (init_Comp, parms_rcpp, beta_full, pop.vector_full, length(t_spinup))
      # print ('Spin-up period finished')
    }
    
    if (vaccination >= 1) {
      
      # Maximum coverage can (obviously) only be 100%
      # To estimate proportion that is vaccinated at each week, we first calculate the number of individuals remaining susceptible
      # Then we calculate the number of individuals that should be vaccinated each week, in order to remain 1 - coverage susceptibles at the end of the timeliness data
      # In essence, this becomes the inverse of the cumulative timeliness curve
      cycov <- c_coverage_routine [year == y & vaccine == "MCV1", coverage]
      
      # Not use the following adjustment for coverage, as it leads to reduced the ...
      # final size of vaccinated population under the assumption of perfect timeliness
      # cycov <- c_coverage_routine [year == y & vaccine == "MCV1", coverage] /
      #   c_timeliness [is.na(age), prop_final_cov]
      
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
      country_year_mcv2 <- c_coverage_routine [year == y & vaccine == "MCV2", coverage]
    } else {
      country_year_mcv2 <- 0
    }
    
    # if ( length (country_year_mcv2) == 0 || is.na(country_year_mcv2) ) {
    #   country_year_mcv2 <- 0
    # }
    
    
    # set up SIA inputs
    cy_coverage_sia <- c_coverage_sia [year == y]
    if ((using_sia >= 1) && (dim(cy_coverage_sia)[1] > 0)) {
      # set up timesteps based on day of the year
      setorder(cy_coverage_sia, mid_day)
      sia_days <- cy_coverage_sia$mid_day
      t_sia_days <- c()
      for (iday in unique(sia_days)){
        t_sia_days <- c(t_sia_days,
                        round(parms_rcpp$tstep*(iday/365))-1 + (1:sum(sia_days == iday)))
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
    
    t_start <- length(t_spinup) + (y-years[1])*parms_rcpp$tstep + 1
    outp <- rcpp_vaccine_oney (out_Comp,
                               parms_rcpp,
                               sia_input,
                               beta_full,
                               pop.vector_full,
                               country_year_timeliness_mcv1_allages,
                               country_year_mcv2,
                               t_start)
    out_Comp <- outp$out_Comp
    
    y_out    [, , (y-years[1])+1] <- out_Comp
    #case_out   [, (y-years[1])+1] <- outp$cases*pop.vector_full
    case0d_out [, (y-years[1])+1] <- outp$cases_0d*pop.vector_full     # new cases among 0-dose
    case1d_out [, (y-years[1])+1] <- outp$cases_1d*pop.vector_full     # new cases among 1-dose
    case2d_out [, (y-years[1])+1] <- outp$cases_2d*pop.vector_full     # new cases among >=2-dose
    pop_out    [, (y-years[1])+1] <- rowSums(out_Comp[, 1:13])*pop.vector_full         # all compartments
    pop0d_out  [, (y-years[1])+1] <- rowSums(out_Comp[, 1:4])*pop.vector_full          # sum of M, S, I, R
    popSus_out [, (y-years[1])+1] <- rowSums(out_Comp[, c(2,5,8,11)])*pop.vector_full  # sum of S, V1S, V2S, V3S
    dose_out   [, (y-years[1])+1] <- outp$doses*pop.vector_full
    reach0_out [, (y-years[1])+1] <- outp$reach_d0*pop.vector_full
    fvp_out    [, (y-years[1])+1] <- outp$fvps*pop.vector_full
    
    # if(y %% 20 == 0) 
    #   message('year ', y, ' finished')
  }
  
  save_name = paste1(save_scenario, iso3, vaccination, using_sia)
  save_file = paste0(pth$output, save_name, ".rds")
  
  saveRDS (list(#cases   = rbind (colSums(  case_out[1:52,]), colSums(  case_out[53:104,]), colSums(  case_out[105:156,]),
    cases0d  = rbind (colSums(case0d_out[1:52,]), colSums(case0d_out[53:104,]), colSums(case0d_out[105:156,]), case0d_out[157:254,]),
    cases1d  = rbind (colSums(case1d_out[1:52,]), colSums(case1d_out[53:104,]), colSums(case1d_out[105:156,]), case1d_out[157:254,]),
    cases2d  = rbind (colSums(case2d_out[1:52,]), colSums(case2d_out[53:104,]), colSums(case2d_out[105:156,]), case2d_out[157:254,]),
    pops     = rbind (colSums(   pop_out[1:52,]), colSums(   pop_out[53:104,]), colSums(   pop_out[105:156,]),    pop_out[157:254,]),
    pops0d   = rbind (colSums( pop0d_out[1:52,]), colSums( pop0d_out[53:104,]), colSums( pop0d_out[105:156,]),  pop0d_out[157:254,]),
    popsSus  = rbind (colSums(popSus_out[1:52,]), colSums(popSus_out[53:104,]), colSums(popSus_out[105:156,]), popSus_out[157:254,]),
    doses    = rbind (colSums(  dose_out[1:52,]), colSums(  dose_out[53:104,]), colSums(  dose_out[105:156,]),   dose_out[157:254,]),
    reachs0d = rbind (colSums(reach0_out[1:52,]), colSums(reach0_out[53:104,]), colSums(reach0_out[105:156,]), reach0_out[157:254,]),
    fvps     = rbind (colSums(   fvp_out[1:52,]), colSums(   fvp_out[53:104,]), colSums(   fvp_out[105:156,]),    fvp_out[157:254,])),
    file = save_file)
  
  writelog (log_name, paste0 (iso3, "; Finished model run & saved outputs"))
}
# end of function -- runCountry_rcpp
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
#' Get burden estimate csv files
#'
#' A function that combines all the RDS files for a selected scenario into csv
#' files. Each RDS file contains case estimates of a single country from a
#' single run using the function \code{runScenario_rcpp()}. If stochastic runs
#' are included (\code{psa > 0}), an additional csv is generated for the mean
#' estimates of each country.
# ------------------------------------------------------------------------------
#' @param scenario_name Name of the vaccination scenario selected or being
#' analysed.
#' @param save_scenario A folder name for saving results from a selected
#' scenario, denoted by a two-digit number. e.g. "scenario08".
#' @param log_name A file name for keeping a log.
#' @param vaccination A numeric indicator that determines vaccination programmes
#'  for children: 0 - No vaccination, 1 - Only MCV1, and 2 - MCV1 and MCV2.
#' @param using_sia A numeric indicator that determines whether supplementary
#' immunisation activities (SIAs) are implemented and how SIAs are distributed
#' between zero-dose and already-vaccinated populations: 0 - no SIA, 1 - SIAs
#' based on a weighted logistic function fitted with Portnoy's data, and 2 -
#' SIAs based on an assumption that 7.7% of the population are never reached by
#' vaccination.
#' @param folder_date Starting date of the simulation, as seen in the folder
#' name for RDS results, in the format "YYYYMMDD".
#' @param sim_years A numeric vector containing calendar years included for
#' model simulation.
#'
#' @import data.table
#'
#' @examples
#' \dontrun{
#' get_burden_estimate (
#'   scenario_name              = "campaign-only-default",
#'   save_scenario               = scenario_number,
#'   log_name                   = "test_log",
#'   vaccination                = 0,
#'   using_sia                  = 1,
#'   folder_date                = "20210930",
#'   sim_years                  = 1980:2100
#'   )
#'   }
get_burden_estimate <- function (
    scenario_name,
    save_scenario,
    log_name,
    vaccination,
    using_sia,
    folder_date,
    sim_years, 
    data
) {
  
  browser()
  
  # ----------------------------------------------------------------------------
  # merge and process results
  # ----------------------------------------------------------------------------
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
  years          <- as.numeric (sim_years)
  ages           <- c(0:100)
  template    	 <- copy (setDT (data_template)) [year %in% years]   # drop rows if simulation period is shorter
  report_years   <- sort (unique (template$year))
  country_names  <- unique (subset (template, select = c("country", "country_name")))
  c_names        <- country_names$country_name
  names(c_names) <- country_names$country
  
  # file name for burden estimates
  burden_estimate_file <- paste0 ("central_burden_estimate_", scenario_name)
  
  # coverage file
  coverage_routine <- copy(fread(paste0(o$pth$coverage,
                                        "routine_",
                                        scenario_name,
                                        ".csv"))
  )[year %in% sim_years]
  
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
  
  # --------------------------------------------------------------------------
  # add columns for remaining life expectancy & MCV1
  # --------------------------------------------------------------------------
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
  
  # --------------------------------------------------------------------------
  # calculate deaths and DALYs
  # --------------------------------------------------------------------------
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
  
  writelog (log_name, "burden estimate csv file generated")
  
  return ()
} # end of function -- merge_case_csv
# ------------------------------------------------------------------------------

