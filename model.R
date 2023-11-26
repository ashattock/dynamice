###########################################################
# MODEL
#
# Main transmission model functionality. 
#
###########################################################

# ------------------------------------------------------------------------------
# Run measles model for a given country and scenario
# ------------------------------------------------------------------------------
run_model = function(sim, data) {
  
  # Global model parameters
  p = prepare_params(data)
  
  # ---- Model set up ----
  
  # country-specific timeliness curve
  country_timeliness <- data$timeliness[!is.na(age), timeliness]
  timeliness_ages    <- data$timeliness[!is.na(age), age]
  
  # expand 0-2 years old to weekly age strata
  s         <- 52 # number of finer stages within an age band (weekly ages, so 52)
  jt        <- 3  # how many ages to expand to s (or to weeks)
  beta_full <- matrix (0, ncol = 254, nrow = 254)
  
  beta_full[(1:(s*jt)), (1:(s*jt))] <- expandMatrix (
    A = data$contact [1:jt, 1:jt]/s,  # needs to be divided by 52 so that the mean total number of contacts stays the same
    expand_rows =  s, expand_cols =  s,
    rescale_rows = FALSE, rescale_cols = FALSE)
  
  beta_full[1:(s*jt),((s*jt)+1):(ncol(beta_full))] <- expandMatrix(
    A = data$contact [1:jt,(jt+1):ncol(data$contact)],
    expand_rows = s, expand_cols = 1,
    rescale_rows = F, rescale_cols = F)
  
  beta_full[((s*jt)+1):(nrow(beta_full)), 1:(s*jt)] <- expandMatrix(
    A = data$contact [(jt+1):nrow(data$contact),1:jt]/s,  # adjust to ensure the mean total number of contacts stays the same
    expand_rows = 1, expand_cols = s,
    rescale_rows = F, rescale_cols = F)
  
  beta_full[((s*jt)+1):(nrow(beta_full)), ((s*jt)+1):(ncol(beta_full))] <-
    data$contact [(jt+1):nrow(data$contact),(jt+1):ncol(data$contact)]
  
  beta_full_unadj <- beta_full
  
  # infection rate under target R0
  beta_tstep  <- data$r0 * p$gamma
  
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
  
  # run model by yearly input
  for (y in o$years) {
    
    message("  - ", y)
    
    # ---- Spin up ----
    
    #old script groups those aged 70-80, but division is by actual population size
    pop.vector <- data$population[year == y, value]
    
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
    
    # Interpretation of set_routine
    #  0: no routine MCV
    #  1: MCV1 only
    #  2: MCV1 + MCV2
    
    if (sim$set_routine >= 1) {
      
      # Maximum coverage can (obviously) only be 100%
      # To estimate proportion that is vaccinated at each week, we first calculate the number of individuals remaining susceptible
      # Then we calculate the number of individuals that should be vaccinated each week, in order to remain 1 - coverage susceptibles at the end of the timeliness data
      # In essence, this becomes the inverse of the cumulative timeliness curve
      cycov <- data$coverage_routine[year == y & vaccine == "MCV1", coverage]
      
      # Not use the following adjustment for coverage, as it leads to reduced the ...
      # final size of vaccinated population under the assumption of perfect timeliness
      # cycov <- data$coverage_routine [year == y & vaccine == "MCV1", coverage] /
      #   data$timeliness [is.na(age), prop_final_cov]
      
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
    
    if (sim$set_routine == 2) {
      country_year_mcv2 <- data$coverage_routine [year == y & vaccine == "MCV2", coverage]
    } else {
      country_year_mcv2 <- 0
    }
    
    # if ( length (country_year_mcv2) == 0 || is.na(country_year_mcv2) ) {
    #   country_year_mcv2 <- 0
    # }
    
    # Interpretation of set_sia
    #  0: no SIA
    #  1: random reach (baseline assumption)
    #  2: 7.7% less likely to be reached at national level
    #  3: zero-dose first
    #  4: already-vaccinated first
    #  5: 7.7% less likely to be reached at subnational level
    
    # set up SIA inputs
    cy_coverage_sia <- data$coverage_sia[year == y]
    if ((sim$set_sia >= 1) && (dim(cy_coverage_sia)[1] > 0)) {
      # set up timesteps based on day of the year
      setorder(cy_coverage_sia, mid_day)
      sia_days <- cy_coverage_sia$mid_day
      t_sia_days <- c()
      for (iday in unique(sia_days)){
        t_sia_days <- c(t_sia_days,
                        round(p$tstep*(iday/365))-1 + (1:sum(sia_days == iday)))
      }
      
      sia_input <- list (
        sia_implement = sim$set_sia, # choose SIA implementation approach
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
  
  # ---- Save model output ----
  
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
  
  # Format into datatable
  output_dt = data.table(
    id       = sim$id, 
    year     = rep(o$years, each  = length(o$ages)),
    age      = rep(o$ages,  times = length(o$years)),
    cases0d  = as.vector(output$cases0d),
    cases1d  = as.vector(output$cases1d),
    cases2d  = as.vector(output$cases2d),
    pops     = as.vector(output$pops),
    pops0d   = as.vector(output$pops0d),
    popsSus  = as.vector(output$popsSus),
    doses    = as.vector(output$doses),
    reachs0d = as.vector(output$reachs0d),
    fvps     = as.vector(output$fvps))
  
  # Save model output
  saveRDS(output_dt, file = paste0(o$pth$sims, sim$id, ".rds"))
}

# ------------------------------------------------------------------------------
# Estimate disease burden from model outcomes
# ------------------------------------------------------------------------------
run_burden = function(sim, data) {
  
  # Load already simulated model output
  model_output = readRDS(paste0(o$pth$sims, sim$id, ".rds"))
  
  # Case fatality rate - some filling needed
  cfr_dt = expand_grid(year = o$years, 
                       age  = o$ages) %>%
    mutate(country = sim$country) %>%
    # Append CFR data...
    left_join(y  = data$cfr_portnoy_21, 
              by = c("year", "age", "country")) %>%
    # Fill missing years...
    group_by(age) %>%
    fill(cfr, .direction = "downup") %>%
    ungroup() %>%
    # Fill missing ages...
    group_by(year) %>%
    fill(cfr, .direction = "downup") %>%
    ungroup() %>%
    # Tidy up...
    select(year, age, cfr) %>%
    as.data.table()
  
  # Life years remaining - bound below by zero
  life_remain_dt = data$life_exp %>%
    filter(year %in% o$years) %>%
    expand_grid(age = o$ages) %>%
    mutate(life_remain = pmax(value - age, 0)) %>%
    select(year, age, life_remain) %>%
    as.data.table()
  
  # Calculate DALYs
  dalys_dt = model_output %>%
    mutate(cases = cases0d + cases1d + cases2d) %>%
    select(id, year, age, cases) %>%
    # Calculate number of deaths...
    left_join(y  = cfr_dt, 
              by = c("year", "age")) %>%
    mutate(deaths = cases * cfr) %>%
    # Calculate years of life lost (YLL)...
    left_join(y  = life_remain_dt, 
              by = c("year", "age")) %>%
    mutate(yll = deaths * life_remain) %>%
    # Calculate years of life with disability (YLD)...
    mutate(yld = o$disability_weight * (cases - deaths)) %>%
    # Finally, calculate DALYs...
    mutate(dalys = yll + yld) %>%
    # Tidy up...
    select(id, year, age, all_of(o$metrics)) %>%
    pivot_longer(cols = o$metrics,
                 names_to = "metric") %>%
    arrange(metric, year, age) %>%
    as.data.table()
  
  # Save disease burden estimates
  saveRDS(dalys_dt, file = paste0(o$pth$burden, sim$id, ".rds"))
}

