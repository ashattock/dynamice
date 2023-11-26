############################################################
# PLOTTING
#
# All plotting functionality in one place.
#
############################################################

# ---------------------------------------------------------
# Plot deaths and DALYs averted over time
# ---------------------------------------------------------
plot_burden_averted = function() {
  
  message(" > Plotting disease burden averted")
  
  # Load central results for all scenarios
  results_dt = load_central_results()
  
  # ---- Disease burden over time ----
  
  # Summarise over all ages
  burden_dt = results_dt %>%
    group_by(country_name, scenario, year, metric) %>%
    summarise(value = sum(value)) %>%
    ungroup() %>%
    arrange(metric, scenario, country_name, year) %>%
    as.data.table()
  
  # Iterate through key metrics
  for (metric in o$metrics) {
    
    # Subset for this metric
    metric_dt = burden_dt %>%
      filter(metric == !!metric)
    
    # Basic plot of disease burden over time
    g = ggplot(metric_dt) +
      aes(x = year, y = value, colour = scenario) +
      geom_line() +
      facet_wrap(~country_name, scales = "free_y")
    
    # Save figure to file
    save_fig(g, "Disease burden", metric)
  }
  
  # ---- Cumulative disease burden averted over time ----
  
  # Metrics averted relative to baseline
  averted_dt = burden_dt %>%
    # First cumulatively sum over time...
    group_by(country_name, scenario, metric) %>%
    mutate(value = cumsum(value)) %>%
    ungroup() %>%
    # Take the difference to baseline...
    group_by(country_name, year, metric) %>%
    mutate(averted = value[scenario == "nomcv"] - value) %>%
    ungroup() %>%
    # Remove baseline scenario...
    filter(scenario != "nomcv") %>%
    as.data.table()
  
  # Iterate through key metrics
  for (metric in o$metrics) {
    
    # Subset for this metric
    metric_dt = averted_dt %>%
      filter(metric == !!metric)
    
    # Disease burden averted over time
    g = ggplot(metric_dt) +
      aes(x = year, y = value, colour = scenario) +
      geom_line() +
      facet_wrap(~country_name, scales = "free_y")
    
    # Save figure to file
    save_fig(g, "Disease burden averted", metric)
  }
  
  # ---- Total disease burden averted ----
  
  # Total impact over all countries
  total_dt = averted_dt %>%
    group_by(scenario, year, metric) %>%
    summarise(value = sum(value)) %>%
    ungroup() %>%
    arrange(metric, scenario, year) %>%
    as.data.table()
  
  # Total disease burden over time for all metrics
  g = ggplot(total_dt) +
    aes(x = year, y = value, colour = scenario) +
    geom_line() + 
    facet_wrap(~metric, scales = "free_y")
  
  # Save figure to file
  save_fig(g, "Disease burden averted total")
}

# ---------------------------------------------------------
# xxx
# ---------------------------------------------------------
plot_coverage_all = function() {
  
  stop("Plotting functionality not yet integrated into this pipeline")

  # TODO: Use countrycode package here
  country_names <- c("India", "Nigeria", "Indonesia", "Ethiopia", "China",
                     "Philippines", "Uganda", "DRC", "Pakistan", "Angola",
                     "Madagascar", "Ukraine", "Malawi", "Somalia")
  names(country_names) <- o$countries

  # folder for saving figures
  res_folder <- paste0 (getwd(), "/previous_res/20230401/")


  # ------------------------------------------------------------------------------
  ## plot coverage trends for both routine immunisation and campaigns
  # ------------------------------------------------------------------------------
  outfile_mcv1_mcv2_sia  <- fread ("coverage/coverage_mcv1-mcv2-sia.csv")
  outfile_mcv1_mcv2alt1  <- fread ("coverage/coverage_mcv1-mcv2alt1.csv")
  outfile_mcv1_mcv2alt2  <- fread ("coverage/coverage_mcv1-mcv2alt2.csv")
  plt_data <- rbind (outfile_mcv1_mcv2_sia,
                     copy (outfile_mcv1_mcv2alt1 [vaccine == "MCV2"])[, vaccine := "MCV2 \n(early intro, fast rollout)"],
                     copy (outfile_mcv1_mcv2alt2 [vaccine == "MCV2"])[, vaccine := "MCV2 \n(early intro, gradual rollout)"])

  # update country names
  plt_data [country_code == "COD", country := "DRC"]
  plt_data [vaccine == "SIA", vaccine := "SIAs"]

  # rank countries by IHME burden
  plt_data [, country := factor (country, levels = country_names[o$countries])]

  # figure 1: MCV1, MCV2, SIAs coverage
  pdf (paste0 (res_folder, "figures/fig1_covall.pdf"), width = 14, height = 7)
  plt_cov <- ggplot (data = plt_data [vaccine %in% c("MCV1", "MCV2") & year >= 2000],
                     aes (x = year, y = coverage, colour = vaccine)) +
    scale_x_continuous (breaks = pretty_breaks ()) +
    geom_line (size = 0.9) +
    facet_wrap (vars(country), ncol = 5) +
    labs (title = " ", x = "Year", y = "Vaccine coverage") +
    theme_bw() +
    theme (legend.position  = c(0.9, 0.08),#"bottom"
           legend.direction = "vertical",
           legend.key.size = unit (1, 'cm'),
           legend.text = element_text (size = 13),
           legend.title = element_text (size = 14),
           axis.title.x = element_text (size = 15, vjust = -0.75),
           axis.title.y = element_text (size = 15, margin = margin(r = 15)),
           axis.text.x = element_text (size = 10, angle = 60, hjust = 1),
           axis.text.y = element_text (size = 10),
           strip.text.x = element_text (size = 12),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           plot.margin = margin (0, 0.2, 0.2, 0.2, "cm")) +
    geom_point (data = plt_data [vaccine == "SIAs" & year >= 2000],
                aes (x = year, y = coverage), size = 1.2) +
    scale_colour_manual ("Delivery strategy",
                         values = c("#42b540", "#00468b",  "#ed0000")) +
    guides(color = guide_legend (override.aes = list (shape = c(NA,NA,16),
                                                      linetype = c(1,1,NA),
                                                      size = c(1,1,1.2))))
  print(plt_cov)
  dev.off()


  # figure S1: MCV2, MCV2 (early intro) coverage
  pdf (paste0 (res_folder, "figures/figS1_mcv2alt-cov.pdf"), width = 12, height = 7)
  ggplot (data = plt_data [vaccine != "SIAs" & year >= 2000],
          aes (x = year, y = coverage, colour = vaccine, linetype = vaccine)) +
    scale_x_continuous (breaks = pretty_breaks ()) +
    geom_line (size = 1) +
    facet_wrap (vars(country), ncol = 5) +
    labs (title = " ", x = "Year", y = "Vaccine coverage") +
    scale_colour_manual (name = "Delivery strategy",
                         values = c("#42b540", "#00468b", "#0099b4", "#00468B99")) +
    scale_linetype_manual ("Delivery strategy", values = c(1,1,1,2)) +
    theme_bw() +
    theme (legend.position  ="bottom",
           legend.key.size = unit (2, 'cm'),
           legend.text = element_text (size = 13),
           legend.title = element_text (size = 14),
           axis.title.x = element_text (size = 15, vjust = -0.75),
           axis.title.y = element_text (size = 15, margin = margin(r = 15)),
           axis.text.x = element_text (size = 10, angle = 60, hjust = 1),
           axis.text.y = element_text (size = 10),
           strip.text.x = element_text (size = 12),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           plot.margin = margin (0, 0.2, 0.2, 0.2, "cm"))
  dev.off()
}

# ---------------------------------------------------------
# xxx
# ---------------------------------------------------------
plot_everything = function() {
  
  stop("Plotting functionality not yet integrated into this pipeline")

  vac_stgs <- c("nomcv",                 # (1) no vaccination
                "mcv1",                  # (2) MCV1 only
                "mcv1-mcv2",             # (3) MCV1 + MCV2
                "mcv1-mcv2-sia",         # (4) MCV1 + MCV2 + SIA
                "mcv1-sia",              # (5) MCV1 + SIA
                "mcv1-mcv2alt1",         # (6) MCV1 + MCV2(early intro, fast rollout)
                "mcv1-mcv2alt1-sia",     # (7) MCV1 + MCV2(early introm fast rollout) + SIA
                "mcv1-mcv2-siaalt1",     # (8) MCV1 + MCV2 + SIA(zero dose first)
                "mcv1-mcv2-siaalt2",     # (9) MCV1 + MCV2 + SIA(already vaccinated first)
                "mcv1-siaalt1",          # (10) MCV1 + SIA(zero dose first)
                "mcv1-siaalt2",          # (11) MCV1 + SIA(already vaccinated first)
                "mcv1-mcv2alt1-siaalt1", # (12) MCV1 + MCV2(early intro) + SIA(zero dose first)
                "mcv1-mcv2alt2"          # (13) MCV1 + MCV2(early intro, gradual rollout)
  )
  vac_stg_names <- c("No vaccination",
                     "MCV1",
                     "MCV1 + MCV2",
                     "MCV1 + MCV2 + SIAs",
                     "MCV1 + SIAs",
                     "MCV1 + MCV2 (early intro, fast rollout)",
                     "MCV1 + MCV2 (early intro, fast rollout) + SIAs",
                     "MCV1 + MCV2 + SIAs (zero-dose first)",
                     "MCV1 + MCV2 + SIAs (vaccinated first)",
                     "MCV1 + SIAs (zero-dose first)",
                     "MCV1 + SIAs (vaccinated first)",
                     "MCV1 + MCV2 (early intro) + SIAs (zero-dose first)",
                     "MCV1 + MCV2 (early intro, gradual rollout)")

  eva_ctries <- c("IND", "NGA", "IDN", "ETH", "CHN",
                  "PHL", "UGA", "COD", "PAK", "AGO",
                  "MDG", "UKR", "MWI", "SOM")

  country_names <- c("India", "Nigeria", "Indonesia", "Ethiopia", "China",
                     "Philippines", "Uganda", "DRC", "Pakistan", "Angola",
                     "Madagascar", "Ukraine", "Malawi", "Somalia")

  names (eva_ctries) <- country_names
  names (country_names) <- eva_ctries

  custom_palette <- c("#00468BFF", "#ED0000FF", "#42B540FF", "#0099B4FF", "#925E9FFF",
                      "#FDAF91FF", "#AD002AFF", "#ADB6B6FF", "#1B1919FF", "#00468B99",
                      "#ED000099", "#42B54099", "#0099B499", "#925E9F99", "#FDAF9199",
                      "#AD002A99", "#ADB6B699", "#1B191999", "#00468B66", "#ED000066")

  res_folder <- paste0 (getwd(), "/previous_res/20230401/")
  output_folder <- paste0 (getwd(), "/tables_figures/")

  # ------------------------------------------------------------------------------
  ## load and combine model outputs
  # ------------------------------------------------------------------------------
  # burden and vaccine doses
  file_burden  <- NULL
  for (scname in vac_stgs){
    scn_burden <- fread (paste0 (res_folder, "siareach_2/central_burden_estimate_", scname, ".csv"))
    scn_burden [, `:=` (cases = cases0d + cases1d + cases2d,
                        deaths = deaths0d + deaths1d + deaths2d)]
    scn_burden [, comp := scname]
    file_burden <- rbind (file_burden, scn_burden)
  }
  file_burden [country == "COD", country_name := "DRC"]

  cum_burden <- file_burden [, lapply (.SD, sum),
                             .SDcols = pops:deaths,
                             by = c("country", "year", "country_name", "comp")]

  # calculate incidence rate (case per million)
  cum_burden [, incrateM := (cases/pops)*1e6]


  # ------------------------------------------------------------------------------
  ## set up a function for averted burden and number needed to vaccinate (NNV)
  # ------------------------------------------------------------------------------
  cal_avtnnv <- function (cum_burden, comp_base, comp_intv){
    total_burden <- cum_burden [, .(total_cases = sum(cases),
                                    total_deaths = sum(deaths),
                                    total_dalys = sum(dalys),
                                    total_doses = sum(doses)),
                                by = c("comp", "country_name", "country")]
    total_burden <- rbind (total_burden,
                           total_burden [,  .(country_name = "Global",
                                              country = "Global",
                                              total_cases = sum(total_cases),
                                              total_deaths = sum(total_deaths),
                                              total_dalys = sum(total_dalys),
                                              total_doses = sum(total_doses)),
                                         by = c("comp")])

    merge_dat <- total_burden [comp == comp_base][total_burden [comp == comp_intv],
                                                  on = .(country_name, country)]
    merge_dat <- merge_dat [, .(comp_set = paste0 (i.comp, "_VS_", comp),
                                country_name, country,
                                avt_cases  = total_cases - i.total_cases,
                                avt_deaths = total_deaths - i.total_deaths,
                                avt_dalys  = total_dalys - i.total_dalys,
                                add_doses  = i.total_doses - total_doses,
                                pr_red_cases = (total_cases - i.total_cases)/total_cases)]
    merge_dat [, nnv := add_doses/avt_cases]
    return(merge_dat)
  }


  # ------------------------------------------------------------------------------
  ## load WHO and IHME data for comparison
  # ------------------------------------------------------------------------------
  # WHO reported cases
  # https://immunizationdata.who.int/pages/incidence/MEASLES.html
  input_WHOcase <- read_excel ("D:/research-data/Measles reported cases and incidence by year (Reported cases) 2023-031-03 0-18 UTC.xlsx")
  input_WHOcase$`Country / Region` [which (input_WHOcase$`Country / Region` == "Democratic Republic of the Congo")] <- "DRC"
  input_WHOcase <- setDT (input_WHOcase) [`Country / Region` %in% country_names, !(`1999`:`1980`)]

  dat_WHOcase <- tidyr::pivot_longer (input_WHOcase [, !c("Disease", "2022", "2021")], cols = `2020`:`2000`,
                                      names_to = "year", values_to = "notifs")
  dat_WHOcase <- setDT (dat_WHOcase) [, year := as.numeric(year)]
  setnames (x = dat_WHOcase, old = c("Country / Region"), new = c("country_name"))

  # adjust data format
  dat_WHOcase [, `:=` (notifs = as.numeric (str_remove_all (notifs, "\\,")),
                       country_name = factor (country_name, levels = country_names),
                       country = factor (eva_ctries [country_name], levels = eva_ctries))]
  dat_WHOcase <- dat_WHOcase [cum_burden [comp == "nomcv" & year <= 2020,
                                          .(comp, country, year, pops)],
                              on = .(country = country, year = year)]
  dat_WHOcase [, incrateM := 1e6*(notifs/pops)]
  setorder (dat_WHOcase, country_name, year)

  # IHME GBD-2019
  # http://ghdx.healthdata.org/gbd-results-tool
  input_IHME <- fread( "D:/research-data/IHME-GBD_2019_DATA-530a6126-1.csv")[measure_name == "Incidence"]
  input_IHME [location_name == "Democratic Republic of the Congo", location_name := "DRC"]

  # # calculate the proportion of global measles burden
  # sum(input_IHME [year %in% 2010:2019 & `location_name` %in% country_names, val])/
  #   sum(input_IHME [year %in% 2010:2019, val])
  # # 78.0%

  dat_IHME <- input_IHME [year %in% 2000:2020, .SD,
                          .SDcols = c("location_name", "year", "val", "upper", "lower")]
  setnames (x = dat_IHME, old = c("location_name", "val"),
            new = c("country_name", "est_cases"))
  dat_IHME [, country_name := factor(country_name, levels = country_names)]

  # adjust data format
  dat_IHME [, `:=` (country_name = factor (country_name, levels = country_names),
                    country = factor (eva_ctries [country_name], levels = eva_ctries))]
  dat_IHME <- dat_IHME [cum_burden [comp == "nomcv" & year < 2020,
                                    .(comp, country, year, pops)],
                        on = .(country = country, year = year)]
  dat_IHME[, incrateM := 1e6*(est_cases/pops)]
  setorder (dat_IHME, country_name, year)


  # ------------------------------------------------------------------------------
  ## set up plotting style functions and data
  # ------------------------------------------------------------------------------
  # add a blank window
  plt_blank <- ggplot() + geom_blank(aes(1,1)) +
    theme (plot.background = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           axis.title.x = element_blank(),
           axis.title.y = element_blank(),
           axis.text.x = element_blank(),
           axis.text.y = element_blank(),
           axis.ticks = element_blank(),
           axis.line = element_blank())

  # show non-scientific numbers for axis and allow decimals for small values only
  format_num_plain <- function (x, ...) {
    format (x, ..., scientific = FALSE, drop0trailing = TRUE)
  }

  # categorise by MCV2 state
  pltdat_mcv2_text <- data.table (country = eva_ctries, country_name = country_names)
  pltdat_mcv2_text  [ , `:=` (mcv2_intro = ifelse (country %in% c("UGA","COD", "SOM"),
                                                   "No MCV2 intro", character(0)),
                              country_name = factor (country_name,
                                                     levels = country_names),
                              comp = factor (vac_stg_names[3],
                                             levels = vac_stg_names [c(1,2,3,5,4)]))]

  # country order by MCV2 introduction
  eva_ctries_mcv2 <- c("IND", "IDN", "CHN", "PHL", "PAK", "AGO", "UKR", "MWI",
                       "NGA","ETH", "MDG",
                       "UGA","COD", "SOM")


  # ------------------------------------------------------------------------------
  ## plot country burden over time
  # ------------------------------------------------------------------------------
  pltdata_ctry_burden <- cum_burden [comp %in% vac_stgs[c(1,2,3,5,4)],
                                     .(comp, year, country_name, cases, deaths)]
  pltdata_ctry_burden <- setDT (pivot_longer (pltdata_ctry_burden,
                                              cols = cases:deaths,
                                              names_to = "measure",
                                              values_to = "value"))
  pltdata_ctry_burden [, `:=` (measure = factor (measure, levels = c("cases", "deaths"),
                                                 labels = c("Cases", "Deaths")),
                               country_name = factor (country_name, levels = country_names[eva_ctries]),
                               comp = factor (comp, levels = vac_stgs[c(1,2,3,5,4)],
                                              labels = vac_stg_names[c(1,2,3,5,4)]))]

  pdf (paste0 (output_folder, "figS3_burden-trend.pdf"), height = 6.5, width = 14)
  ggplot (data = pltdata_ctry_burden,
          aes(x = year, y = value/1e6, fill = country_name)) +
    geom_area () +
    facet_grid (rows = vars(measure), cols = vars(comp), scales = "free_y") +
    labs (x = "Year", y = "Estimated health burden (millions)", fill = "Country") +
    scale_fill_manual (values = custom_palette) +
    #scale_x_continuous (labels = label_number (accuracy = 1)) +
    theme_bw () +
    theme (legend.text = element_text (size = 14),
           legend.title = element_text (size = 15),
           axis.title.x = element_text (size = 14, margin = margin (t = 10)),
           axis.title.y = element_text (size = 14, margin = margin (r = 10)),
           axis.text.x = element_text (size = 10, angle = 60, hjust = 1),
           axis.text.y = element_text (size = 10),
           strip.text.x = element_text (size = 14),
           strip.text.y = element_text (size = 14),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           plot.margin = margin (0.2, 0.2, 0.2, 0.2, "cm"))
  dev.off()


  # ------------------------------------------------------------------------------
  ## plot incidence trends by scenarios
  # ------------------------------------------------------------------------------
  # adjust scenario order to allow 'MCV1 only' at the top
  pltdat_incrateM <- cum_burden [comp %in% vac_stgs[c(1,2,3,5,4)]]
  pltdat_incrateM [, `:=` (country_name = factor (country_name,
                                                  levels = country_names[eva_ctries]),
                           comp = factor (comp, levels = vac_stgs [c(1,2,3,5,4)],
                                          labels = vac_stg_names [c(1,2,3,5,4)]))]
  # manually specify text location for the note on MCV2 introduction
  pltdat_mcv2_text [, `:=` (inc_x = c(NA, NA, NA, NA, NA,
                                      NA, 2014, 2012.5, NA, NA,
                                      NA, NA, NA, 2015),
                            inc_y = c(NA, NA, NA, NA, NA,
                                      NA, 2e4, 2.35e4, NA, NA,
                                      NA, NA, NA, 2.85e4))]

  pdf (paste0 (output_folder, "fig2_incrateM.pdf"), width = 14, height = 7)
  ggplot (data = pltdat_incrateM,
          aes(x = year, y = incrateM, colour = comp)) +
    geom_line (size = 1, alpha = 1) +
    facet_wrap (vars(country_name), scales = "free", ncol = 5) +
    scale_colour_manual (values = c("grey 50", custom_palette[c(3,1,2,5)]),
                         guide = guide_legend (reverse = FALSE)) +
    geom_text (data = pltdat_mcv2_text, size = 4,
               aes (x = inc_x, y = inc_y, label = mcv2_intro, colour = comp),
               show.legend = FALSE) +
    labs (x = "Year", y = "Estimated measles incidence rate \nper million population",
          colour = "Delivery strategies") +
    theme_bw () +
    theme (legend.position = c(0.9, 0.15),
           legend.text = element_text (size = 12),
           legend.title = element_text (size = 12),
           legend.key.width = unit (2.2, "line"),
           axis.title.x = element_text (size = 15, margin = margin(t = 15)),
           axis.title.y = element_text (size = 15, margin = margin(r = 15)),
           axis.text.x = element_text (size = 10.5),
           axis.text.y = element_text (size = 10.5),
           strip.text.x = element_text (size = 12))
  dev.off()


  # Compare the 'optimal-use' scenarios
  pltdat_incrateM_opt <- cum_burden [comp %in% vac_stgs[c(4,7,8,12)]]
  pltdat_incrateM_opt [, `:=` (country_name = factor (country_name,
                                                      levels = country_names[eva_ctries]),
                               comp = factor (comp, levels = vac_stgs [c(4,7,8,12)],
                                              labels = vac_stg_names [c(4,7,8,12)]))]

  pdf (paste0 (output_folder, "figS5_incrateM_optim.pdf"), width = 14, height = 8)
  ggplot (data = pltdat_incrateM_opt,
          aes(x = year, y = incrateM, colour = comp)) +
    geom_line (size = 0.9, alpha = 1) +
    facet_wrap (vars(country_name), scales = "free", ncol = 5) +
    scale_colour_manual (name = "MCV2 introduction and SIA delivery strategies",
                         values = custom_palette[c(5,4,6,8)],
                         guide = guide_legend (order = 1, ncol = 2, byrow = T)) +
    scale_y_log10 (labels = format_num_plain) +
    labs (x = "Year", y = "Estimated measles incidence rate \nper million population (log scale)") +
    theme_bw () +
    theme (legend.position = "bottom",
           legend.text = element_text (size = 11.5),
           legend.title = element_text (size = 12),
           plot.margin = unit (c(0.1, 0.25, 0.1, 0.1), "cm"),
           axis.title.x = element_text (size = 15, margin = margin(t = 15)),
           axis.title.y = element_text (size = 15, margin = margin(r = 10)),
           axis.text.x = element_text (size = 11),
           axis.text.y = element_text (size = 11),
           strip.text.x = element_text (size = 12.5))
  dev.off()


  # ------------------------------------------------------------------------------
  ## plot susceptible population and birth cohort by scenarios
  # ------------------------------------------------------------------------------
  # birth cohort: not changed by scenario
  pltdat_sus_age0 <- copy (file_burden [age == 0 & comp == vac_stgs[1],
                                        c("country_name", "year", "comp", "pops")]) [, comp := "birth"]
  pltdat_sus_age0  [, `:=` (country_name = factor (country_name,
                                                   levels = country_names[eva_ctries]),
                            data_type = "Annual birth cohort")]
  pltdat_sus <- file_burden [age < 5 & comp %in% vac_stgs[c(1,2,3,5,4)],
                             lapply (.SD, sum), .SDcols = "popsSus",
                             by = c("country", "year", "country_name", "comp")]
  pltdat_sus [, `:=` (country_name = factor (country_name,
                                             levels = country_names[eva_ctries]),
                      comp = factor (comp, levels = c(vac_stgs[c(1,2,3,5,4)]),
                                     labels = c(vac_stg_names[c(1,2,3,5,4)])))]
  # manually specify text position for the note on MCV2 introduction
  pltdat_mcv2_text [, sus_y := c(NA, NA, NA, NA, NA,
                                 NA, 2.9, 6.25, NA, NA,
                                 NA, NA, NA, 1.05)]

  # plot
  pdf (paste0 (output_folder, "fig3_sus-birth.pdf"), width = 14, height = 7)
  ggplot () +
    geom_line (data = pltdat_sus, size = 1,
               aes(x = year, y = popsSus/1e6, colour = comp)) +
    facet_wrap (vars(country_name), scales = "free", ncol = 5) +
    scale_colour_manual (values = c("grey 50", custom_palette[c(3,1,2,5)]),
                         guide = guide_legend (reverse = FALSE, order = 1)) +
    geom_text (data = pltdat_mcv2_text, size = 4,
               aes (x = 2005, y = sus_y, label = mcv2_intro, colour = comp),
               show.legend = FALSE) +
    labs (x = "Year", y = "Estimated number of susceptible population < 5 y/o (millions)",
          colour = "Delivery strategies") +
    theme_bw () +
    theme (legend.position = c(0.91, 0.12),
           legend.text = element_text (size = 12),
           legend.title = element_text (size = 12),
           legend.key.width = unit (2.2, "line"),
           axis.title.x = element_text (size = 15, margin = margin(t = 15)),
           axis.title.y = element_text (size = 15, margin = margin(r = 15)),
           axis.text.x = element_text (size = 10.5),
           axis.text.y = element_text (size = 10.5),
           strip.text.x = element_text (size = 12)) +
    ggnewscale::new_scale("linetype") +
    geom_line (data = pltdat_sus_age0,
               aes(x = year, y = pops/1e6, linetype = data_type),
               size = 0.9, alpha = 0.85) +
    scale_linetype_manual (name = NULL, values = 3,
                           guide = guide_legend (order = 2))
  dev.off()


  # calculate number of years with susceptibles > birth cohort
  compdat_sus <- pltdat_sus [pltdat_sus_age0 [, !c("comp")], on = c("country_name", "year")]
  compdat_sus [, outbreak := ifelse(popsSus >= pops, 1, 0)]
  tab_sus <- compdat_sus [, .(pr_outbreak = 1-sum(outbreak)/21), by = c("country", "comp")]
  setorder (tab_sus, -comp)
  tab_sus <- setDT (pivot_wider (tab_sus, values_from = pr_outbreak, names_from = comp))
  tab_sus [, country := factor (country, levels = eva_ctries)]
  setorder (tab_sus, country)

  # get median and 25th-75th percentile
  tab_sus <- rbind (tab_sus,
                    cbind (data.table (country = c("pr25", "median", "pr75")),
                           tab_sus [, lapply (.SD, quantile, prob = c(.25, .5, .75),
                                              na.rm = TRUE), .SDcols = !"country"]))
  fwrite (x = tab_sus,
          file = paste0 (output_folder, "tabs3_yr-sus-outbreak.csv"))


  # ------------------------------------------------------------------------------
  ## calculate averted burden and number needed to vaccinate (NNV)
  # ------------------------------------------------------------------------------
  # check absolute case reduction
  sel_vacc_impact <- cal_avtnnv (cum_burden, "nomcv", "mcv1")
  setorder (sel_vacc_impact, pr_red_cases)
  sel_vacc_impact [country != "Global"]
  setorder (sel_vacc_impact, avt_cases)
  sel_vacc_impact [country != "Global"]

  cal_avtnnv (cum_burden, "mcv1-mcv2", "mcv1-mcv2alt1") [country == "Global", avt_cases] # 96804681
  cal_avtnnv (cum_burden, "mcv1-mcv2", "mcv1-mcv2alt2") [country == "Global", avt_cases] # 74521266


  # combine results of all comparison pairs
  all_avtnnv <- rbind (cal_avtnnv (cum_burden, "nomcv", "mcv1"),
                       cal_avtnnv (cum_burden, "nomcv", "mcv1-sia"),
                       cal_avtnnv (cum_burden, "nomcv", "mcv1-mcv2"),
                       cal_avtnnv (cum_burden, "nomcv", "mcv1-mcv2-sia"),
                       cal_avtnnv (cum_burden, "mcv1", "mcv1-sia"),
                       cal_avtnnv (cum_burden, "mcv1", "mcv1-mcv2"),
                       cal_avtnnv (cum_burden, "mcv1", "mcv1-mcv2-sia"),
                       cal_avtnnv (cum_burden, "mcv1-sia", "mcv1-mcv2-sia"),
                       cal_avtnnv (cum_burden, "mcv1-mcv2", "mcv1-mcv2-sia"),
                       cal_avtnnv (cum_burden, "mcv1", "mcv1-mcv2alt1"),
                       cal_avtnnv (cum_burden, "mcv1", "mcv1-mcv2alt2"),
                       cal_avtnnv (cum_burden, "mcv1", "mcv1-siaalt1"),
                       cal_avtnnv (cum_burden, "mcv1", "mcv1-siaalt2"))

  # averted cases, deaths, and dalys
  output_case <- setDT (pivot_wider (all_avtnnv [, .(comp_set, country_name, country, avt_cases)],
                                     values_from = avt_cases,
                                     names_from = comp_set))
  output_case [, country := factor (country, levels = c(eva_ctries, "Global"))]
  setorder (output_case, country)

  output_death <- setDT (pivot_wider (all_avtnnv [, .(comp_set, country_name, country, avt_deaths)],
                                      values_from = avt_deaths,
                                      names_from = comp_set))
  output_death [, country := factor (country, levels = c(eva_ctries, "Global"))]
  setorder (output_death, country)

  output_daly <- all_avtnnv [, .(comp_set, country_name, country, avt_dalys)]
  output_daly [, avt_dalys_K := avt_dalys/1000]
  output_daly <- setDT (pivot_wider (output_daly [, !c("avt_dalys")],
                                     values_from = avt_dalys_K,
                                     names_from = comp_set))
  output_daly [, country := factor (country, levels = c(eva_ctries, "Global"))]
  setorder (output_daly, country)


  # Table 1: Averted cases
  value_cols <- names(output_case)[-c(1:2)]
  sel_cols_tab1 <- c("country_name", "mcv1-mcv2-sia_VS_nomcv", "mcv1_VS_nomcv",
                     "mcv1-sia_VS_mcv1", "mcv1-mcv2_VS_mcv1", "mcv1-mcv2-sia_VS_mcv1")

  output_case [, (value_cols) := lapply (.SD,
                                         function (avt_burden){
                                           return (avt_burden/1000)}),
               .SDcols = value_cols]
  output_case [, country := factor (country, levels = c(eva_ctries_mcv2, "Global"))]
  setorder (output_case, country)
  fwrite (x = output_case [, ..sel_cols_tab1],
          file = paste0 (output_folder, "tab1_avtcase.csv"))


  # Table S4: Averted deaths
  output_death [, (value_cols) := lapply (.SD, function (avt_burden){
    return (avt_burden/1000)}),
    .SDcols = value_cols]
  output_death [, country := factor (country, levels = c(eva_ctries_mcv2, "Global"))]
  setorder (output_death, country)
  fwrite (x = output_death[, ..sel_cols_tab1],
          file = paste0 (output_folder, "tabs4_avtdeath.csv"))

  # Table 2: NNV
  # not calculated for countries have not yet introduced MCV2 when MCV1 only is the comparator
  output_nnv <- setDT (pivot_wider (all_avtnnv [comp_set %in% c("mcv1_VS_nomcv",
                                                                "mcv1-sia_VS_mcv1",
                                                                "mcv1-mcv2_VS_mcv1",
                                                                "mcv1-mcv2-sia_VS_mcv1-mcv2",
                                                                "mcv1-mcv2-sia_VS_mcv1-sia"),
                                                .(comp_set, country_name, country, nnv)],
                                    values_from = nnv,
                                    names_from = comp_set))
  output_nnv [, country := factor (country, levels = c(eva_ctries_mcv2, "Global"))]
  setorder (output_nnv, country)

  output_nnv <- rbind (output_nnv,
                       output_nnv [country_name != "Global",
                                   lapply (.SD, function(x) median (x, na.rm = T)),
                                   .SDcols = `mcv1_VS_nomcv`:`mcv1-mcv2-sia_VS_mcv1-mcv2`],
                       fill = TRUE)
  output_nnv [is.na(country), `:=` (country_name = "Median", country = "Median")]
  fwrite (x = output_nnv [country != "Global"],
          file = paste0 (output_folder, "tab2_nnv.csv"))

  # report IQRs for NNV estimates
  quantile (all_avtnnv [country_name != "Global" & comp_set == "mcv1-mcv2_VS_mcv1", nnv],
            prob = c(.25, .5, .75), na.rm = TRUE)
  quantile (all_avtnnv [country_name != "Global" & comp_set == "mcv1-sia_VS_mcv1", nnv],
            prob = c(.25, .5, .75))
  quantile (all_avtnnv [country_name != "Global" & comp_set == "mcv1-mcv2-sia_VS_mcv1-mcv2", nnv],
            prob = c(.25, .5, .75))
  quantile (all_avtnnv [country_name != "Global" & comp_set == "mcv1-mcv2-sia_VS_mcv1-sia", nnv],
            prob = c(.25, .5, .75), na.rm = TRUE)
  quantile (all_avtnnv [country_name != "Global" & comp_set == "mcv1-mcv2alt1_VS_mcv1", nnv],
            prob = c(.25, .5, .75), na.rm = TRUE)
  quantile (all_avtnnv [country_name != "Global" & comp_set == "mcv1-siaalt1_VS_mcv1", nnv],
            prob = c(.25, .5, .75))
  quantile (all_avtnnv [country_name != "Global" & comp_set == "mcv1-siaalt2_VS_mcv1", nnv],
            prob = c(.25, .5, .75))


  # ------------------------------------------------------------------------------
  ## plot dose distribution
  # ------------------------------------------------------------------------------
  # process data for doses
  pltdat_dose <- copy (cum_burden [, .(comp, country_name, country, year,
                                       doses, reachs0d, fvps)])
  pltdat_dose [, `:=` (doses0 = reachs0d,
                       doses1 = fvps,
                       doses2 = doses-reachs0d-fvps)]
  pltdat_dose <- setDT (pivot_longer (pltdat_dose [, !c("doses", "fvps", "reachs0d")],
                                      col = doses0:doses2,
                                      names_to = "measure", values_to = "value"))
  pltdat_dose [, `:=` (country_name = factor (country_name,
                                              levels = country_names[eva_ctries]),
                       measure = factor (measure,
                                         levels = c("doses0", "doses1", "doses2"),
                                         labels = c("zero-dose", "single-dose", "multi-dose")))]

  # plot cumulative doses over 2000-2020 by scenarios
  pltdat_dose_sum <- pltdat_dose [, .(total_dose = sum(value)),
                                  by = c("comp", "country_name", "country", "measure")]
  plt_dose_sum <- function (sel_ctries, sel_scns){
    pltdat_tmp <- pltdat_dose_sum [country %in% sel_ctries & comp %in% vac_stgs[sel_scns]
                                   & total_dose > 0] # remove bar borders in the plot for zero values
    pltdat_tmp [, comp := factor (comp, levels = vac_stgs[sel_scns],
                                  labels = vac_stg_names[sel_scns])]
    plt <- ggplot (data = pltdat_tmp, aes(x = comp, y = total_dose/1e6)) +
      geom_col (aes(fill = measure), width = 0.8, colour = "white", size = 0,
                position = position_stack (reverse = TRUE)) +
      facet_wrap (vars(country_name), scale = "free_y", ncol = 5) +
      labs (x = "", y = "Estimated number of doses over 2000-2020 (millions)") +
      scale_fill_manual ("Predicted vaccination state of population reached",
                         values = c("#42b540", "#00468b", "#ed0000")) +
      theme_bw () +
      theme (legend.position  = "top",
             legend.text = element_text (size = 14),
             legend.title = element_text (size = 15),
             axis.title.y = element_text (size = 15, vjust = 2),
             axis.text.x = element_text (size = 10, angle = 60, hjust = 1),
             axis.text.y = element_text (size = 10),
             strip.text.x = element_text (size = 12),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.background = element_blank(),
             plot.margin = margin (0.1, 0.2, 0, 0.6, "cm"))
    return(plt)
  }
  pdf (paste0 (output_folder, "figS4_dose-sum.pdf"), width = 11, height = 7)
  print (plt_dose_sum (eva_ctries, c(2,3,5,4)))
  dev.off()


  # ------------------------------------------------------------------------------
  ## plot averted cases and NNV for sensitivity analysis
  # ------------------------------------------------------------------------------
  all_avtnnv_senanl <- all_avtnnv [comp_set %in% c("mcv1-mcv2_VS_mcv1",
                                                   "mcv1-mcv2alt1_VS_mcv1",
                                                   "mcv1-mcv2alt2_VS_mcv1",
                                                   "mcv1-sia_VS_mcv1",
                                                   "mcv1-siaalt1_VS_mcv1",
                                                   "mcv1-siaalt2_VS_mcv1")]
  all_avtnnv_senanl [, `:=` (country_name = factor (country_name,
                                                    levels = country_names[eva_ctries]),
                             comp_set = factor (comp_set,
                                                levels = c("mcv1-mcv2_VS_mcv1",
                                                           "mcv1-mcv2alt1_VS_mcv1",
                                                           "mcv1-mcv2alt2_VS_mcv1",
                                                           "mcv1-sia_VS_mcv1",
                                                           "mcv1-siaalt1_VS_mcv1",
                                                           "mcv1-siaalt2_VS_mcv1"),
                                                labels = c("MCV1 + MCV2",
                                                           "MCV1 + MCV2 (early intro, fast rollout)",
                                                           "MCV1 + MCV2 (early intro, gradual rollout)",
                                                           "MCV1 + SIAs",
                                                           "MCV1 + SIAs (zero-dose first)",
                                                           "MCV1 + SIAs (vaccinated first)")),
                             avt_cases_M = ifelse (avt_cases <= 0.0001, NA, avt_cases/1e6))]

  # plot
  plt_senanl_bar <- function (plt_data, sel_mea, sel_ylab,
                              sel_lgdpos, sel_lgdtitle, sel_fill_values){
    plt <- ggplot (data = plt_data [country != "Global"],
                   aes (x = country_name, y = get(sel_mea))) +
      geom_col (aes (fill = comp_set, colour = comp_set), width = 0.7,  size = 0,
                position = position_dodge()) +
      labs (x = "", y = sel_ylab, fill = sel_lgdtitle, colour = sel_lgdtitle) +
      scale_fill_manual (values = sel_fill_values,
                         guide = guide_legend (ncol = 3, byrow = T)) +
      scale_colour_manual (values = sel_fill_values,
                           guide = guide_legend (ncol = 3, byrow = T)) +
      theme_bw () +
      theme (legend.position = sel_lgdpos,
             legend.text = element_text (size = 13),
             legend.title = element_text (size = 14),
             axis.title.y = element_text (size = 14, vjust = 2),
             axis.text.x = element_text (size = 12, angle = 60, hjust = 1),
             axis.text.y = element_text (size = 10),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.background = element_blank(),
             plot.margin = margin (0, 0.2, 0, 0.6, "cm"))
    return(plt)
  }

  ggsave (paste0 (output_folder, "fig4_nnv-avtcase-senanl.pdf"),
          ggarrange (plt_blank,
                     plt_senanl_bar (all_avtnnv_senanl, "avt_cases_M",
                                     "Estimated number of \naverted cases (millions)",
                                     c(0.5,1.33), "Vaccination strategy compared to MCV1 only",
                                     custom_palette [c(1,4,10,2,6,7)]),
                     plt_senanl_bar (all_avtnnv_senanl, "nnv",
                                     "Estimated number\nneeded to vaccinate",
                                     "none", NA,
                                     custom_palette [c(1,4,10,2,6,7)]),
                     ncol = 1, vjust = -1.2,
                     labels = c("", "A", "B"), heights = c(1.5, 4, 4)),
          height = 8, width = 11)


  # ------------------------------------------------------------------------------
  ## plot case and NNV estimates by different SIA delivery methods
  # ------------------------------------------------------------------------------
  # load data
  all_siareach_senanl <- NULL
  for (isia in c(1,2,5)) {
    for (iscn in vac_stgs[c(2,3,4,5)]){
      outputs <- fread (paste0 (res_folder, "siareach_", isia,
                                "/central_burden_estimate_", iscn, ".csv"))
      all_siareach_senanl <- rbind (all_siareach_senanl,
                                    outputs [, `:=` (scenario = iscn,
                                                     siareach = isia)])
      remove(outputs)
    }
  }
  siareach_methods <- c("7.7% less-likely-to-be-reached (national level)",
                        "7.7% less-likely-to-be-reached (subnational level)",
                        "random reach")

  all_siareach_senanl [, `:=` (country_name = factor (country_names[country], levels = country_names),
                               cases = cases0d + cases1d + cases2d,
                               deaths = deaths0d + deaths1d + deaths2d,
                               siareach = factor (siareach, levels = c(2,5,1), labels = siareach_methods))]
  cum_siareach_senanl <- all_siareach_senanl [, lapply (.SD, sum),
                                              .SDcols = c("pops", "cases", "deaths", "dalys", "doses"),
                                              by = c("year", "country", "country_name", "scenario", "siareach")]
  setorder (cum_siareach_senanl, country, year, siareach)

  # calculate NNV and add results of the main scenarios
  nnv_siareach_senanl <- NULL
  for (isia in siareach_methods){
    output_nnv <- cum_siareach_senanl [siareach == isia & scenario %in% vac_stgs[c(2,3,5)]]
    output_nnv [, comp := scenario]
    nnv_siareach_senanl <- rbind (nnv_siareach_senanl,
                                  cal_avtnnv (output_nnv, "mcv1", "mcv1-mcv2")[, siareach := isia],
                                  cal_avtnnv (output_nnv, "mcv1", "mcv1-sia")[, siareach := isia])
    remove (output_nnv)
  }
  nnv_siareach_senanl [, `:=` (comp_set = factor (comp_set, labels = c("MCV1 + MCV2", "MCV1 + SIAs")),
                               siareach = factor (siareach, levels = siareach_methods))]

  # plot case estimates
  pdf (paste0 (output_folder, "figS6_case-siareach.pdf"), width = 14, height = 8.5)
  ggplot(data = cum_siareach_senanl [scenario == "mcv1-mcv2-sia"],
         aes(x = year, y = cases/1e3, colour = siareach)) +
    scale_x_continuous (breaks = scales::pretty_breaks ()) +
    scale_colour_manual (name = "SIA dose delivery methods",
                         values = custom_palette[5:7]) +
    geom_line (size = 0.9) +
    facet_wrap (vars(country_name), nrow = 3, scales = "free") +
    labs (x = "Year", y = "Cases (thousand)") +
    theme_bw () +
    theme (legend.position = "bottom",
           legend.direction = "vertical",
           legend.text = element_text (size = 11.5),
           legend.title = element_text (size = 12),
           plot.margin = unit (c(0.1, 0.25, 0.1, 0.1), "cm"),
           axis.title.x = element_text (size = 15, margin = margin(t = 15)),
           axis.title.y = element_text (size = 15, margin = margin(r = 10)),
           axis.text.x = element_text (size = 11),
           axis.text.y = element_text (size = 11),
           strip.text.x = element_text (size = 12.5))
  dev.off()

  # plot NNV estimates
  pdf (paste0 (output_folder, "figS7_nnv-siareach.pdf"), width = 9, height = 7)
  ggplot (data = nnv_siareach_senanl [country != "Global"],
          aes(x = country_name, y = nnv)) +
    geom_col (aes(fill = comp_set), width = 0.8, colour = NA, size = 0,
              position = position_dodge()) +
    facet_wrap (vars(siareach), ncol = 1) +
    labs (x = "", y = "Estimated number needed to vaccinate") +
    scale_fill_manual ("Vaccination strategies compared to MCV1 only",
                       values = custom_palette[1:2]) +
    theme_bw () +
    theme (legend.position = "top",
           legend.text = element_text (size = 12),
           legend.title = element_text (size = 13),
           axis.title.y = element_text (size = 13.5, vjust = 2),
           strip.text.x = element_text (size = 12),
           axis.text.x = element_text (size = 12, angle = 60, hjust = 1),
           axis.text.y = element_text (size = 10),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.background = element_blank(),
           plot.margin = margin (0, 0.2, 0, 0.6, "cm"))
  dev.off()


  # ------------------------------------------------------------------------------
  ## plot case and nnv estimates by R0
  # ------------------------------------------------------------------------------
  # load data
  all_r0_senanl <- NULL
  for (ir0 in seq(6,26,2)) {
    if (ir0 %in% c(6,26)) {
      scns <- vac_stgs[c(2,3,4,5)]
    } else {
      scns <- vac_stgs[4]
    }
    for (iscn in scns){
      outputs <- fread (paste0 (res_folder, "siareach_2/senanl_r0/central_burden_estimate_",
                                iscn, "_r0-", ir0, ".csv"))
      all_r0_senanl <- rbind (all_r0_senanl, outputs [, `:=` (r0 = ir0, scenario = iscn)])
      remove(outputs)
    }
  }
  # calculate burden
  all_r0_senanl [, `:=` (country_name = factor (country_names[country], levels = country_names),
                         cases = cases0d + cases1d + cases2d,
                         deaths = deaths0d + deaths1d + deaths2d)]
  cum_r0_senanl <- all_r0_senanl [, lapply (.SD, sum),
                                  .SDcols = c("pops", "cases", "deaths", "dalys", "doses"),
                                  by = c("year", "country", "country_name", "scenario", "r0")]
  setorder (cum_r0_senanl, r0, country, year)

  # calculate NNV and add results of the main scenarios
  nnv_r0_senanl <- NULL
  for (ir0 in c(6,26)){
    output_nnv <- cum_r0_senanl [r0 == ir0 & scenario %in% vac_stgs[c(2,3,5)]]
    output_nnv [, comp := scenario]
    nnv_r0_senanl <- rbind (nnv_r0_senanl,
                            cal_avtnnv (output_nnv, "mcv1", "mcv1-mcv2")[, r0 := ir0],
                            cal_avtnnv (output_nnv, "mcv1", "mcv1-sia")[, r0 := ir0])
    remove (output_nnv)
  }
  nnv_r0_senanl <- rbind (nnv_r0_senanl,
                          all_avtnnv [comp_set %in% c("mcv1-mcv2_VS_mcv1", "mcv1-sia_VS_mcv1")][, r0 := 15.9])
  nnv_r0_senanl [, `:=` (comp_set = factor (comp_set, labels = c("MCV1 + MCV2", "MCV1 + SIAs")),
                         plt_r0 = ifelse (r0 == 26, "R[0] *\" = 26\"",
                                          ifelse (r0 == 6, "R[0] *\" = 6\"",
                                                  "R[0] *\" = 15.9 (main scenario)\"")))]


  # plot case estimates
  pdf (paste0 (output_folder, "figS8_case-r0.pdf"), width = 14, height = 7)
  ggplot (data = cum_r0_senanl[scenario == "mcv1-mcv2-sia"] ) +
    geom_line (size = 0.9, aes (x = year, y = cases/1e3, colour = r0,
                                group = paste0 (country, "-", r0))) +
    scale_x_continuous (breaks = scales::pretty_breaks ()) +
    # scale_y_log10 (labels = format_num_plain) +
    scale_colour_distiller (name = expression(R['0']), palette = "Blues", direction = 1) +
    facet_wrap (vars(country_name), nrow = 3, scales = "free") +
    labs (x = "Year", y = "Cases (thousands)") +
    theme_bw () +
    theme (legend.position = "bottom",
           legend.text = element_text (size = 11.5),
           legend.title = element_text (size = 12),
           plot.margin = unit (c(0.1, 0.25, 0.1, 0.1), "cm"),
           axis.title.x = element_text (size = 15, margin = margin(t = 15)),
           axis.title.y = element_text (size = 15, margin = margin(r = 10)),
           axis.text.x = element_text (size = 11),
           axis.text.y = element_text (size = 11),
           strip.text.x = element_text (size = 12.5)) +
    ggnewscale::new_scale_color () +
    geom_point (data = dat_WHOcase [notifs != 0][, `Data type` := "WHO reported cases"],
                aes(x = year, y = notifs/1e3, shape = `Data type`),
                size = 2, colour = custom_palette[7]) +
    scale_shape_manual ("", values = 2) +
    ggnewscale::new_scale_color () +
    geom_pointrange (data = dat_IHME [, `Data type` := "IHME estimates"],
                     shape = 20, size = 0.4,
                     aes (x = year, y = est_cases/1e3,
                          ymin = lower/1e3, ymax = upper/1e3, colour = `Data type`)) +
    scale_colour_manual (name = " ", values = "#bdbdbd")
  dev.off()
}

# ---------------------------------------------------------
# Save a ggplot figure to file with default settings
# ---------------------------------------------------------
save_fig = function(g, ..., dir = NULL) {
  
  # Collapse inputs into vector of strings
  fig_name_parts = unlist(list(...))
  
  # Construct file name to concatenate with file path
  save_name = paste(fig_name_parts, collapse = " - ")
  
  # Construct path to save file to
  save_path = o$pth$figures
  if (!is.null(dir)) {
    ext_path  = paste(dir, collapse = file_sep())
    save_path = paste0(save_path, ext_path, file_sep())
  }
  
  # Create directory if it exists
  if (!dir.exists(save_path))
    dir.create(save_path)
  
  # Repeat the saving process for each image format in figure_format
  for (fig_format in o$figure_format) {
    full_path = paste0(save_path, save_name, ".", fig_format)
    
    # Save figure (size specified in options.R)
    ggsave(full_path, 
           plot   = g, 
           device = fig_format, 
           dpi    = o$save_resolution, 
           width  = o$save_width, 
           height = o$save_height, 
           units  = o$save_units)
  }
}

