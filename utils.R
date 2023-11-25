# utils.R
# Functions for supporting utilities in the DynaMICE model (Dynamic Measles Immunisation Calculation Engine)

# ------------------------------------------------------------------------------
#' Create vaccine coverage files by scenarios
#'
#' Generates files for routine programmes (MCV1, MCV2) and supplementary
#' immunisation activities (SIAs) from the vaccine coverage files.
# ------------------------------------------------------------------------------
#' @param scenario_name Name of the vaccination scenario selected or being analysed.
#'
#' @import data.table
#'
#' @examples
#' \dontrun{
#'   create_vaccine_coverage_routine_sia("campaign-only-bestcase")
#'   }
create_vaccine_coverage_routine_sia = function(scenario_name) {
  
  # Coverage file for this scenario - to be split between routine and SIA
  file = list(
    coverage = paste0(o$pth$coverage, "coverage_", scenario_name, ".csv"), 
    routine  = paste0(o$pth$coverage, "routine_",  scenario_name, ".csv"), 
    sia      = paste0(o$pth$coverage, "sia_",      scenario_name, ".csv"))
  
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


# ------------------------------------------------------------------------------
#' Expand the matrix to a different dimension
#
#' This function returns an expanded contact matrix with the specified age
#' structure for model inputs.
# ------------------------------------------------------------------------------
#' @param A A matrix to be expanded.
#' @param expand_rows Number of times to repeat each row.
#' @param expand_cols Number of times to repeat each column.
#' @param rescale_rows A logical variable to control whether to re-scale the
#' expanded rows.
#' @param rescale_cols A logical variable to control whether to re-scale the
#' expanded columns.
#' @return An expanded matrix with re-scaling if applicable.
#' @examples
#' expandMatrix (matrix(1:9,3,3), 2, 1, FALSE, FALSE)
expandMatrix <- function (A,
                          expand_rows  = 1,
                          expand_cols  = 1,
                          rescale_rows = F,
                          rescale_cols = F) {
  
  if(!is.matrix(A)){
    stop("A is not a matrix")
  }
  
  matvals <- numeric(0)
  rows <- nrow(A)
  cols <- ncol(A)
  
  for(c in 1:cols) {
    matvals <- c(
      matvals,
      rep(
        A[,c],
        expand_cols,
        each = expand_rows
      )
    )
  }
  
  B <- matrix (matvals,
               nrow = rows * expand_rows,
               ncol = cols * expand_cols)
  
  if(rescale_rows & rescale_cols){
    B <- B/(expand_rows*expand_cols)
  } else if(rescale_rows){
    B <- B/expand_rows
  } else if(rescale_cols){
    B <- B/expand_cols
  }
  
  return (B)
  
} # end of function -- expandMatrix
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
#' Tailor the data structure for life expectancy by age and year
#'
#' Tailor the \code{\link{data_lexp_remain}} data to the format for processing burden
#' estimates. Linear interpolation between calender years was applied.
# ------------------------------------------------------------------------------
#' @param sel_countries ISO-3 codes for countries included for evaluation. If
#' "all", all countries in the original data are selected.
#'
#' @import data.table
#'
#' @examples
#' lexp_remain <- tailor_data_lexp_remain (sel_countries = c("AGO","BGD"))
tailor_data_lexp_remain <- function (sel_countries = "all"){
  
  lexp_remain <- setDT (data_lexp_remain)
  lexp_remain <- lexp_remain [year >= 1980]
  
  if (sel_countries[[1]] != "all") {
    lexp_remain <- lexp_remain [country_code %in% sel_countries]
  }
  
  # copy values for year 2095 to year 2100
  lexp_remain <- rbind (lexp_remain, copy (lexp_remain [year == 2095])[, year := 2100])
  
  # calculate the difference between years
  setorder (lexp_remain, country_code, age_from, year)
  lexp_remain [ , diffy := value - shift(value) , by = .(country_code, age_from)]
  
  # interpolate a linear trend for between-years
  lexp_remain_yr <- copy (lexp_remain) [year == 1980]
  for (btwyr in 0:4) {
    dt <- copy (lexp_remain) [year != 1980]
    dt [, `:=` (year = year - btwyr,
                value = value - btwyr*(diffy/5))]
    
    lexp_remain_yr <- rbind (lexp_remain_yr, dt)
  }
  
  # interpolate a linear trend for between-ages
  setorder (lexp_remain_yr, country_code, year, age_from)
  lexp_remain_yr [ , age_mean := (age_from + age_to)/2]
  
  lexp_remain_full <- lexp_remain_yr [, .(value = approx(age_mean, value, xout = 0:100)$y),
                                      by = .(country_code, country, year)]
  lexp_remain_full [, age := rep(0:100, length(unique(year))*length(unique(country_code)))]
  
  return (lexp_remain_full)
  
} # end of function -- tailor_data_lexp_remain
# ------------------------------------------------------------------------------

