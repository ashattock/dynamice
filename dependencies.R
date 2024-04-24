###########################################################
# DEPENDENCIES
#
# Deal with all package dependencies in one place.
#
###########################################################

# ---- R version check ----

# R versions for which this project has been tested and is stable
stable_versions = "4.3.0"

# R versions for which this project is stable (as a string)
stable_str = paste(stable_versions, collapse = ", ")

# Get details of R version currently running
version_info = R.Version()

# Construct version number from list details
version_num = paste0(version_info$major, ".",  version_info$minor)

# Throw an error if this R version is unsuitable
if (!version_num %in% stable_versions)
  stop("This software is stable with R version(s): ", stable_str,
       " (currently running ", version_num, ")")

# Clear global environment
rm(list = ls())

# ---- Define packages ----

message("* Installing required packages")

# Complete list of all R packages required for this project
packages = c(
  "tidyverse",      # Includes ggplot2, dplyr, tidyr (www.tidyverse.org/packages/)
  "data.table",     # Next generation dataframes
  "Rcpp",           # C++ integration
  "magrittr",       # Additional pipe operators, such as %<>%
  "wrapr",          # Convenience functions (eg qc)
  "countrycode",    # Country name <-> code transformation
  "readxl",         # Data loading functionality
  "lubridate",      # Data formatting functionality
  "naniar",         # Data formatting functionality
  "ggpubr",         # Plotting functionality
  "scales",         # Plotting functionality
  "pals",           # Colour palettes
  "colorspace")     # Colour palettes

# ---- Install and/or load R packages with pacman ----

# Check whether pacman itself has been installed
pacman_installed = "pacman" %in% rownames(installed.packages())

# If not, install it
if (!pacman_installed)
  install.packages("pacman")

# Load pacman
library(pacman)

# Load all required packages, installing them if required
pacman::p_load(char = packages)

# ---- Source files ----

# Scripts that should bot be sourced
no_src = c("launch.R", "dependencies.R", "submit.R")

# All R files, and those to source
all_files = list.files(pattern = ".+\\.R$")
src_files = setdiff(all_files, no_src)

# Source each of these files
for (file in src_files)
  source(file)

# Source all C++ files
all_cpp_files = list.files(pattern = ".+\\.cpp$")
for (file in all_cpp_files)
  sourceCpp(file)

# ---- Redefine or unmask particular functions ----

# Unmask certain commonly used functions
select  = dplyr::select
filter  = dplyr::filter
rename  = dplyr::rename
recode  = dplyr::recode
count   = dplyr::count
predict = stats::predict

# ---- Clean up ----

# Tidy up console after package loading
if (interactive()) clf()  # Close figures
if (interactive()) clc()  # Clear console

