# ============================================
# Comprehensive R Script: CFTC and Macroeconomic Data Processing
# Author: Simon-Pierre Boucher, Gabriel Power, Marie-Hélène Gagnon
# Date: 2024-12-04
# Description: This script processes CFTC data and macroeconomic indicators,
#              performs regression analyses with weighted least squares (KUROV and ANDERSON methods),
#              estimates GARCH models for volatility, and saves the results for further analysis.
# ============================================

# ------------------------------
# 1. Loading Required Libraries
# ------------------------------

# List of required packages
required_packages <- c(
  "Quandl", 
  "tidyverse", 
  "xts", 
  "lubridate", 
  "highfrequency", 
  "tseries", 
  "texreg", 
  "nlme", 
  "stargazer"
)

# Install any missing packages
installed_packages <- rownames(installed.packages())
for(p in required_packages){
  if(!(p %in% installed_packages)){
    install.packages(p, dependencies = TRUE)
  }
}

# Load the libraries
library(Quandl)
library(tidyverse)
library(xts)
library(lubridate)
library(highfrequency)
library(tseries)
library(texreg)
library(nlme)
library(stargazer)

# ------------------------------
# 2. Defining Utility Functions
# ------------------------------

#' Shift Dates by a Given Number of Days
#'
#' @param data A data frame containing a 'Date' column.
#' @param shift_days An integer specifying the number of days to shift. Positive for future, negative for past.
#' @return A data frame with shifted dates.
shift_dates <- function(data, shift_days) {
  data %>%
    mutate(Date = Date + days(shift_days))
}

#' Process CFTC Data for a Given Commodity
#'
#' This function processes either PC (Commitments of Traders) or COT (Commitments of Traders) data for a specified commodity.
#' It applies date shifts and calculates relevant metrics.
#'
#' @param pc_data A data frame containing CFTC PC data.
#' @param shift_values A numeric vector of days to shift the dates for creating lagged/future observations.
#' @param commodity_type A string indicating the type of data ('PC' or 'COT').
#' @return A processed data frame with calculated metrics.
process_cftc_data <- function(pc_data, shift_values, commodity_type) {
  # Initialize with original data
  tt <- pc_data
  
  # Apply date shifts and bind rows
  for (shift in shift_values) {
    shifted_tt <- shift_dates(pc_data, shift)
    tt <- bind_rows(tt, shifted_tt)
  }
  
  # Calculate metrics
  tt <- tt %>%
    mutate(
      MSCT = (`Noncommercial Long` + `Noncommercial Short`) / (2 * `Open Interest`),
      NLS = (`Noncommercial Long` - `Noncommercial Short`) / `Open Interest`,
      WT1 = 1 + (`Noncommercial Short` / (`Commercial Long` + `Commercial Short`)),
      WT2 = 1 + (`Noncommercial Long` / (`Commercial Long` + `Commercial Short`)),
      WT = ifelse(`Commercial Short` > `Commercial Long`, WT1, WT2)
    )
  
  # Additional metrics for COT data
  if (commodity_type == "COT") {
    tt <- tt %>%
      mutate(
        MM_NLS = (`Money Manager Longs` - `Money Manager Shorts`) / `Open Interest`,
        SWAP_NLS = (`Swap Dealer Longs` - `Swap Dealer Shorts`) / `Open Interest`
      )
  }
  
  return(tt)
}

#' Process 5-Minute Commodity Data
#'
#' This function reads 5-minute interval data for a commodity, calculates returns, and merges shifted data.
#'
#' @param file_path A string specifying the path to the 5-minute CSV data file.
#' @param shift_minutes An integer specifying the number of minutes to shift the time for creating lagged observations.
#' @return A merged data frame with original and shifted data.
process_commodity_data <- function(file_path, shift_minutes = 5) {
  # Read CSV data
  data <- read.csv(file_path, header = FALSE, sep = ";") %>%
    rename(
      Date = V1,
      Time = V2,
      Open = V3,
      High = V4,
      Low = V5,
      Close = V6,
      Volume = V7
    ) %>%
    mutate(
      DAY = as.numeric(substr(Date, 1, 2)),
      MONTH = as.numeric(substr(Date, 4, 5)),
      YEAR = as.numeric(substr(Date, 7, 10)),
      HOUR = as.numeric(substr(Time, 1, 2)),
      MINUTE = as.numeric(substr(Time, 4, 5)),
      r_t1 = (Close - Open) / Open,
      t = HOUR * 60 + MINUTE
    )
  
  # Create shifted time
  data_shifted <- data %>%
    mutate(t = t + shift_minutes) %>%
    select(YEAR, MONTH, DAY, t, r_t1)
  
  # Merge original and shifted data
  data_merged <- merge(data, data_shifted, by = c("YEAR", "MONTH", "DAY", "t"), all.x = TRUE)
  
  return(data_merged)
}

#' Process Macroeconomic Data
#'
#' This function reads multiple macroeconomic CSV files, assigns categories, calculates surprises, and standardizes them.
#'
#' @param macro_files_path A string specifying the directory path containing macroeconomic CSV files.
#' @return A processed data frame with standardized surprises and temporal variables.
process_macro_data <- function(macro_files_path) {
  # List all macroeconomic CSV files
  macro_files <- list.files(macro_files_path, pattern = "^data_\\d{4}\\.csv$", full.names = TRUE)
  
  # Read and combine all macroeconomic data
  macro_data <- macro_files %>%
    map_dfr(~ read_csv(.x, 
                       col_types = cols(
                         date = col_date(format = "%d/%m/%Y"),
                         time = col_time(format = "%H:%M"),
                         actual = col_double(),
                         forecast = col_double(),
                         previous = col_double()
                       )))
  
  # Define event categories
  events <- list(
    IJC = "Initial Jobless Claims",
    ADP = "ADP Nonfarm Employment Change",
    CBCC = "CB Consumer Confidence",
    ARS = "Retail Sales (MoM)",
    BP = "Building Permits",
    ConstS = "Construction Spending (MoM)",
    ConsC = "Consumer Credit",
    CPI = "CPI (MoM)",
    DGO = "Durable Goods Orders (MoM)",
    EHS = "Existing Home Sales",
    FO = "Factory Orders (MoM)",
    GDP = "GDP (QoQ)",
    HS = "Housing Starts",
    IP = "Industrial Production (MoM)",
    MCS = "Michigan Consumer Sentiment",
    NHS = "New Home Sales",
    NFPR = "Nonfarm Payrolls",
    PHS = "Pending Home Sales (MoM)",
    PersoS = "Personal Spending (MoM)",
    PersoI = "Personal Income (MoM)",
    PPI = "PPI (MoM)",
    TB = "Trade Balance",
    NGI = "Natural Gas Storage",
    AWCOS = "API Weekly Crude Oil Stock"
  )
  
  # Assign MNA categories based on event names
  macro_filtered <- macro_data %>%
    mutate(
      MNA = case_when(
        str_detect(event, paste0("^", events$IJC, ".*")) ~ "IJC",
        str_detect(event, paste0("^", events$ADP, ".*")) ~ "ADP",
        str_detect(event, paste0("^", events$CBCC, ".*")) ~ "CBCC",
        str_detect(event, paste0("^", events$ARS, ".*")) ~ "ARS",
        str_detect(event, paste0("^", events$BP, ".*")) ~ "BP",
        str_detect(event, paste0("^", events$ConstS, ".*")) ~ "ConstS",
        str_detect(event, paste0("^", events$ConsC, ".*")) ~ "ConsC",
        str_detect(event, paste0("^", events$CPI, ".*")) ~ "CPI",
        str_detect(event, paste0("^", events$DGO, ".*")) ~ "DGO",
        str_detect(event, paste0("^", events$EHS, ".*")) ~ "EHS",
        str_detect(event, paste0("^", events$FO, ".*")) ~ "FO",
        str_detect(event, paste0("^", events$GDP, ".*")) ~ "GDP",
        str_detect(event, paste0("^", events$HS, ".*")) ~ "HS",
        str_detect(event, paste0("^", events$IP, ".*")) ~ "IP",
        str_detect(event, paste0("^", events$MCS, ".*")) ~ "MCS",
        str_detect(event, paste0("^", events$NHS, ".*")) ~ "NHS",
        str_detect(event, paste0("^", events$NFPR, ".*")) ~ "NFPR",
        str_detect(event, paste0("^", events$PHS, ".*")) ~ "PHS",
        str_detect(event, paste0("^", events$PersoS, ".*")) ~ "PersoS",
        str_detect(event, paste0("^", events$PersoI, ".*")) ~ "PersoI",
        str_detect(event, paste0("^", events$PPI, ".*")) ~ "PPI",
        str_detect(event, paste0("^", events$TB, ".*")) ~ "TB",
        str_detect(event, paste0("^", events$NGI, ".*")) ~ "NGI",
        str_detect(event, paste0("^", events$AWCOS, ".*")) ~ "AWCOS",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(MNA)) %>%
    mutate(
      SURPRISE = actual - forecast,
      DATE = as.Date(paste(year(date), month(date), day(date), sep = "-"), "%Y-%m-%d"),
      HOUR = hour(time),
      MINUTE = minute(time),
      t = HOUR * 60 + MINUTE
    )
  
  # Calculate standard deviation per MNA category
  ecart_type_par_categorie <- macro_filtered %>%
    group_by(MNA) %>%
    summarise(ecart_type = sd(SURPRISE, na.rm = TRUE))
  
  # Merge and standardize surprises
  macro_final <- macro_filtered %>%
    left_join(ecart_type_par_categorie, by = "MNA") %>%
    mutate(
      SS = SURPRISE / ecart_type
    )
  
  # Create standardized surprise columns per MNA category
  unique_cats <- unique(macro_final$MNA)
  for(cat in unique_cats){
    var_name <- paste0("SS", cat)
    macro_final[[var_name]] <- ifelse(macro_final$MNA == cat, macro_final$SS, 0)
    macro_final[[var_name]] <- as.numeric(macro_final[[var_name]])
  }
  
  # Remove rows with any missing values
  macro_final <- macro_final %>%
    drop_na()
  
  # Create additional temporal variables
  macro_final <- macro_final %>%
    mutate(
      YEAR = year(DATE),
      MONTH = month(DATE),
      DAY = day(DATE),
      TIMEF = factor(HOUR),
      t = HOUR * 60 + MINUTE
    )
  
  # Save the processed macroeconomic data
  save(macro_final, file = "MNA_DATA_ALL.Rdata")
  
  return(macro_final)
}

#' Merge Commodity Data with Macroeconomic and COT Data
#'
#' This function merges 5-minute commodity data with macroeconomic indicators and COT data.
#'
#' @param commodity_symbol A string representing the commodity symbol (e.g., "CL" for Crude Oil).
#' @param cot_data_file A string specifying the path to the COT data file.
#' @param pc_data_db A data frame containing the processed PC data for the commodity.
#' @return A merged data frame containing commodity, macroeconomic, and COT data.
merge_commodity_data <- function(commodity_symbol, cot_data_file, pc_data_db) {
  # Define the file path for 5-minute data
  file_path <- paste0(tolower(commodity_symbol), "-5m.csv")
  commodity_data <- process_commodity_data(file_path)
  
  # Load macroeconomic data
  load("MNA_DATA_ALL.Rdata")
  
  # Merge with macroeconomic data on YEAR, MONTH, DAY, and t
  commodity_merged <- merge(commodity_data, macro_final, by = c("YEAR", "MONTH", "DAY", "t"), all.x = TRUE)
  
  # Load COT data
  load(cot_data_file)
  
  # Prepare COT data for merging
  cot_data_db <- get(gsub(".Rda$", "", cot_data_file)) %>%
    mutate(
      YEAR = as.numeric(substr(Date, 1, 4)),
      MONTH = as.numeric(substr(Date, 6, 7)),
      DAY = as.numeric(substr(Date, 9, 10))
    )
  
  # Merge with COT data on YEAR, MONTH, DAY
  commodity_final <- merge(commodity_merged, cot_data_db, by = c("YEAR", "MONTH", "DAY"), all.x = TRUE)
  
  # Replace missing values with 0
  commodity_final[is.na(commodity_final)] <- 0
  
  # Sort the data
  commodity_final <- commodity_final %>%
    arrange(YEAR, MONTH, DAY, t)
  
  # Save the merged data
  save(commodity_final, file = paste0(commodity_symbol, "_PC_DB.Rdata"))
  
  return(commodity_final)
}

#' Create Data Subsets Based on Specific Periods
#'
#' This function filters the data for specified date ranges and saves the subsets.
#'
#' @param db A data frame containing the merged data.
#' @param period A string indicating the period name (e.g., "ZLB", "COVID").
#' @param start_date A string representing the start date in "YYYY-MM-DD" format.
#' @param end_date A string representing the end date in "YYYY-MM-DD" format.
#' @param suffix A string to append to the subset data frame's name.
#' @return None. The subset is saved as an RData file.
create_subsets <- function(db, period, start_date, end_date, suffix) {
  # Create a DATE column if not already present
  if(!"DATE" %in% colnames(db)){
    db <- db %>%
      mutate(DATE = as.Date(paste(YEAR, MONTH, DAY, sep = "-"), "%Y-%m-%d"))
  }
  
  # Filter the data for the specified period
  subset_db <- db %>%
    filter(DATE >= as.Date(start_date) & DATE <= as.Date(end_date))
  
  # Generate the new variable name
  new_var_name <- paste0(deparse(substitute(db)), "_", suffix)
  
  # Assign the subset to a new variable in the global environment
  assign(new_var_name, subset_db, envir = .GlobalEnv)
  
  # Save the subset data frame
  save(list = new_var_name, file = paste0(new_var_name, ".Rdata"))
}

# ------------------------------
# 3. Downloading and Processing CFTC Data
# ------------------------------

# Set the working directory
setwd("/Volumes/LaCie/PHD/CH1")

# Set Quandl API key
Quandl.api_key("e9aZeuvn-X6isXXx9TCk")

# Define CFTC data codes for different commodities
cftc_codes <- list(
  GC_COT = "CFTC/088691_F_ALL",
  GC_PC = "CFTC/088691_F_L_ALL",
  CL_COT = "CFTC/067651_F_ALL",
  CL_PC = "CFTC/067651_F_L_ALL",
  HG_COT = "CFTC/085692_F_ALL",
  HG_PC = "CFTC/085692_F_L_ALL",
  SI_COT = "CFTC/084691_F_ALL",
  SI_PC = "CFTC/084691_F_L_ALL",
  NG_COT = "CFTC/023651_F_ALL",
  NG_PC = "CFTC/023651_F_L_ALL",
  PA_COT = "CFTC/075651_F_ALL",
  PA_PC = "CFTC/075651_F_L_ALL"
)

# Download and assign CFTC data
for(name in names(cftc_codes)){
  assign(name, Quandl(cftc_codes[[name]]))
}

# Define shift values in days
shift_values_days <- c(-3, -2, -1, 1, 2, 3)

# List of commodities for PC and COT data
commodities_pc <- c("CL", "GC", "SI", "HG", "PA", "NG")
commodities_cot <- c("CL", "GC", "SI", "HG", "PA", "NG")

# Process and save PC data for each commodity
for(commodity in commodities_pc){
  pc_data <- get(paste0("CFTC_", commodity, "_PC"))
  processed_pc <- process_cftc_data(pc_data, shift_values_days, "PC")
  assign(paste0(commodity, "_PC"), processed_pc)
  
  # Save the processed PC data
  save(list = paste0(commodity, "_PC"), file = paste0(commodity, "_PC.Rda"))
}

# Process and save COT data for each commodity
for(commodity in commodities_cot){
  cot_data <- get(paste0("CFTC_", commodity, "_COT"))
  processed_cot <- process_cftc_data(cot_data, shift_values_days, "COT")
  assign(paste0(commodity, "_COT"), processed_cot)
  
  # Save the processed COT data
  save(list = paste0(commodity, "_COT"), file = paste0(commodity, "_COT.Rda"))
}

# ------------------------------
# 4. Processing Macroeconomic Data
# ------------------------------

# Define the path to macroeconomic CSV files
macro_files_path <- "/Volumes/LaCie/PHD/CH1"

# Process macroeconomic data
macro_final <- process_macro_data(macro_files_path)

# ------------------------------
# 5. Merging Commodity Data with Macroeconomic and COT Data
# ------------------------------

# List of all commodities
commodities <- c("CL", "GC", "SI", "HG", "PA", "NG")

# Initialize a list to store merged data frames
commodity_dbs <- list()

for(commodity in commodities){
  # Define the file paths
  cot_file <- paste0(commodity, "_COT.Rda")
  pc_db_file <- paste0(commodity, "_PC.Rda")
  
  # Merge data
  merged_db <- merge_commodity_data(
    commodity_symbol = commodity,
    cot_data_file = cot_file,
    pc_data_db = get(paste0(commodity, "_PC"))
  )
  
  # Store the merged data in the list
  commodity_dbs[[commodity]] <- merged_db
}

# ------------------------------
# 6. Creating Data Subsets for Specific Periods
# ------------------------------

# Define specific periods with start and end dates
periods <- list(
  ZLB = list(start_date = "2008-12-22", end_date = "2015-12-21"),
  COVID = list(start_date = "2020-01-31", end_date = "2022-06-10"),
  SUBPRIME = list(start_date = "2007-08-01", end_date = "2009-06-30"),
  YEAR_2011_2019 = list(start_date = "2011-01-01", end_date = "2019-12-31")
)

# Create and save subsets for each commodity and period
for(commodity in commodities){
  db <- commodity_dbs[[commodity]]
  
  for(period_name in names(periods)){
    create_subsets(
      db = db,
      period = period_name,
      start_date = periods[[period_name]]$start_date,
      end_date = periods[[period_name]]$end_date,
      suffix = period_name
    )
  }
}

# ------------------------------
# 7. Defining Custom Regression and GARCH Functions
# ------------------------------

#' Custom Linear Regression using Maximum Likelihood Estimation
#'
#' @param formula An object of class "formula" (or one that can be coerced to that class).
#' @param data A data frame containing the variables in the model.
#' @param weights An optional vector of weights to be used in the fitting process.
#' @param initial_params Optional initial parameters for optimization.
#' @return A list containing estimated coefficients, sigma, residuals, log-likelihood, convergence status, and message.
custom_lm_mle <- function(formula, data, weights = NULL, initial_params = NULL) {
  model_matrix <- model.matrix(formula, data)
  response <- model.response(model.frame(formula, data))
  n <- length(response)
  
  # Negative Log-Likelihood Function
  neg_log_likelihood <- function(params) {
    beta <- params[1:(length(params) - 1)]
    log_sigma <- params[length(params)]
    sigma <- exp(log_sigma)
    predictions <- model_matrix %*% beta
    residuals <- response - predictions
    
    if (!is.null(weights)) {
      residuals <- residuals * sqrt(weights)
      sigma <- sigma / sqrt(mean(weights))
    }
    
    nll <- 0.5 * n * log(2 * pi) + n * log(sigma) + (1 / (2 * sigma^2)) * sum(residuals^2)
    return(nll)
  }
  
  # Initial parameter estimates
  if (is.null(initial_params)) {
    initial_beta <- rep(0, ncol(model_matrix))
    initial_log_sigma <- log(sd(response))
    initial_params <- c(initial_beta, initial_log_sigma)
  }
  
  # Optimization using BFGS method
  optim_result <- optim(
    par = initial_params, 
    fn = neg_log_likelihood, 
    method = "BFGS", 
    control = list(maxit = 1000)
  )
  
  # Check for convergence
  if (optim_result$convergence != 0) {
    warning("Optimization did not converge.")
  }
  
  # Extract estimates
  beta_est <- optim_result$par[1:(length(optim_result$par) - 1)]
  log_sigma_est <- optim_result$par[length(optim_result$par)]
  sigma_est <- exp(log_sigma_est)
  
  # Calculate residuals
  predictions <- model_matrix %*% beta_est
  residuals <- response - predictions
  
  # Return results as a list
  list(
    coefficients = beta_est,
    sigma = sigma_est,
    residuals = residuals,
    logLik = -optim_result$value,
    convergence = optim_result$convergence,
    message = optim_result$message
  )
}

#' Calculate Weights for Weighted Least Squares (KUROV Method)
#'
#' @param residuals A numeric vector of residuals from an initial regression.
#' @return A numeric vector of weights.
calculate_weights <- function(residuals) {
  w <- numeric(length = length(residuals) + 1)
  w[1] <- abs(residuals[1])
  for (i in 1:(length(residuals)-1)) {
    w[i + 1] <- 0.1 * abs(residuals[i]) + 0.9 * w[i]
  }
  w <- w[1:length(residuals)]
  weights <- 1 / (w^2)
  return(weights)
}

#' Custom GARCH(1,1) Model Estimation using Maximum Likelihood
#'
#' @param residuals A numeric vector of residuals from a regression.
#' @param initial_params Optional initial parameters for optimization (omega, alpha, beta).
#' @return A list containing estimated GARCH parameters, sigma squared, log-likelihood, convergence status, and message.
custom_garch_mle <- function(residuals, initial_params = NULL) {
  n <- length(residuals)
  
  # Negative Log-Likelihood for GARCH(1,1)
  neg_log_likelihood_garch <- function(params) {
    omega <- params[1]
    alpha <- params[2]
    beta <- params[3]
    
    # Constraints to ensure positivity and stationarity
    if (omega <= 0 || alpha < 0 || beta < 0 || (alpha + beta) >= 1) {
      return(1e10)  # Penalize invalid parameters
    }
    
    sigma2 <- numeric(n)
    sigma2[1] <- var(residuals)  # Initialize with variance of residuals
    
    # Recursively calculate sigma squared
    for (t in 2:n) {
      sigma2[t] <- omega + alpha * residuals[t - 1]^2 + beta * sigma2[t - 1]
    }
    
    # Calculate log-likelihood
    log_likelihood <- -0.5 * sum(log(2 * pi * sigma2) + (residuals^2) / sigma2)
    return(-log_likelihood)  # Return negative log-likelihood for minimization
  }
  
  # Initial parameter estimates
  if (is.null(initial_params)) {
    omega_init <- 0.1 * var(residuals)
    alpha_init <- 0.05
    beta_init <- 0.90
    initial_params <- c(omega_init, alpha_init, beta_init)
  }
  
  # Optimization using L-BFGS-B method with bounds
  optim_result <- optim(
    par = initial_params, 
    fn = neg_log_likelihood_garch, 
    method = "L-BFGS-B",
    lower = c(1e-6, 0, 0),
    upper = c(var(residuals), 1, 1),
    control = list(maxit = 1000)
  )
  
  # Check for convergence
  if (optim_result$convergence != 0) {
    warning("GARCH optimization did not converge.")
  }
  
  # Extract estimates
  omega_est <- optim_result$par[1]
  alpha_est <- optim_result$par[2]
  beta_est <- optim_result$par[3]
  
  # Calculate sigma squared recursively
  sigma2 <- numeric(n)
  sigma2[1] <- var(residuals)
  for (t in 2:n) {
    sigma2[t] <- omega_est + alpha_est * residuals[t - 1]^2 + beta_est * sigma2[t - 1]
  }
  
  # Return results as a list
  list(
    omega = omega_est,
    alpha = alpha_est,
    beta = beta_est,
    sigma2 = sigma2,
    logLik = -optim_result$value,
    convergence = optim_result$convergence,
    message = optim_result$message
  )
}

# ------------------------------
# 8. Defining Regression Functions
# ------------------------------

#' Perform Regression Analyses (KUROV and ANDERSON Methods)
#'
#' This function performs regression analyses using both KUROV (weighted least squares) and ANDERSON methods,
#' and estimates GARCH models for volatility. It saves each model as an RDS file.
#'
#' @param data A data frame containing the necessary variables for regression.
#' @param commodity A string representing the commodity symbol (e.g., "SI").
#' @param var A string representing the explanatory variable (e.g., "MM_NLS" or "SWAP_NLS").
#' @return None. Models are saved as RDS files.
perform_regressions <- function(data, commodity, var){
  # Define the main regression formula
  formula_main <- as.formula(paste0(
    "r_t1 ~ (SSIJC + SSADP + SSCBCC + SSARS + SSBP + SSConstS + SSConsC + SSCPI + SSDGO + ",
    "SSEHS + SSFO + SSGDP + SSHS + SSIP + SSMCS + SSNHS + SSNFPR + SSPHS + SSPersoS + ",
    "SSPersoI + SSPPI + SSTB) * ", var, " + r_m1"
  ))
  
  # Initial regression without weights
  initial_regress <- custom_lm_mle(formula_main, data = data)
  
  # Calculate weights using KUROV method
  regress_weights_KUROV <- calculate_weights(initial_regress$residuals)
  
  # Weighted regression (KUROV)
  regress_KUROV <- lm(formula_main, data = data, weights = regress_weights_KUROV)
  model_KUROV_name <- paste0("regress_", commodity, "_", var, "_KUROV")
  assign(model_KUROV_name, regress_KUROV, envir = .GlobalEnv)
  
  # Save KUROV model
  saveRDS(regress_KUROV, paste0(model_KUROV_name, ".rds"))
  
  # Prepare for ANDERSON regression by adding absolute residuals
  data$abs_residuals <- abs(initial_regress$residuals)
  
  # Define the Anderson regression formula
  formula_anderson <- as.formula(paste0(
    "abs_residuals ~ (SSIJC + SSADP + SSCBCC + SSARS + SSBP + SSConstS + SSConsC + SSCPI + SSDGO + ",
    "SSEHS + SSFO + SSGDP + SSHS + SSIP + SSMCS + SSNHS + SSNFPR + SSPHS + SSPersoS + ",
    "SSPersoI + SSPPI + SSTB) * ", var, " + TIMEF"
  ))
  
  # Anderson regression
  regress_ANDERSON <- lm(formula_anderson, data = data, weights = regress_weights_KUROV)
  model_ANDERSON_name <- paste0("regress_", commodity, "_", var, "_ANDERSON")
  assign(model_ANDERSON_name, regress_ANDERSON, envir = .GlobalEnv)
  
  # Save ANDERSON model
  saveRDS(regress_ANDERSON, paste0(model_ANDERSON_name, ".rds"))
  
  # GARCH model estimation on KUROV residuals
  garch_model <- custom_garch_mle(regress_KUROV$residuals)
  model_VOL_name <- paste0("regress_", commodity, "_", var, "_VOL")
  assign(model_VOL_name, garch_model, envir = .GlobalEnv)
  
  # Save VOL model
  saveRDS(garch_model, paste0(model_VOL_name, ".rds"))
}

# ------------------------------
# 9. Processing Specific Commodity Data (e.g., SI_COT_DB)
# ------------------------------

# Set working directory for specific commodity processing
setwd("/Volumes/LaCie/PHD/CH1")

# Load the specific commodity COT data
load("SI_COT_DB.Rdata")
X <- SI_COT_DB

# Prepare explanatory variables with 'SS' prefix
X <- X %>%
  mutate(
    SSIJC = SS_IJC,
    SSADP = SS_ADP,
    SSCBCC = SS_CBCC,
    SSARS = SS_ARS,
    SSBP = SS_BP,
    SSConstS = SS_ConstS,
    SSConsC = SS_ConsC,
    SSCPI = SS_CPI,
    SSDGO = SS_DGO,
    SSEHS = SS_EHS,
    SSFO = SS_FO,
    SSGDP = SS_GDP,
    SSHS = SS_HS,
    SSIP = SS_IP,
    SSMCS = SS_MCS,
    SSNHS = SS_NHS,
    SSNFPR = SS_NFPR,
    SSPHS = SS_PHS,
    SSPersoS = SS_PersoS,
    SSPersoI = SS_PersoI,
    SSPPI = SS_PPI,
    SSTB = SS_TB
  )

# Create binary variables S1 to S22 based on SS_* variables
binary_vars <- paste0("S", 1:22)
events_list <- events_names()

for(i in 1:22){
  var_name <- paste0("S", i)
  ss_var <- paste0("SS", events_list[i])
  X[[var_name]] <- ifelse(X[[ss_var]] == 0, 0, 1)
}

# Perform regressions for MM_NLS and SWAP_NLS using the custom regression functions
perform_regressions(X, "SI", "MM_NLS")
perform_regressions(X, "SI", "SWAP_NLS")


