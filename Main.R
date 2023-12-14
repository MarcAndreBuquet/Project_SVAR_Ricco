
# _______________________________________________________________________   
#' @Title : Code SVAR Weale Wieladek
#' @author: Marc-André Buquet
#' @date: 26/11/2023
#  Mail: marc-andre.buquet@ensae.fr  
# _______________________________________________________________________       


## clean environment

rm(list = ls())
setwd("C:/Users/mabuq/Documents/M2_ENSAE/VAR-LP/Projet")
set.seed(42)

## Library

library(remotes)
library(minqa)
library(mvnfast)
library(HI) #install directly from the R archives as not available anymore from CRAN direct installation 

# remotes::install_github("https://github.com/cran/VARsignR/tree/master")
# remotes::install_github("roootra/ZerosignR")
# 
# library(VARsignR)
# library(ZerosignR)
# library(BVAR)



## Source functions from external repositories in different files

#source("C:/Users/mabuq/Documents/M2_ENSAE/VAR-LP/Projet/function_package_SVAR_sign.R")
source("C:/Users/mabuq/Documents/M2_ENSAE/VAR-LP/Projet/code_R/Fonctions_packages.R")

## Retrieve custom functions to make the code run

source("C:/Users/mabuq/Documents/M2_ENSAE/VAR-LP/Projet/All_functions_VARsignR_modified.R")
source("C:/Users/mabuq/Documents/M2_ENSAE/VAR-LP/Projet/RWZaccept_modified.R")
source("C:/Users/mabuq/Documents/M2_ENSAE/VAR-LP/Projet/RWZreject_modified.R")
source("C:/Users/mabuq/Documents/M2_ENSAE/VAR-LP/Projet/plot_figure_new.R")

## Preliminary info for plot computation

var_names = c("Real GDP (in log)", "CPI (in log)", "Asset Purchases", "Long rate", "Real Equity Prices")
IRFs_horizon = 36
Confidence_bands = c(16,84)
Pspar = 0.675 # smoothness parameter for IRF. 0= no smoothing, 1= lots of smoothing, originally = 0.675
SL = 11 # letter sizes ggplots
WW = 14 # width graphs
HH = 11 # height graphs


## Data import and first cleaning

Datastream_set = rio::import("C:/Users/mabuq/Documents/M2_ENSAE/VAR-LP/Projet/Data_Time_Series_Weale_Wieladeck.xlsx") %>% 
  slice(-c(1:4)) %>%
  lapply(., as.numeric) %>%
  as.data.frame() %>%
  dplyr::rename(Period = `Variable.code`) %>%
  mutate(Period_date = base::as.Date(as.numeric(Period), origin = "1899-12-30"))

## Some transformations to match transformations operated in the original paper

Datastream_set_transformed = create_log_vars(Datastream_set, vars = c("CPI_US","CPI_UK","Real_GDP_US","Real_GDP_UK","Real_Equity_Prices_US","Real_Equity_Prices_UK","HHUNC_US"))

### First implementation of code with the first identification scheme (standard Cholesky decomposition) for US -----

VAR_data_1st = Datastream_set_transformed %>% select(CPI_US_log, Real_GDP_US_log, Asset_purchases_US, Yield_US_10Y , Real_Equity_Prices_US_log  )

# plot(x = 1: length(VAR_data_1st$CPI_US_log), y =  VAR_data_1st$CPI_US_log) # possibly trend
# plot(x = 1: length(VAR_data_1st$CPI_US_log), y =  VAR_data_1st$Real_GDP_US_log) # possibly trend
# plot(x = 1: length(VAR_data_1st$CPI_US_log), y =  VAR_data_1st$Asset_purchases_US)
# plot(x = 1: length(VAR_data_1st$CPI_US_log), y =  VAR_data_1st$Yield_US_10Y)
# plot(x = 1: length(VAR_data_1st$CPI_US_log), y =  VAR_data_1st$Real_Equity_Prices_US_log)

## Compute the RF VAR first

VARselect(VAR_data_1st, lag.max = 12, type = c("both"))
          

RF_VAR = vars::VAR(VAR_data_1st, type = "trend" )

RF_VAR

# Setup the Cholesky matrix with ordering as done in their paper

Cholesky_matrix = diag(1,5)
Cholesky_matrix[lower.tri(Cholesky_matrix)] <- NA

SVAR_Cholesky = SVAR(RF_VAR, Bmat = Cholesky_matrix )

IRFs_Cholesky_US = vars::irf(SVAR_Cholesky, ci = 0.68) # IRFs computed for Cholesky for US (do the same for UK)

IRFs_Cholesky_median = IRFs_Cholesky_US$irf$Asset_purchases_US[,3]
IRFs_Cholesky_lower = IRFs_Cholesky_US$Lower$Asset_purchases_US[,3]
IRFs_Cholesky_upper = IRFs_Cholesky_US$Upper$Asset_purchases_US[,3]

### Compute SVAR for the the 3 remaining schemes ------




## Implementation 

# Data

SVAR_Fig_2_data = Datastream_set_transformed %>% 
  filter(Period > 39845) %>%
  select(CPI_US_log, Real_GDP_US_log, Asset_purchases_US, Yield_US_10Y, Real_Equity_Prices_US) %>%
  as.matrix()

## Specify sign restrictions

# Second identification scheme

Sign_restrictions_2nd_scheme = list(matrix(c(-1, 1, NA, NA, NA, 1, 1, NA, NA, NA, NA, NA, 1 , NA, NA, 1 , 1 , -1, NA, NA, 1 , 1 , 1 , NA, NA), nrow = 5), # Period 0
                                    matrix(c(-1, 1, NA, NA, NA, 1, 1, NA, NA, NA, NA, NA, 1 , NA, NA, 1 , 1 , -1, NA, NA, 1 , 1 , 1 , NA, NA), nrow = 5)) # Period 1
                                   #  matrix(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 1 , NA, NA, NA , NA , NA, NA, NA, NA , NA , NA , NA, NA), nrow = 5),
                                   #  matrix(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 1 , NA, NA, NA , NA , NA, NA, NA, NA , NA , NA , NA, NA), nrow = 5),
                                   #  matrix(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 1 , NA, NA, NA , NA , NA, NA, NA, NA , NA , NA , NA, NA), nrow = 5),
                                   #  matrix(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 1 , NA, NA, NA , NA , NA, NA, NA, NA , NA , NA , NA, NA), nrow = 5))

SVAR_model_2nd_scheme = RWZreject_modified(Y = SVAR_Fig_2_data, nlags = 2, draws = 10000, subdraws = 10000, nkeep = 1000, zero = FALSE, constrained = Sign_restrictions_2nd_scheme, constant =  TRUE, steps =  10)

IRFs_2nd_scheme = SVAR_model_2nd_scheme$IRFS

IRFs = apply(IRFs_2nd_scheme, MARGIN = c(2,3,4) , FUN = median)
IRFs_lower_bound_start = apply(IRFs_2nd_scheme, MARGIN = c(2,3,4) , FUN = quantile, probs =  Confidence_bands[1] / 100)
IRFs_upper_bound_start = apply(IRFs_2nd_scheme, MARGIN = c(2,3,4) , FUN = quantile, probs =  Confidence_bands[2] / 100)

IRFs_median = Data_retrieval(IRFs, nb_var = 5, nb_shocks = 5, horizon = 10)[3,]
IRFs_lower_bound = Data_retrieval(IRFs_lower_bound_start, nb_var = 5, nb_shocks = 5, horizon = 10)[3,]
IRFs_upper_bound = Data_retrieval(IRFs_upper_bound_start, nb_var = 5, nb_shocks = 5, horizon = 10)[3,]

Plot_IRFs_2nd_scheme = plot_figure_new(IRFs_median, IRFs_lower_bound, IRFs_upper_bound, horizon = 10, nb_var = 5 , shock = "Asset Purchase shock" , var_names = var_names   )

# Third identification scheme

Sign_restrictions_3rd_scheme = list(matrix(c(-1, 1, NA, NA, NA, 1, 1, NA, NA, NA, 0, 0, 1 , 1, NA, NA , NA , NA, NA, NA, NA , NA , 1 , -1, NA), nrow = 5), # Period 0
                                    matrix(c(-1, 1, NA, NA, NA, 1, 1, NA, NA, NA, 0, 0, 1 , 1, NA, NA , NA , NA, NA, NA, NA , NA , 1 , -1, NA), nrow = 5)) # Period 1
#  matrix(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 1 , NA, NA, NA , NA , NA, NA, NA, NA , NA , NA , NA, NA), nrow = 5),
#  matrix(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 1 , NA, NA, NA , NA , NA, NA, NA, NA , NA , NA , NA, NA), nrow = 5),
#  matrix(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 1 , NA, NA, NA , NA , NA, NA, NA, NA , NA , NA , NA, NA), nrow = 5),
#  matrix(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 1 , NA, NA, NA , NA , NA, NA, NA, NA , NA , NA , NA, NA), nrow = 5))

debug(RWZreject_modified)

SVAR_model_3rd_scheme = RWZreject_modified(Y = SVAR_Fig_2_data, nlags = 2, draws = 10000, subdraws = 10000, nkeep = 1000, zero = TRUE, constrained = Sign_restrictions_3rd_scheme, constant =  TRUE, steps =  10)

IRFs_3rd_scheme = SVAR_model_3rd_scheme$IRFS

IRFs = apply(IRFs_3rd_scheme, MARGIN = c(2,3,4) , FUN = median)
IRFs_lower_bound_start = apply(IRFs_3rd_scheme, MARGIN = c(2,3,4) , FUN = quantile, probs =  Confidence_bands[1] / 100)
IRFs_upper_bound_start = apply(IRFs_3rd_scheme, MARGIN = c(2,3,4) , FUN = quantile, probs =  Confidence_bands[2] / 100)

IRFs_median = Data_retrieval(IRFs, nb_var = 5, nb_shocks = 5, horizon = 10)[3,]
IRFs_lower_bound = Data_retrieval(IRFs_lower_bound_start, nb_var = 5, nb_shocks = 5, horizon = 10)[3,]
IRFs_upper_bound = Data_retrieval(IRFs_upper_bound_start, nb_var = 5, nb_shocks = 5, horizon = 10)[3,]

Plot_IRFs_3rd_scheme = plot_figure_new(IRFs_median, IRFs_lower_bound, IRFs_upper_bound, horizon = 10, nb_var = 5 , shock = "Asset Purchase shock" , var_names = var_names   )

# Fourth identification scheme

Sign_restrictions_4th_scheme = list(matrix(c(-1, 1, NA, NA, NA, 1, 1, NA, NA, NA, NA, NA, 1 , NA, NA, NA , NA , NA, NA, NA, NA , NA , NA , NA, NA), nrow = 5), # Period 0
                                    matrix(c(-1, 1, NA, NA, NA, 1, 1, NA, NA, NA, NA, NA, 1 , NA, NA, NA , NA , NA, NA, NA, NA , NA , NA , NA, NA), nrow = 5)) # Period 1
#  matrix(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 1 , NA, NA, NA , NA , NA, NA, NA, NA , NA , NA , NA, NA), nrow = 5),
#  matrix(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 1 , NA, NA, NA , NA , NA, NA, NA, NA , NA , NA , NA, NA), nrow = 5),
#  matrix(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 1 , NA, NA, NA , NA , NA, NA, NA, NA , NA , NA , NA, NA), nrow = 5),
#  matrix(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 1 , NA, NA, NA , NA , NA, NA, NA, NA , NA , NA , NA, NA), nrow = 5))

FEVD_check = list(matrix(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0, 0, 1 , NA, NA, NA , NA , NA, NA, NA, NA , NA , NA , NA, NA), nrow = 5),
                  matrix(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0, 0, 1 , NA, NA, NA , NA , NA, NA, NA, NA , NA , NA , NA, NA), nrow = 5),
                  matrix(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 0, 0, 1 , NA, NA, NA , NA , NA, NA, NA, NA , NA , NA , NA, NA), nrow = 5))

SVAR_model_4th_scheme = RWZreject_modified(Y = SVAR_Fig_2_data, nlags = 2, draws = 10000, subdraws = 10000, nkeep = 1000, zero = TRUE, constrained = Sign_restrictions_3rd_scheme, FEVD_check = FEVD_check, constant =  TRUE, steps =  10)

