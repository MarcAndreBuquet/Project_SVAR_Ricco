
# _______________________________________________________________________   
#' @Title : Code SVAR Weale Wieladek
#' @author: Marc-André Buquet
#' @date: 26/11/2023
#  Mail: marc-andre.buquet@ensae.fr  
# _______________________________________________________________________       


## clean environment

rm(list = ls())
setwd("C:/Users/mabuq/Documents/M2_ENSAE/VAR-LP/Projet")

## Library

library(remotes)
library(minqa)
library(mvnfast)
library(HI) #install directly from the R archives as not available anymore from CRAN direct installation 

remotes::install_github("https://github.com/cran/VARsignR/tree/master")
remotes::install_github("roootra/ZerosignR")

library(VARsignR)
library(ZerosignR)
library(BVAR)



## Source functions from external repositories in different files

source("C:/Users/mabuq/Documents/M2_ENSAE/VAR-LP/Projet/function_package_SVAR_sign.R")
source("C:/Users/mabuq/Documents/M2_ENSAE/VAR-LP/Projet/code_R/Fonctions_packages.R")

## Retrieve custom functions to make the code run

source("C:/Users/mabuq/Documents/M2_ENSAE/VAR-LP/Projet/All_functions_VARsignR_modified.R")

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

plot(x = 1: length(VAR_data_1st$CPI_US_log), y =  VAR_data_1st$CPI_US_log) # possibly trend
plot(x = 1: length(VAR_data_1st$CPI_US_log), y =  VAR_data_1st$Real_GDP_US_log) # possibly trend
plot(x = 1: length(VAR_data_1st$CPI_US_log), y =  VAR_data_1st$Asset_purchases_US)
plot(x = 1: length(VAR_data_1st$CPI_US_log), y =  VAR_data_1st$Yield_US_10Y)
plot(x = 1: length(VAR_data_1st$CPI_US_log), y =  VAR_data_1st$Real_Equity_Prices_US_log)

## Compute the RF VAR first

VARselect(VAR_data_1st, lag.max = 12, type = c("both"))
          

RF_VAR = vars::VAR(VAR_data_1st, type = "trend" )

RF_VAR

# Setup the Cholesky matrix with ordering as done in their paper

Cholesky_matrix = diag(1,5)
Cholesky_matrix[lower.tri(Cholesky_matrix)] <- NA

SVAR_Cholesky = SVAR(RF_VAR, Bmat = Cholesky_matrix )

IRFs_Cholesky_US = vars::irf(SVAR_Cholesky, ci = 0.68) # IRFs computed for Cholesky for US (do the same for UK)


### Compute SVAR for the second scheme ------


## Some tests

# Set impulse responses to a horizon of 20 time periods and enable FEVD
# (Identification is performed via Cholesky decomposition)
bv_irf(horizon = 20, fevd = TRUE)

# Set up structural impulse responses using sign restrictions
signs <- matrix(c(NA, NA, NA, NA, NA, -1, -1, 1, 1), nrow = 3)
bv_irf(sign_restr = signs)

# Set up structural impulse responses using zero and sign restrictions
zero_signs <- matrix(c(NA, 0, NA, NA, NA, NA, -1, 1, 1), nrow = 3)
test_2 = bv_irf(sign_restr = zero_signs)

# Prepare to estimate unidentified impulse responses
bv_irf(identification = FALSE)


# Access a subset of the fred_qd dataset
data <- fred_qd[, c("CPIAUCSL", "UNRATE", "FEDFUNDS")]
# Transform it to be stationary
data <- fred_transform(data, codes = c(5, 5, 1), lag = 4)

# Estimate a BVAR using one lag, default settings and very few draws
x <- bvar(data, lags = 1, n_draw = 1000L, n_burn = 200L, verbose = FALSE, irf = bv_irf(sign_restr = zero_signs) )

# Calculate and store forecasts and impulse responses
predict(x) <- predict(x, horizon = 8)
irf(x) <- irf(x, horizon = 8, fevd = FALSE)

## Not run: 
# Check convergence of the hyperparameters with a trace and density plot
plot(x)
# Plot forecasts and impulse responses
plot(predict(x))
plot(irf(x))
# Check coefficient values and variance-covariance matrix
summary(x)




### Compute SVAR for the third scheme -------


## Some tests using built-in data and functions

set.seed(12345)
data(uhligdata)

# variable labels for plots
vl <- c("GDP","GDP Deflator","Comm.Pr.Index","Fed Funds Rate",
        "NB Reserves", "Total Reserves")

# sign restrictions
# shock of interest enters first.
# you MUST provide a restriction for the shock of interest
# restriction variable 4 is >0
# 2nd, 3rd, and 5th variable are <0.
# 1st and 6th variable are unrestricted

constr <- c(+4,-3,-2,-5)

# estimates the model
mode$l3 <- VARsignR::rwz.reject(Y=uhligdata, nlags=12, draws=200, subdraws=200, nkeep=1000,
                     KMIN=1, KMAX=6, constrained=constr, constant=FALSE, steps=60)

# get posterior draws
irfs0 <- model3$IRFS

# plot impulse response functions

vl <- c("GDP","GDP Deflator","Comm.Pr.Index","Fed Funds Rate",
        "NB Reserves", "Total Reserves")


irfplot(irfdraws=irfs0, type="mean", labels=vl, save=FALSE, bands=c(0.16, 0.84),
        grid=TRUE, bw=FALSE)



## Test sign restrictions to implement


test_sign_restrictions = list(matrix(c(NA, 0, NA, NA, NA, NA, -1, 1, 1), nrow = 3), matrix(c(NA, 0, NA, 1, NA, NA, -1, 0, 1), nrow = 3))

#### Test Fevd

rdm_col = as.matrix(runif(5)) 
rdm_col[2] 
test_true = if (rdm_col[3] == max(rdm_col)) {T} else {F}
test_true




######## Test avec mes propres fonctions 

zero = TRUE
test_sign_restrictions = list(matrix(c(NA, 0, NA, NA, NA, NA, -1, 1, 1), nrow = 3), matrix(c(NA, 0, NA, 1, NA, NA, -1, 0, 1), nrow = 3))
constrained = test_sign_restrictions
sign_restr = test_sign_restrictions
a <- matrix(rnorm(nvar^2, mean = 0, sd = 1), nvar, nvar)
qr_object <- qr(matrix(rnorm(nvar^2, 0, 1), nvar, nvar))
Q <- qr.Q(qr_object)
shock = t(fevd0[3,,] %*% Q^2)
test_sum = sum(shock[,6])
test_sum

View(shock)
shock_vec <- as.vector(shock)
View(shock_vec)
shock_vec[which(shock_vec > 0)] <- 1
View(sign_vec)

FEVD_check = list(matrix(c(NA, NA, NA, NA, NA , NA, NA, NA, NA, NA,0,0,1,NA,NA,NA, NA, NA, NA, NA , NA, NA, NA, NA, NA), nrow = 5))
View(fevd_vec)
fevd_sign_vec = as.vector(FEVD_check[[1]])
View(fevd_sign_vec)
restricted <- which(!is.na(fevd_sign_vec))
View(fevd_vec_restricted)
