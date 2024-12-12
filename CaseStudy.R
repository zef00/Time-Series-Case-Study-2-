rm(list=ls())

#libraries 
library(ggplot2) # to make nice looking graphs
library(tidyverse) # keeping things tidy
library(vars) # for the VARS 
library(svars) # for the SVARS
library(patchwork) # for nice graphs (used for irf plots)
library(expm) #multiplying matrices more easily
library(tseries) #residual diagnostics
library(FinTS) #residual diagnostics
library(lmtest) #Granger causality
library(devtools) # used to extract the exctract_irf function
source_url("https://raw.githubusercontent.com/anguyen1210/var-tools/master/R/extract_varirf.R") 

# Data set used
FRED_url <- url("https://files.stlouisfed.org/files/htdocs/fred-md/monthly/current.csv")
FRED_MD <- read.csv(FRED_url)

data <- FRED_MD
data <- data[-1, ]
data$sasdate <- as.Date(data$sasdate, format = "%m/%d/%Y")

data <- data %>% 
  dplyr::select(sasdate, INDPRO, CPIAUCSL, FEDFUNDS) %>%
  dplyr::filter(sasdate <= "2019-12-01")

# Plotting the data before transformation
ggplot(data, aes(x = sasdate, y = INDPRO)) +
  geom_line(color = "red") +
  theme_minimal() +
  labs(
    title = "INDPRO (Original)",
    x = "Date",
    y = "INDPRO"
  )

ggplot(data, aes(x = sasdate, y = CPIAUCSL)) +
  geom_line(color = "blue") +
  theme_minimal() +
  labs(
    title = "CPIAUCSL (Original)",
    x = "Date",
    y = "CPIAUCSL"
  )

ggplot(data, aes(x = sasdate, y = FEDFUNDS)) +
  geom_line(color = "green") +
  theme_minimal() +
  labs(
    title = "FEDFUNDS (Original)",
    x = "Date",
    y = "FEDFUNDS"
  )

# Transforming the data
    # INPDPRO: code 5 = delta(log(xt)); FEDFUNDS: code 2 = delta(xt); CPIAUCSL: code 6 = delta^2(log(xt))
data <- data %>%
  mutate(log_diff_indpro = c(NA, diff(log(data$INDPRO))),
         diff_fedfunds = c(NA, diff(data$FEDFUNDS)),
         log_diff2_cpi = c(NA, NA, diff(log(data$CPIAUCSL), differences = 2))) 

ggplot(data, aes(x = sasdate)) +
  geom_line(aes(y = log_diff_indpro), color = "red") +
  theme_minimal() +
  labs(
    title = "INDPRO (Transformed)",
    x = "Date",
    y = "Log difference of INDPRO"
  )

ggplot(data, aes(x = sasdate)) +
  geom_line(aes(y = log_diff2_cpi), color = "blue") +
  theme_minimal() +
  labs(
    title = "CPIAUCSL (Transformed)",
    x = "Date",
    y = "Second log difference of CPIAUCSL"
  )

ggplot(data, aes(x = sasdate)) +
  geom_line(aes(y = diff_fedfunds), color = "green") +
  theme_minimal() +
  labs(
    title = "FEDFUNDS (Transformed)",
    x = "Date",
    y = "First difference of FEDFUNDS"
  )


# VAR Model 
var_data <- dplyr::select(data, log_diff_indpro, diff_fedfunds, log_diff2_cpi)
var_data <- var_data[-(1:2), ]

lag_selection <- VARselect(var_data, lag.max = 10)
lag_selection #Lag selection gave us the following lag orders: 3 (SC), 9 (HQ), and 10 (AIC/FPE).

#Fitting var models for each lag order.
#LAG ORDER = 3

#Fitting the model
var_model3 <- VAR(var_data, p = 3, type = "const")
summary(var_model3)

#Residual diagnostics entire VAR
serial.test(var_model3, lags.pt = 16, type = "PT.asymptotic")
normality.test(var_model3)
arch.test(var_model3, lags.multi = 5)

#Residual diagnostics per equation
residuals_indpro <- residuals(var_model3)[, "log_diff_indpro"]
Box.test(residuals_indpro, lag = 16, type = "Ljung-Box", fitdf = 3)
jarque.bera.test(residuals_indpro)
ArchTest(residuals_indpro, lags = 5)

residuals_cpiaucsl <- residuals(var_model3)[, "log_diff2_cpi"]
Box.test(residuals_cpiaucsl, lag = 16, type = "Ljung-Box", fitdf = 3)
jarque.bera.test(residuals_cpiaucsl)
ArchTest(residuals_cpiaucsl, lags = 5)

residuals_fedfunds <- residuals(var_model3)[, "diff_fedfunds"]
Box.test(residuals_fedfunds, lag = 16, type = "Ljung-Box", fitdf = 3)
jarque.bera.test(residuals_fedfunds)
ArchTest(residuals_fedfunds, lags = 5)

#acf and pacf of residuals
combined_residuals3 <- as.vector(residuals(var_model3))
acf(combined_residuals3, main = "ACF of Combined Residuals for p=3")
pacf(combined_residuals3, main = "PACF of Combined Residuals for p=3")
acf(residuals(var_model3)[, "log_diff_indpro"], main = "ACF of Residuals: INDPRO for p=3")
pacf(residuals(var_model3)[, "log_diff_indpro"], main = "PACF of Residuals: INDPRO for p=3")
acf(residuals(var_model3)[, "log_diff2_cpi"], main = "ACF of Residuals: CPIAUCSL for p=3")
pacf(residuals(var_model3)[, "log_diff2_cpi"], main = "PACF of Residuals: CPIAUCSL for p=3")
acf(residuals(var_model3)[, "diff_fedfunds"], main = "ACF of Residuals: FEDFUNDS for p=3")
pacf(residuals(var_model3)[, "diff_fedfunds"], main = "PACF of Residuals: FEDFUNDS for p=3")

#LAG ORDER = 9
var_model9 <- VAR(var_data, p = 9, type = "const")
summary(var_model9)

#Residual diagnostics entire VAR
serial.test(var_model9, lags.pt = 16, type = "PT.asymptotic")
normality.test(var_model9)
arch.test(var_model9, lags.multi = 5)

#Residual diagnostics per equation
residuals_indpro <- residuals(var_model9)[, "log_diff_indpro"]
Box.test(residuals_indpro, lag = 16, type = "Ljung-Box", fitdf = 9)
jarque.bera.test(residuals_indpro)
ArchTest(residuals_indpro, lags = 5)

residuals_cpiaucsl <- residuals(var_model9)[, "log_diff2_cpi"]
Box.test(residuals_cpiaucsl, lag = 16, type = "Ljung-Box", fitdf = 9)
jarque.bera.test(residuals_cpiaucsl)
ArchTest(residuals_cpiaucsl, lags = 5)

residuals_fedfunds <- residuals(var_model9)[, "diff_fedfunds"]
Box.test(residuals_fedfunds, lag = 16, type = "Ljung-Box", fitdf = 9)
jarque.bera.test(residuals_fedfunds)
ArchTest(residuals_fedfunds, lags = 5)

#acf of residuals
combined_residuals9 <- as.vector(residuals(var_model9))
acf(combined_residuals9, main = "ACF of Combined Residuals for p=9")
pacf(combined_residuals9, main = "PACF of Combined Residuals for p=9")
acf(residuals(var_model9)[, "log_diff_indpro"], main = "ACF of Residuals: INDPRO for p=9")
pacf(residuals(var_model9)[, "log_diff_indpro"], main = "PACF of Residuals: INDPRO for p=9")
acf(residuals(var_model9)[, "log_diff2_cpi"], main = "ACF of Residuals: CPIAUCSL for p=9")
pacf(residuals(var_model9)[, "log_diff2_cpi"], main = "PACF of Residuals: CPIAUCSL for p=9")
acf(residuals(var_model9)[, "diff_fedfunds"], main = "ACF of Residuals: FEDFUNDS for p=9")
pacf(residuals(var_model9)[, "diff_fedfunds"], main = "PACF of Residuals: FEDFUNDS for p=9")

#LAG ORDER = 10

#Fitting the model
var_model10 <- VAR(var_data, p = 10, type = "const")
summary(var_model10)

#Residual diagnostics entire VAR
serial.test(var_model10, lags.pt = 16, type = "PT.asymptotic")
normality.test(var_model10)
arch.test(var_model10, lags.multi = 5)

#Residual diagnostics per equation
residuals_indpro <- residuals(var_model10)[, "log_diff_indpro"]
Box.test(residuals_indpro, lag = 16, type = "Ljung-Box", fitdf = 10)
jarque.bera.test(residuals_indpro)
ArchTest(residuals_indpro, lags = 5)

residuals_cpiaucsl <- residuals(var_model10)[, "log_diff2_cpi"]
Box.test(residuals_cpiaucsl, lag = 16, type = "Ljung-Box", fitdf = 10)
jarque.bera.test(residuals_cpiaucsl)
ArchTest(residuals_cpiaucsl, lags = 5)

residuals_fedfunds <- residuals(var_model10)[, "diff_fedfunds"]
Box.test(residuals_fedfunds, lag = 16, type = "Ljung-Box", fitdf = 10)
jarque.bera.test(residuals_fedfunds)
ArchTest(residuals_fedfunds, lags = 5)

#acf of residuals
combined_residuals10 <- as.vector(residuals(var_model10))
acf(combined_residuals10, main = "ACF of Combined Residuals for p=10")
pacf(combined_residuals10, main = "PACF of Combined Residuals for p=10")
acf(residuals(var_model10)[, "log_diff_indpro"], main = "ACF of Residuals: INDPRO for p=10")
pacf(residuals(var_model10)[, "log_diff_indpro"], main = "PACF of Residuals: INDPRO for p=10")
acf(residuals(var_model10)[, "log_diff2_cpi"], main = "ACF of Residuals: CPIAUCSL for p=10")
pacf(residuals(var_model10)[, "log_diff2_cpi"], main = "PACF of Residuals: CPIAUCSL for p=10")
acf(residuals(var_model10)[, "diff_fedfunds"], main = "ACF of Residuals: FEDFUNDS for p=10")
pacf(residuals(var_model10)[, "diff_fedfunds"], main = "PACF of Residuals: FEDFUNDS for p=10")

#Mean square error of each model
overall_mse3 <- mean(combined_residuals3^2)
print(overall_mse3)
overall_mse9 <- mean(combined_residuals9^2)
print(overall_mse9)
overall_mse10 <- mean(combined_residuals10^2)
print(overall_mse10)
#Mean square error of each equation per model
var_model <- VAR(var_data, p=3, type = "const")
residuals <- as.matrix(residuals(var_model))
mse_per_equation <- colMeans(residuals^2)
print(mse_per_equation)
var_model <- VAR(var_data, p=9, type = "const")
residuals <- as.matrix(residuals(var_model))
mse_per_equation <- colMeans(residuals^2)
print(mse_per_equation)
var_model <- VAR(var_data, p=10, type = "const")
residuals <- as.matrix(residuals(var_model))
mse_per_equation <- colMeans(residuals^2)
print(mse_per_equation)

#Selected model: VAR(3)
optimal_lag <- lag_selection$selection['SC(n)']
var_model <- VAR(var_data, p=optimal_lag, type = "const", ic="SC")
summary(var_model)

#Coefficients
coefficients <- coef(var_model)
print(coefficients)

#Stability of model
stability <- stability(var_model)
plot(stability) 
stability_roots <- roots(var_model)  
print(stability)


# Granger causality

#INDPRO
#Does INDPRO Granger cause CPIAUCSl and FEDFUNDS?
#Joint test:
granger_indpro <- causality(var_model, cause="log_diff_indpro")
print(granger_indpro)
#Pairwise test 1: Does INDPRO Granger cause CPIAUCSL?
grangertest(log_diff2_cpi ~ log_diff_indpro, order = 3, data = var_data)
#Pairwise test 2:Does INDPRO Granger cause FEDFUNDS?
grangertest(diff_fedfunds ~ log_diff_indpro, order = 3, data = var_data)

#CPIAUCSL
#Does CPIAUCSL Granger cause INDPRO and FEDFUNDS?
#Joint test:
granger_cpiaucsl <- causality(var_model, cause="log_diff2_cpi")
print(granger_cpiaucsl)
#Pairwise test 1: Does CPIAUCSL Granger cause INDPRO?
grangertest(log_diff_indpro ~ log_diff2_cpi, order = 3, data = var_data)
#Pairwise test 2:Does CPIAUCSL Granger cause FEDFUNDS?
grangertest(diff_fedfunds ~ log_diff2_cpi, order = 3, data = var_data)

#FEDFUNDS
#Does FEDFUNDS Granger cause INDPRO and CPIAUCSL?
#Joint test:
granger_fedfunds <- causality(var_model, cause="diff_fedfunds")
print(granger_fedfunds)
#Pairwise test 1: Does FEDFUNDS Granger cause INDPRO?
grangertest(log_diff_indpro ~ diff_fedfunds, order = 3, data = var_data)
#Pairwise test 2:Does FEDFUNDS Granger cause CPIAUCSL?
grangertest(log_diff2_cpi ~ diff_fedfunds, order = 3, data = var_data)

granger_test_cpi <- causality(var_model, cause = "log_diff2_cpi")
print(granger_test_cpi)


# IRF reduced form

# For ez finding of the IRFs (these are the same as below but there they are seperated per variable)
IRF_log_diff_indpro <- irf(var_model,impulse = "log_diff_indpro", n.ahead = 24,
                           ortho = FALSE, cumulative = FALSE, boot = TRUE, ci = 0.95,
                           runs = 100)
irf_values_ld_ip <- extract_varirf(IRF_log_diff_indpro) # extract all the values

IRF_diff_fedfunds <- irf(var_model,impulse = "diff_fedfunds", n.ahead = 24,
                         ortho = FALSE, cumulative = FALSE, boot = TRUE, ci = 0.95,
                         runs = 100)
irf_values_d_ff <- extract_varirf(IRF_diff_fedfunds)

IRF_log_diff2_cpi <- irf(var_model,impulse = "log_diff2_cpi", n.ahead = 24,
                         ortho = FALSE, cumulative = FALSE, boot = TRUE, ci = 0.95,
                         runs = 100)
irf_values_l2d_cpi <- extract_varirf(IRF_log_diff2_cpi)

# IRF log differences INDPRO 
IRF_log_diff_indpro <- irf(var_model,impulse = "log_diff_indpro", n.ahead = 24,
                           ortho = FALSE, cumulative = FALSE, boot = TRUE, ci = 0.95,
                           runs = 100)
plot(IRF_log_diff_indpro) # can't see the axis well so need to adjust 
irf_values_ld_ip <- extract_varirf(IRF_log_diff_indpro) # extract all the values

# Make the plots of all the IRF's of INDPRO in ggplot 
ld_ip_ld_ip <- irf_values_ld_ip %>% # INDPRO - INDPRO
  ggplot(aes(x=period, y=irf_log_diff_indpro_log_diff_indpro, ymin=lower_log_diff_indpro_log_diff_indpro, 
             ymax=upper_log_diff_indpro_log_diff_indpro)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Impulse response, INDPRO - INDPRO")+
  ylab("log(∆INDPRO)")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

ld_ip_d_ff <- irf_values_ld_ip %>% # INDPRO - FEDFUNDS
  ggplot(aes(x=period, y=irf_log_diff_indpro_diff_fedfunds, ymin=lower_log_diff_indpro_diff_fedfunds, 
             ymax=upper_log_diff_indpro_diff_fedfunds)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Impulse response, INDPRO - FEDFUNDS")+
  ylab("∆FEDFUNDS")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

ld_ip_l2d_cpi <- irf_values_ld_ip %>% # INDPRO - CPIAUCSL
  ggplot(aes(x=period, y=irf_log_diff_indpro_log_diff2_cpi, ymin=lower_log_diff_indpro_log_diff2_cpi, 
             ymax=upper_log_diff_indpro_log_diff2_cpi)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Impulse response, INDPRO - CPI")+
  ylab("log(∆^2(CPI))")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

ld_ip_irf <- ld_ip_ld_ip / ld_ip_d_ff/ ld_ip_l2d_cpi # Combine the IRF plots in a nice single plot
ld_ip_irf # The final reduced-form IRF of log differenced INDPRO 

# IRF for differenced FEDFUNDS 
IRF_diff_fedfunds <- irf(var_model,impulse = "diff_fedfunds", n.ahead = 24,
                         ortho = FALSE, cumulative = FALSE, boot = TRUE, ci = 0.95,
                         runs = 100)
irf_values_d_ff <- extract_varirf(IRF_diff_fedfunds)

# Make the plots of all the IRF's of FEDFUNDS in ggplot 
d_ff_ld_ip <- irf_values_d_ff %>% 
  ggplot(aes(x=period, y=irf_diff_fedfunds_log_diff_indpro, ymin=lower_diff_fedfunds_log_diff_indpro, 
             ymax=upper_diff_fedfunds_log_diff_indpro)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Impulse response, FEDFUNDS - INDPRO")+
  ylab("log(∆INDPRO)")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

d_ff_d_ff <- irf_values_d_ff %>% 
  ggplot(aes(x=period, y=irf_diff_fedfunds_diff_fedfunds, ymin=lower_diff_fedfunds_diff_fedfunds, 
             ymax=upper_diff_fedfunds_diff_fedfunds)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Impulse response, FEDFUNDS - FEDFUNDS")+
  ylab("∆FEDFUNDS")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

d_ff_l2d_cpi <- irf_values_d_ff %>% 
  ggplot(aes(x=period, y=irf_diff_fedfunds_log_diff2_cpi, ymin=lower_diff_fedfunds_log_diff2_cpi, 
             ymax=upper_diff_fedfunds_log_diff2_cpi)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Impulse response, FEDFUNDS - CPI")+
  ylab("log(∆^2(CPI))")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

d_ff_irf <- d_ff_ld_ip / d_ff_d_ff / d_ff_l2d_cpi 
d_ff_irf # The final reduced-form IRF of differenced FEDFUNDS 

# IRF for the log twice differenced CPIAUSL 
IRF_log_diff2_cpi <- irf(var_model,impulse = "log_diff2_cpi", n.ahead = 24,
                         ortho = FALSE, cumulative = FALSE, boot = TRUE, ci = 0.95,
                         runs = 100)
irf_values_l2d_cpi <- extract_varirf(IRF_log_diff2_cpi)

l2d_cpi_ld_ip <- irf_values_l2d_cpi %>% # CPIAUCSL - INDPRO
  ggplot(aes(x=period, y=irf_log_diff2_cpi_log_diff_indpro, ymin=lower_log_diff2_cpi_log_diff_indpro, 
             ymax=upper_log_diff2_cpi_log_diff_indpro)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Impulse response, CPI - INDPRO")+
  ylab("log(∆INDPRO)")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

l2d_cpi_d_ff <- irf_values_l2d_cpi %>% # CPIAUCSL - FEDFUNDS
  ggplot(aes(x=period, y=irf_log_diff2_cpi_diff_fedfunds, ymin=lower_log_diff2_cpi_diff_fedfunds, 
             ymax=upper_log_diff2_cpi_diff_fedfunds)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Impulse response, CPI - FEDFUNDS")+
  ylab("∆FEDFUNDS")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

l2d_cpi_l2d_cpi <- irf_values_l2d_cpi %>% # CPIAUCSL - CPIAUCSL
  ggplot(aes(x=period, y=irf_log_diff2_cpi_log_diff2_cpi, ymin=lower_log_diff2_cpi_log_diff2_cpi, 
             ymax=upper_log_diff2_cpi_log_diff2_cpi)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Impulse response, CPI - CPI")+
  ylab("log(∆^2(CPI))")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

l2d_cpi_irf <- l2d_cpi_ld_ip / l2d_cpi_d_ff / l2d_cpi_l2d_cpi
l2d_cpi_irf # The final reduced-form IRF of log twice differenced CPIAUCSL


# IRF of the original variables

# Tansforming starting from 1
transformation <- function(x){ # x is the irf values dataframe 
  irf <- data_frame(rep(1:nrow(x))) 
  for (i in 1:ncol(x)){
    column <- x[,i]
    if (i == 1){
      column = period
    } 
    else if (i %in% c(2, 5, 8)) { # Transform for the log differenced INDPRO 
      indpro <- exp(cumsum(column))
      irf <- cbind(irf, indpro)
    } 
    else if (i %in% c(3, 6, 9)) { # Transform for the differenced FEDFUNDS
      fedfunds <- cumsum(column)
      irf <- cbind(irf, fedfunds)
    }
    else if (i %in% c(4, 7, 10)) { # Transform for the log twice differenced CPIAUCSL
      cpi <- exp(cumsum(c(0,cumsum(column))))
      cpi <- cpi[-26]
      irf <- cbind(irf, cpi)
    }
  }
  colnames(irf) <- c("Period","indpro", "fedfunds", "cpi",    # Putting the right names on the column
                     "lower_indpro", "lower_ff", "lower_cpi",
                     "upper_indpro", "upper_ff", "upper_cpi")
  return(irf) # Return the IRF with the pre-transformed values 
}

# The back transformed IRF's (level equation)
irf_indpro <- transformation(irf_values_ld_ip)
irf_fedfunds <- transformation(irf_values_d_ff)
irf_cpi <- transformation(irf_values_l2d_cpi)

# Plotting the back transformed IRF's (level equation)
# INDPRO
indpro_indpro <- irf_indpro %>% # INDPRO - INDPRO
  ggplot(aes(x=Period, y=indpro, ymin=lower_indpro, ymax=upper_indpro)) +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Impulse response, INDPRO - INDPRO")+
  ylab("INDPRO")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

indpro_fedfunds <- irf_indpro %>% # INDPRO - FEDFUNDS
  ggplot(aes(x=Period, y=fedfunds, ymin=lower_ff, ymax=upper_ff)) +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Impulse response, INDPRO - FEDFUNDS")+
  ylab("FEDFUNDS")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

indpro_cpi <- irf_indpro %>% # INDPRO - CPI
  ggplot(aes(x=Period, y=cpi, ymin=lower_cpi, ymax=upper_cpi)) +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Impulse response, INDPRO - CPI")+
  ylab("CPI")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

irf_plot_indpro <- indpro_indpro / indpro_fedfunds / indpro_cpi

# FEDFUNDS
fedfunds_indpro <- irf_fedfunds %>% # INDPRO - INDPRO
  ggplot(aes(x=Period, y=indpro, ymin=lower_indpro, ymax=upper_indpro)) +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Impulse response, FEDFUNDS - INDPRO")+
  ylab("INDPRO")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

fedfunds_fedfunds <- irf_fedfunds %>% # FEDFUNDS - FEDFUNDS
  ggplot(aes(x=Period, y=fedfunds, ymin=lower_ff, ymax=upper_ff)) +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Impulse response, FEDFUNDS - FEDFUNDS")+
  ylab("FEDFUNDS")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

fedfunds_cpi <- irf_fedfunds %>% # FEDFUNDS - CPI
  ggplot(aes(x=Period, y=cpi, ymin=lower_cpi, ymax=upper_cpi)) +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Impulse response, FEDFUNDS - CPI")+
  ylab("CPI")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

irf_plot_fedfunds <- fedfunds_indpro / fedfunds_fedfunds / fedfunds_cpi

# CPIAUCSL
cpi_indpro <- irf_cpi %>% # CPI - INDPRO
  ggplot(aes(x=Period, y=indpro)) +
  geom_ribbon(aes(ymin=lower_indpro, ymax=upper_indpro),fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Impulse response, CPI - INDPRO")+
  ylab("INDPRO")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

cpi_fedfunds <- irf_cpi %>% # CPI - FEDFUNDS
  ggplot(aes(x=Period, y=fedfunds, ymin=lower_ff, ymax=upper_ff)) +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Impulse response, CPI - FEDFUNDS")+
  ylab("FEDFUNDS")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

cpi_cpi <- irf_cpi %>% # FEDFUNDS - CPI
  ggplot(aes(x=Period, y=cpi)) +
  geom_ribbon(aes(ymin=lower_cpi, ymax=upper_cpi),fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Impulse response, CPI - CPI")+
  ylab("CPI")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

irf_plot_cpi <- cpi_indpro / cpi_fedfunds / cpi_cpi


# IRF of Recursive form VAR

# ez access to the irfs (same as with the reduced-form var irfs)
IRF_log_diff_indpro <- irf(var_model,impulse = "log_diff_indpro", n.ahead = 24,
                           ortho = TRUE, cumulative = FALSE, boot = TRUE, ci = 0.95,
                           runs = 100) # put cumulative = T if you want levels
irf_values_ld_ip <- extract_varirf(IRF_log_diff_indpro) # extract all the values

IRF_diff_fedfunds <- irf(var_model,impulse = "diff_fedfunds", n.ahead = 24,
                         ortho = TRUE, cumulative = FALSE, boot = TRUE, ci = 0.95,
                         runs = 100)
irf_values_d_ff <- extract_varirf(IRF_diff_fedfunds)

IRF_log_diff2_cpi <- irf(var_model,impulse = "log_diff2_cpi", n.ahead = 24,
                         ortho = TRUE, cumulative = FALSE, boot = TRUE, ci = 0.95,
                         runs = 100)
irf_values_l2d_cpi <- extract_varirf(IRF_log_diff2_cpi)

# IRF log differences INDPRO
IRF_log_diff_indpro <- irf(var_model,impulse = "log_diff_indpro", n.ahead = 24,
                           ortho = TRUE, cumulative = FALSE, boot = TRUE, ci = 0.95,
                           runs = 100) # put cumulative = T if you want levels
irf_values_ld_ip <- extract_varirf(IRF_log_diff_indpro) # extract all the values

# Make the plots of all the IRF's of INDPRO in ggplot 
ld_ip_ld_ip <- irf_values_ld_ip %>% # INDPRO - INDPRO
  ggplot(aes(x=period, y=irf_log_diff_indpro_log_diff_indpro, ymin=lower_log_diff_indpro_log_diff_indpro, 
             ymax=upper_log_diff_indpro_log_diff_indpro)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Orthogonal impulse response, INDPRO - INDPRO")+
  ylab("log(∆INDPRO)")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

ld_ip_d_ff <- irf_values_ld_ip %>% # INDPRO - FEDFUNDS
  ggplot(aes(x=period, y=irf_log_diff_indpro_diff_fedfunds, ymin=lower_log_diff_indpro_diff_fedfunds, 
             ymax=upper_log_diff_indpro_diff_fedfunds)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Orthogonal impulse response, INDPRO - FEDFUNDS")+
  ylab("∆FEDFUNDS")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

ld_ip_l2d_cpi <- irf_values_ld_ip %>% # INDPRO - CPIAUCSL
  ggplot(aes(x=period, y=irf_log_diff_indpro_log_diff2_cpi, ymin=lower_log_diff_indpro_log_diff2_cpi, 
             ymax=upper_log_diff_indpro_log_diff2_cpi)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Orthogonal impulse response, INDPRO - CPI")+
  ylab("log(∆^2(CPI))")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

ld_ip_irf <- ld_ip_ld_ip / ld_ip_d_ff/ ld_ip_l2d_cpi # Combine the IRF plots in a nice single plot
ld_ip_irf # The final orthogonalized IRF of log differenced INDPRO 

# IRF for differenced FEDFUNDS 
IRF_diff_fedfunds <- irf(var_model,impulse = "diff_fedfunds", n.ahead = 24,
                         ortho = TRUE, cumulative = FALSE, boot = TRUE, ci = 0.95,
                         runs = 100)
irf_values_d_ff <- extract_varirf(IRF_diff_fedfunds)

# Make the plots of all the IRF's of FEDFUNDS in ggplot 
d_ff_ld_ip <- irf_values_d_ff %>% 
  ggplot(aes(x=period, y=irf_diff_fedfunds_log_diff_indpro, ymin=lower_diff_fedfunds_log_diff_indpro, 
             ymax=upper_diff_fedfunds_log_diff_indpro)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Orthogonal impulse response, FEDFUNDS - INDPRO")+
  ylab("log(∆INDPRO)")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

d_ff_d_ff <- irf_values_d_ff %>% 
  ggplot(aes(x=period, y=irf_diff_fedfunds_diff_fedfunds, ymin=lower_diff_fedfunds_diff_fedfunds, 
             ymax=upper_diff_fedfunds_diff_fedfunds)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Orthogonal impulse response, FEDFUNDS - FEDFUNDS")+
  ylab("∆FEDFUNDS")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

d_ff_l2d_cpi <- irf_values_d_ff %>% 
  ggplot(aes(x=period, y=irf_diff_fedfunds_log_diff2_cpi, ymin=lower_diff_fedfunds_log_diff2_cpi, 
             ymax=upper_diff_fedfunds_log_diff2_cpi)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Orthogonal impulse response, FEDFUNDS - CPI")+
  ylab("log(∆^2(CPI))")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

d_ff_irf <- d_ff_ld_ip / d_ff_d_ff / d_ff_l2d_cpi 
d_ff_irf # The final orthogonalized IRF of differenced FEDFUNDS 

# IRF for the log twice differenced CPIAUSL 
IRF_log_diff2_cpi <- irf(var_model,impulse = "log_diff2_cpi", n.ahead = 24,
                         ortho = TRUE, cumulative = FALSE, boot = TRUE, ci = 0.95,
                         runs = 100)
irf_values_l2d_cpi <- extract_varirf(IRF_log_diff2_cpi)

l2d_cpi_ld_ip <- irf_values_l2d_cpi %>% # CPIAUCSL - INDPRO
  ggplot(aes(x=period, y=irf_log_diff2_cpi_log_diff_indpro, ymin=lower_log_diff2_cpi_log_diff_indpro, 
             ymax=upper_log_diff2_cpi_log_diff_indpro)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Orthogonal impulse response, CPI - INDPRO")+
  ylab("log(∆INDPRO)")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

l2d_cpi_d_ff <- irf_values_l2d_cpi %>% # CPIAUCSL - FEDFUNDS
  ggplot(aes(x=period, y=irf_log_diff2_cpi_diff_fedfunds, ymin=lower_log_diff2_cpi_diff_fedfunds, 
             ymax=upper_log_diff2_cpi_diff_fedfunds)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Orthogonal impulse response, CPI - FEDFUNDS")+
  ylab("∆FEDFUNDS")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

l2d_cpi_l2d_cpi <- irf_values_l2d_cpi %>% # CPIAUCSL - CPIAUCSL
  ggplot(aes(x=period, y=irf_log_diff2_cpi_log_diff2_cpi, ymin=lower_log_diff2_cpi_log_diff2_cpi, 
             ymax=upper_log_diff2_cpi_log_diff2_cpi)) +
  geom_hline(yintercept = 0, color="red") +
  geom_ribbon(fill="grey", alpha=.2, color="grey50", linetype="dashed") +
  geom_line() +
  theme_light() +
  ggtitle("Orthogonal impulse response, CPI - CPI")+
  ylab("log(∆^2(CPI))")+
  xlab("Period") +
  theme(plot.title = element_text(size = 11, hjust=0.5),
        axis.title.y = element_text(size=11))

l2d_cpi_irf <- l2d_cpi_ld_ip / l2d_cpi_d_ff / l2d_cpi_l2d_cpi
l2d_cpi_irf # The final orthogonalized IRF of log twice differenced CPIAUCSL


# Structural VAR 
svar_data <- data.frame(data$log_diff_indpro, data$log_diff2_cpi, data$diff_fedfunds)
svar_data <- svar_data[-(1:2), ]
var_model <- VAR(svar_data, p = 3, type = "const")
svar_chol <- id.chol(var_model)
irf_chol <- irf(svar_chol, impulse = "log_diff_indpro", n.ahead = 24)
plot(irf_chol)


# Extra 

# IRF of the original variables 

# Setting up the last value used in the back transformation 
initial_value_indpro <- tail(data$INDPRO, 1)
initial_value_ff <- tail(data$FEDFUNDS, 1)
initial_value_cpi <- tail(data$CPIAUCSL, 1)
initial_diff <- diff(log(data$CPIAUCSL), differences = 1)[1]

# Function that transforms it back into the original values 
transformation <- function(x){ # x is the irf values dataframe 
  irf <- data_frame(x[,1]) 
  for (i in 1:ncol(x)){
    column <- x[,i]
    if (i == 1){
      column = period
    } 
    else if (i %in% c(2, 5, 8)) { # Transform for the log differenced INDPRO 
      indpro <- initial_value_indpro * exp(cumsum(column))
      irf <- cbind(irf, indpro)
    } 
    else if (i %in% c(3, 6, 9)) { # Transform for the differenced FEDFUNDS
      fedfunds <- initial_value_ff + cumsum(column)
      irf <- cbind(irf, fedfunds)
    }
    else if (i %in% c(4, 7, 10)) { # Transform for the log twice differenced CPIAUCSL
      cpi <- initial_value_cpi + 
        initial_value_cpi*exp(cumsum(c(initial_diff, initial_diff+cumsum(column))))
      cpi <- cpi[-26]
      irf <- cbind(irf, cpi)
    }
  }
  colnames(irf) <- c("Period","indpro", "fedfunds", "cpi",    # Putting the right names on the column
                     "lower_indpro", "lower_ff", "lower_cpi",
                     "upper_indpro", "upper_ff", "upper_cpi")
  return(irf) # Return the IRF with the pre-transformed values 
}

# The back transformed IRF's (from initial value)
irf_indpro <- transformation(irf_values_ld_ip)
irf_fedfunds <- transformation(irf_values_d_ff)
irf_cpi <- transformation(irf_values_l2d_cpi)

# The making of the MA representation of the VAR 
coefficients <- as.data.frame(Bcoef(var_model)) # extract the coefficents
coefficients <- coefficients[,-10] # remove the constant 
coefficients <- as.matrix(coefficients) # put it in a matrix for ease
K <- var_model$K # variables of interest
p <- var_model$p # lags 

I3 <- diag(1,K) # Identity (3x3)-matrix 
nul <-matrix(c(rep(0,3), # 0 (3x3)-matrix
               rep(0,3),
               rep(0,3)), nrow = 3, ncol = 3, byrow = T)

A <- rbind(coefficients, # companion matrix
           cbind(I3, nul, nul),
           cbind(nul, I3, nul))

J <- matrix(c(1,rep(0,8), # Matrix needed to put it into MA form 
              0, 1, rep(0,7),
              0,0,1, rep(0,6)), nrow = K, ncol = p*K, byrow = T)
J_t <- t(J)

psi = J%*%(A)%*%J_t # The MA form coefficient matrix 

# tried to make cumulated effect but i dont think it works 
ir = 0 
for (i in 1:1000){
  ir = as.data.frame(ir + J%*%(A%^%i)%*%J_t)
}

# Epsilon variables (could maybe be useful one day)
e1 = matrix(c(1,0,0), nrow = 3, ncol = 1)
e2 = matrix(c(0,1,0), nrow = 3, ncol = 1)
e3 = matrix(c(0,0,1), nrow = 3, ncol = 1)
