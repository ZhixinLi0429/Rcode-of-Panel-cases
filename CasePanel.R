########################
######panel case########
########################
setwd("C:\\Users\\PaulRomer\\Desktop\\Rworkfile\\Panel") 
if(!require(tseries)){install.packages("tseries")}
if(!require(LSTS)){install.packages("LSTS")}
if(!require(knitr)){install.packages("knitr")}
if(!require(stargazer)){install.packages("stargazer")}
if(!require(xtable)){install.packages("xtable")}
if(!require(rmarkdown)){install.packages("rmarkdown")}
if(!require(plm)){install.packages("plm")}
if(!require(gplots)){install.packages("gplots")}

## required packages
library(stargazer)
library(xtable)
library(knitr)
library(tools)
library(rmarkdown)
library(plm)
library(gplots)
source("TimeSeriesFunctions.R")

panel_data <- read.table("C:/Users/PaulRomer/Desktop/Rworkfile/Panel/deleted_paneldata.csv", header = TRUE,sep = ",")
panel_dataFrame = pdata.frame(panel_data, index = c ("CountryID","TimeID"))
View(panel_dataFrame)

CountryID <- panel_dataFrame$CountryID
TimeID    <- panel_dataFrame$TimeID
lnGDP     <- panel_dataFrame$lnGDP
D         <- panel_dataFrame$D
class(lnGDP)

###data handling###
GDPgrowth <- diff(lnGDP,differences = 1)
panel_dataFrame[, "GDPgrowth"] <- GDPgrowth

###
GDPgrowth_lag4Matrix <- data.frame(matrix(NA, nrow = nrow(panel_data), ncol = 4))
for (i in 1:4) {
  GDPgrowth_lag4Matrix[, i] <- c(rep(NA, i), lag(GDPgrowth, i)[1:(nrow(panel_data) - i)])
}
colnames(GDPgrowth_lag4Matrix) <- paste0("GDPgrowth_lag", 1:4)
panel_dataFrame <- cbind(panel_dataFrame, GDPgrowth_lag4Matrix)

###
lnGDP_lag4Matrix <- data.frame(matrix(NA, nrow = nrow(panel_data), ncol = 4))
for (i in 1:4) {
  lnGDP_lag4Matrix[, i] <- c(rep(NA, i), lag(lnGDP, i)[1:(nrow(panel_data) - i)])
}
colnames(lnGDP_lag4Matrix) <- paste0("lnGDP_lag", 1:4)
panel_dataFrame <- cbind(panel_dataFrame, lnGDP_lag4Matrix)

############################################
###static panel data model,GDP in levels###
###########################################

POLS = plm(lnGDP ~ D, data = panel_dataFrame, model = "pooling")
summary(POLS)

FE   = plm(lnGDP ~ D, data = panel_dataFrame, model = "within")
summary(FE)

RE   = plm(lnGDP ~ D, data = panel_dataFrame, model = "random")
summary(RE)

pFtest(FE, POLS)
#=> fe

phtest(FE,RE)
#=> fe

#breusch-Godfrey-Wooldridge Autocorreltion test
pbgtest(FE, order = 1)



############################################
###static panel data model,GDP in Growth###
###########################################


POLS_inGrowth = plm(GDPgrowth ~ D, data = panel_dataFrame, model = "pooling")
summary(POLS_inGrowth)

FE_inGrowth   = plm(GDPgrowth ~ D, data = panel_dataFrame, model = "within")
summary(FE_inGrowth)

RE_inGrowth   = plm(GDPgrowth ~ D, data = panel_dataFrame, model = "random")
summary(RE_inGrowth)


pFtest(FE_inGrowth, POLS_inGrowth)

phtest(FE_inGrowth,RE_inGrowth)


#breusch-Godfrey-Wooldridge Autocorreltion test
pbgtest(RE_inGrowth, order = 1)


######################################
###dynamic model for GDP in levels###
#####################################

RE_residual <- RE$residuals
acf(RE_residual)# why we plot the acf of the residual here? what's the intuition behind?
FE_residual <- FE$residuals
acf(FE_residual)#完全没必要干这个
lag1 <- panel_dataFrame$lnGDP_lag1
lag2 <- panel_dataFrame$lnGDP_lag2
lag3 <- panel_dataFrame$lnGDP_lag3
lag4 <- panel_dataFrame$lnGDP_lag4


POLS__inlevel_lag <- plm(lnGDP ~ D+lag1+lag2+lag3+lag4, data = panel_dataFrame, model = "pooling")
summary(POLS__inlevel_lag)

FE_inlevel_lag <- plm(lnGDP ~ D+lag1+lag2+lag3+lag4, data = panel_dataFrame, model = "within")
summary(FE_inlevel_lag)

RE_inlevel_lag <- plm(lnGDP ~ D+lag1+lag2+lag3+lag4, data = panel_dataFrame, model = "random")
summary(RE_inlevel_lag)


pFtest(FE_inlevel_lag, POLS__inlevel_lag)


phtest(FE_inlevel_lag, RE_inlevel_lag)

pbgtest(FE_inlevel_lag, order = 1)


FE_inlevel_lag_adjust <- plm(lnGDP ~ D+lag1+lag2, data = panel_dataFrame, model = "within")
summary(FE_inlevel_lag_adjust)

pbgtest(FE_inlevel_lag_adjust, order = 1)


TFE_inlevel_lag_adjust <- plm(lnGDP ~ D+lag1+lag2, data = panel_dataFrame, model = "within", effect= "twoways")
summary(TFE_inlevel_lag_adjust)

pFtest(TFE_inlevel_lag_adjust, FE_inlevel_lag_adjust)

pbgtest(TFE_inlevel_lag_adjust, order = 1)


TFE_inlevel_lag_adjust2 <- plm(lnGDP ~ D+lag1+lag2+lag3+lag4, data = panel_dataFrame, model = "within", effect= "twoways")
summary(TFE_inlevel_lag_adjust2)

pbgtest(TFE_inlevel_lag_adjust2, order = 1)

CCEP<- pcce(lnGDP ~ D+lag1+lag2, data = panel_dataFrame, model = "p")
summary(CCEP)

pbgtest(CCEP, order = 1)

CCEMG<- pcce(lnGDP ~ D+lag1+lag2, data = panel_dataFrame, model = "mg")
summary(CCEMG)

pbgtest(CCEMG, order = 1)

######################################
###dynamic model for GDPgrowth    ###
#####################################

lag_1 <- panel_dataFrame$GDPgrowth_lag1
lag_2 <- panel_dataFrame$GDPgrowth_lag2
lag_3 <- panel_dataFrame$GDPgrowth_lag3
lag_4 <- panel_dataFrame$GDPgrowth_lag4


POLS_GDPgrowth_lag <- plm(GDPgrowth ~ D+lag_1+lag_2, data = panel_dataFrame, model = "pooling")
summary(POLS_GDPgrowth_lag)

FE_GDPgrowth_lag <- plm(GDPgrowth ~ D+lag_1+lag_2, data = panel_dataFrame, model = "within")
summary(FE_GDPgrowth_lag)

RE_GDPgrowth_lag <- plm(GDPgrowth ~ D+lag_1+lag_2, data = panel_dataFrame, model = "random")
summary(RE_GDPgrowth_lag)

pFtest(FE_GDPgrowth_lag, POLS_GDPgrowth_lag)

phtest(FE_GDPgrowth_lag, RE_GDPgrowth_lag)

pbgtest(FE_GDPgrowth_lag)

FE_GDPgrowth_lag_adjust <- plm(GDPgrowth ~ D+lag_1+lag_2+lag_3+lag_4, data = panel_dataFrame, model = "within")
summary(FE_GDPgrowth_lag_adjust)

pbgtest(FE_GDPgrowth_lag_adjust)


TFE_GDPgrowth_lag_adjust <- plm(GDPgrowth ~ D+lag1+lag2, data = panel_dataFrame, model = "within", effect= "twoways")
summary(TFE_GDPgrowth_lag_adjust)

pbgtest(TFE_GDPgrowth_lag_adjust)




#######################
###local projections###
#######################

install.packages("lpirfs")
library(lpirfs)

DeltaD <- diff(D,differences = 1)
panel_dataFrame[, "DeltaD"] <- DeltaD

lp_GDPlong_growth_DeltaD <- lp_lin_panel(panel_dataFrame, 
                                         endog_data = "GDPgrowth", shock = "DeltaD", l_exog_data = "GDPgrowth", 
                                         lags_exog_data = 2, diff_shock = FALSE, confint = 1.96, hor = 30)

plot(lp_GDPlong_growth_DeltaD)




lp_GDP_DeltaD <- lp_lin_panel(panel_dataFrame, 
                                         endog_data = "lnGDP", shock = "DeltaD", l_exog_data = "lnGDP", 
                                         lags_exog_data = 2, diff_shock = FALSE, confint = 1.96, hor = 30)

plot(lp_GDP_DeltaD)





#######################################################
### Panel unit root test and cointegration analysis ###
#######################################################


purtest(panel_dataFrame$lnGDP, exo= "trend",test = "levinlin")
purtest(panel_dataFrame$lnGDP, exo= "trend", test = "ips",index = c ("CountryID","TimeID"))
purtest(panel_dataFrame$lnGDP,pmax = 4, exo = "trend", test = "madwu")



purtest(panel_dataFrame[, "GDPgrowth"],exo = "intercept", test = "levinlin")
purtest(panel_dataFrame[, "GDPgrowth"], test = "ips", exo = "intercept")
purtest(panel_dataFrame[, "GDPgrowth"],pmax = 4, exo = "intercept", test = "madwu")





pcdtest(lnGDP)
pcdtest(GDPgrowth)



CCEP_forcointegration<- pcce(lnGDP ~ D, data = panel_dataFrame, model = "p")
summary(CCEP_forcointegration)

residuals_cce <- residuals(CCEP_forcointegration)
residuals_df_cce <- data.frame(individual_id=panel_dataFrame$CountryID, time_id=panel_dataFrame$TimeID, 
                               residuals=residuals_cce)
paneldata_residuals_cce <- pdata.frame(residuals_df_cce, index = c ("individual_id","time_id"))

purtest(paneldata_residuals_cce$residuals, test = "madwu")
