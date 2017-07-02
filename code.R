if(!require("strucchange")) install.packages("strucchange"); library("strucchange")
if(!require("forecast")) install.packages("forecast");library("forecast")
if(!require("lmtest")) install.packages("lmtest");library("lmtest")
if(!require("ggplot2")) install.packages("ggplot2");library("ggplot2")
if(!require("tseries")) install.packages("tseries");library("tseries")


#Loading the data----
input <- read.csv("data/df.csv", header = TRUE, sep = ",")
input$date <- as.Date(input$date, "%m/%d/%Y")
input$X.1 <- NULL


#Detecting structural breaks in the time series----
#Adapted code from: https://www.r-bloggers.com/endogenously-detecting-structural-breaks-in-a-time-series-implementation-in-r/

ts_crix <- as.ts(input$crix)
ts_svi <- as.ts(input$svi)

bp_ts_crix <- breakpoints(ts_crix~1)
bp_ts_svi <- breakpoints(ts_svi~1)


summary(bp_ts_crix)
summary(bp_ts_svi)

#plot the breakpoints 
plot(input$date, input$crix, xaxt = "n")
axis.Date(side = 1, input$date, format = "%d/%m/%Y")
lines(bp_ts_crix)


plot(input$date, input$svi, xaxt = "n", type="n")
axis.Date(side = 1, input$date, format = "%d/%m/%Y")
lines(input$date, input$svi)
lines(bp_ts_svi)


input$crix <- NULL

#Granger causality test----
#Adapted code from: https://www.r-bloggers.com/chicken-or-the-egg-granger-causality-for-the-masses/


#Here, we will separate ts into three series, according to the structural breaks
#in svi time series


input_1 <- input[1:91,]
input_2 <- input[92:124,]
input_3 <- input[125:146,]
input_list <- list()
input_list[[1]] <- input_1
input_list[[2]] <- input_2
input_list[[3]] <- input_3


for(i in 1:3){
  input <- input_list[[i]]
  
  #Checking for the stationarity of the series
  ts_returns <- ts(input$returns)
  ndiffs (ts_returns, alpha=0.05, test= "kpss")
  ndiffs (ts_returns, alpha=0.05, test= "adf")
  dreturns <- diff(ts_returns)
  
  
  ts_svi <- ts(input$svi)
  ndiffs(ts_svi, alpha=0.05, test= "kpss")
  ndiffs(ts_svi, alpha=0.05, test= "adf")
  dsvi <- diff(ts_svi)
  
  
  #testing for granger causality
  #the output is a table for each of the periods 
  #identified by the structural break test, respectively 
  granger.caus_p.values <- matrix(nrow = 9, ncol = 2)
  colnames(granger.caus_p.values) <- c("svi_runs_returns", "returns_run_svi")
  
  for(i in 1:9){
    gt1 <- grangertest(dreturns ~ dsvi, order=i)
    granger.caus_p.values[i,1] <- gt1$`Pr(>F)`[2]
    
    gt2 <- grangertest(dsvi ~ dreturns, order=i)
    granger.caus_p.values[i,2] <- gt2$`Pr(>F)`[2]
    
  }
  
  View (granger.caus_p.values)
  
  
}


#Regression analysis----
input$svi <- log(input$svi)
input$volume <- log(input$volume)
input$epuix <- log(input$epuix)

pred_input <- read.csv("data/pred_df.csv", header = TRUE, sep = ",")
pred_input_1 <- pred_input[1:91,]
pred_input_2 <- pred_input[92:124,]
pred_input_3 <- pred_input[125:146,]
pred_input_list <- list()
pred_input_list[[1]] <- pred_input_1
pred_input_list[[2]] <- pred_input_2
pred_input_list[[3]] <- pred_input_3

cont_reg_list <- vector(mode="list")
pred_reg_list <- vector(mode="list")
bp_cont <- vector(mode="list")
bp_pred <- vector(mode="list")

for(i in 1:3){
  input <- input_list[[i]]
  pred_input <- pred_input_list[[i]]
  #Contempraneous regression
  cont_reg <- lm(returns~., data = input[,-(1:2)])
  cont_reg_list[[i]] <- cont_reg
  bp_cont[[i]] <- bptest(cont_reg)
  #Predictive regression
  pred_reg <- lm(returns~., data = pred_input[,-(1:2)])
  pred_reg_list[[i]] <- pred_reg
  bp_pred[[i]] <- bptest(pred_reg)
  
}