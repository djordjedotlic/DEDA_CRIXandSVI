rm(list = ls(all = TRUE))
graphics.off()

if(!require("strucchange")) install.packages("strucchange"); library("strucchange")
if(!require("forecast")) install.packages("forecast");library("forecast")
if(!require("lmtest")) install.packages("lmtest");library("lmtest")
if(!require("ggplot2")) install.packages("ggplot2");library("ggplot2")
if(!require("tseries")) install.packages("tseries");library("tseries")


#Loading the data----
input <- read.csv("data/df.csv", header = TRUE, sep = ",")
input$date <- as.Date(input$date, "%m/%d/%Y")
input$X.1 <- NULL

ggplot( data = input, aes( date, svi )) + geom_line() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                              panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggplot( data = input, aes( date, crix )) + geom_line() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                               panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggplot( data = input, aes( date, returns )) + geom_line() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                  panel.background = element_blank(), axis.line = element_line(colour = "black"))




#Detecting structural breaks in the time series----
#Adapted code from: https://www.r-bloggers.com/endogenously-detecting-structural-breaks-in-a-time-series-implementation-in-r/

ts_crix <- as.ts(input$crix)
ts_svi <- as.ts(input$svi)

bp_ts_crix <- breakpoints(ts_crix~1)
bp_ts_svi <- breakpoints(ts_svi~1)


summary(bp_ts_crix)
summary(bp_ts_svi)

p <- array()
p[1] <- bp_ts_svi$breakpoints[1]
p[2] <- bp_ts_svi$breakpoints[2]
input$date[p] #dates when structural breaks happened in svi time series

c <- array()
c[1] <- bp_ts_crix$breakpoints[1]
c[2] <- bp_ts_crix$breakpoints[2]
input$date[c] #dates when structural breaks happened in svi time series

#plot the breakpoints 
bp_date_crix <- array()
bp_date_crix[1] <- input$date[c[1]]
bp_date_crix[1] <- as.Date(bp_date_crix[1], format = "%Y-%m-%d")
bp_date_crix[2] <- input$date[c[2]]
bp_date_crix[2] <- as.Date(bp_date_crix[2], format = "%Y-%m-%d")
ggplot( data = input, aes( date, crix )) + geom_line() + geom_vline(xintercept = c(bp_date_crix))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                         panel.background = element_blank(), axis.line = element_line(colour = "black"))



bp_date <- array()
bp_date[1] <- input$date[p[1]]
bp_date[1] <- as.Date(bp_date[1], format = "%Y-%m-%d")
bp_date[2] <- input$date[p[2]]
bp_date[2] <- as.Date(bp_date[2], format = "%Y-%m-%d")
ggplot( data = input, aes( date, svi )) + geom_line() + geom_vline(xintercept = c(bp_date))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                                   panel.background = element_blank(), axis.line = element_line(colour = "black"))


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


granger.caus_p.values <- matrix(nrow = 6, ncol = 2)
gc_results <- vector(mode = "list",length = 3)


for(l in 1:3){
  input <- input_list[[l]]
  
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
 
  
  for(i in 1:6){
    gt1 <- grangertest(dreturns ~ dsvi, order= i)
    granger.caus_p.values[i,1] <- gt1$`Pr(>F)`[2]
    
    gt2 <- grangertest(dsvi ~ dreturns, order=i)
    granger.caus_p.values[i,2] <- gt2$`Pr(>F)`[2]
    
  }

  gc_results[[l]] <- granger.caus_p.values
  
}

gc_results 
#Gives p-values for granger causality tests, 
#across the number of lags, across 3 periods


#Regression analysis----

input <- read.csv("data/df.csv", header = TRUE, sep = ",")
input$date <- as.Date(input$date, "%m/%d/%Y")
input$X.1 <- NULL
input$crix <- NULL
input$svi <- log(input$svi)
input$volume <- log(input$volume)
input$epuix <- log(input$epuix)
colnames(input)[3:9] <- c("log(svi)","log_returns", "log_sp500_returns",
                     "log(vix)", "log(volume)", "volatility", "log(epuix)")

input_1 <- input[1:91,]
input_2 <- input[92:124,]
input_3 <- input[125:146,]
input_list <- list()
input_list[[1]] <- input_1
input_list[[2]] <- input_2
input_list[[3]] <- input_3

pred_input <- read.csv("data/pred_df.csv", header = TRUE, sep = ",")
colnames(pred_input)[3:9] <- c("log(svi)","log_returns", "log_sp500_returns",
                               "log(vix)", "log(volume)", "volatility", "log(epuix)")

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
  cont_reg <- lm(log_returns~., data = input[,-(1:2)])
  cont_reg_list[[i]] <- cont_reg
  bp_cont[[i]] <- bptest(cont_reg)
  #Predictive regression
  pred_reg <- lm(log_returns~., data = pred_input[,-(1:2)])
  pred_reg_list[[i]] <- pred_reg
  bp_pred[[i]] <- bptest(pred_reg)
  
}
