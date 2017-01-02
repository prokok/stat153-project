library(quantmod)
dev.off()
setwd("C:/Users/vlfgn/Desktop/Stat153")

t1 = as.numeric(unlist(read.csv("2DS.csv", as.is = TRUE, header = FALSE)))
t3 = read.csv("2DS_withWeek.csv", header = FALSE)

t1.1=xts(t3$V2, as.Date(t3$V1, format = '%Y-%m-%d'))

##1. Differencing
plot(t1,type = 'l')
#plot(t1.1, type = "l")

# the resulting plot shows periodicity

#Taking log
t1.log = log(t1)
plot(t1.log, type='l')
#t1.1.log = log(t1.1)
#plot(t1.1.log, type='l')


#Taking difference to remove the trend
t1.log.d = diff(t1.log)
plot(t1.log.d, type = "l")
#t1.1.log.d = diff(t1.1.log)
#plot(t1.1.log.d, type = "l")


#Taking differnece to remove the period 
t1.log.d.dp = diff(t1.log.d, lag=52)
plot(na.omit(t1.log.d.dp), type= 'l')
#t1.1.log.d.dp = diff(t1.1.log.d, lag=52)
#plot(na.omit())
#now plot seems to be white noise


acf(na.omit(t1.log.d.dp), lag.max = 100)$acf
#This suggest that the data mgiht folloW ARMA(0,1)*(0,1)_52
pacf(na.omit(t1.log.d.dp), lag.max = 100)
#This does not reveal much information but there might ARMA(1,0)*(0,1)_52 OR ARMA(2,0)*(0,1)_52.

#However, just to be sure, let's check if there is a ar component.
####################################################################################
a1.0.1 = arima(t1.log, order=c(1,1,0), seasonal = list(order=c(0,1,1), period=52))
tsdiag(a1.0.1)
AIC(a1.0.1) #-286.3196
BIC(a1.0.1) #-277.1701
MSE.a1.0.1 = computeCVmse(c(1,1,0),c(0,1,1))
mean(MSE.a1.0.1) #134.8648
MSE2.a1.0.1 = computeCVmse2(c(1,1,0),c(0,1,1))
mean(na.omit(MSE2.a1.0.1)) #85.50809
tail(MSE2.a1.0.1) #23.10114 33.85171 33.39530 55.16469 56.37186 31.43467

a1.1.1 = arima(t1.log, order=c(1,1,1), seasonal = list(order=c(0,1,1),period=52))
tsdiag(a1.1.1)
AIC(a1.1.1) #-304.3865
BIC(a1.1.1) #-292.1871
MSE.a1.1.1 = computeCVmse(c(1,1,1),c(0,1,1))
mean(MSE.a1.1.1) #203.4652
MSE2.a1.1.1 = computeCVmse2(c(1,1,1),c(0,1,1))
mean(na.omit(MSE2.a1.1.1)) #59.30519
tail(MSE2.a1.1.1) #36.573520 42.631222 44.454132 63.162130 72.829159  4.554825

a2.1.1 = arima(t1.log, order=c(2,1,1), seasonal = list(order=c(0,1,1),period=52))
tsdiag(a2.1.1)
AIC(a2.1.1) #-292.9847
BIC(a2.1.1) #-280.7853
MSE.a2.1.1 = computeCVmse(c(2,1,1),c(0,1,1))
mean(MSE.a2.1.1) #102.4728
MSE2.a2.1.1 = computeCVmse2(c(2,1,1),c(0,1,1))
mean(na.omit(MSE2.a2.1.1)) #57.78728
tail(MSE2.a2.1.1) #38.890615 43.228464 44.153853 62.478671 71.806414  3.984438

a2.1.0 = arima(t1.log, order=c(2,1,0), seasonal = list(order=c(0,1,1),period=52))
tsdiag(a2.1.0)
AIC(a2.1.0) #-398.9074
BIC(a2.1.0) #-383.6581
MSE.a2.1.0 = computeCVmse(c(2,1,0),c(0,1,1))
mean(MSE.a2.1.0) #105.3383
MSE2.a2.1.0 = computeCVmse2(c(2,1,0),c(0,1,1))
mean(na.omit(MSE2.a2.1.0)) #70.15789
tail(MSE2.a2.1.0) #34.37890 33.61641 35.22100 55.73129 59.29191 13.43151


##If there's a ar componenet a1.1.1 is the best model.
##########################################################################################

m1.1 = arima(t1.log, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1), period = 52))
m1.2 = arima(t1.log, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 2), period = 52))
m1.3 = arima(t1.log, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 3), period = 52))

tsdiag(m1.1)
tsdiag(m1.2)
tsdiag(m1.3)

BIC(m1.1) #-297.1351
BIC(m1.2) #-292.38
BIC(m1.3) #-287.383

AIC(m1.1) #-306.2846
AIC(m1.2) #-304.5794
AIC(m1.3) #-302.6323
#Result: Among these models, m1.1 is selected as a best model.

m2.1 = arima(t1.log, order = c(0, 1, 2), seasonal = list(order = c(0, 1, 1), period = 52))
m3.1 = arima(t1.log, order = c(0, 1, 3), seasonal = list(order = c(0, 1, 1), period = 52))
m4.1 = arima(t1.log, order = c(0, 1, 4), seasonal = list(order = c(0, 1, 1), period = 52))
m5.1 = arima(t1.log, order = c(0, 1, 5), seasonal = list(order = c(0, 1, 1), period = 52))
m6.1 = arima(t1.log, order = c(0, 1, 6), seasonal = list(order = c(0, 1, 1), period = 52))

tsdiag(m2.1)
tsdiag(m3.1)
tsdiag(m4.1)
tsdiag(m5.1)
tsdiag(m6.1)

BIC(m2.1) #-292.2093
BIC(m3.1) #-287.7974
BIC(m4.1) #-287.2553
BIC(m5.1) #-283.2242
BIC(m6.1) #-286.6051
#AIC would choose m2.1.

AIC(m2.1) #-304.4087
AIC(m3.1) #-303.0467
AIC(m4.1) #-305.5544
AIC(m5.1) #-304.5732
AIC(m6.1) #-311.0039
#BIC would choose m6.1.
#Result: Among these models, m2.1, m6.1 are selected as preferable models.

m2.2 = arima(t1.log, order = c(0, 1, 2), seasonal = list(order = c(0, 1, 2), period = 52))
m2.3 = arima(t1.log, order = c(0, 1, 2), seasonal = list(order = c(0, 1, 3), period = 52))
m2.4 = arima(t1.log, order = c(0, 1, 2), seasonal = list(order = c(0, 1, 4), period = 52))

tsdiag(m2.2)
tsdiag(m2.3)
tsdiag(m2.4)
#According to Ljung Box statistic, above three methods show no difference.

BIC(m2.2)#-287.4233
BIC(m2.3)#-282.4086
BIC(m2.4)#-277.3587
#BIC would choose m2.2

AIC(m2.2)#-302.6726
AIC(m2.3)#-300.7077
AIC(m2.4)#-298.7077
#AIC would choose m2.2
#Result: Among these models, m2.2 is selected as a best model.

#What is the best model among m1.1, m2.1, m6.1, m2.2
##Let's use CV to find out 

len = length(t1.log)

computeCVmse = function(order.totry,seasorder.totry){
  MSE = numeric()
  for(k in 2:1){
    train.dt = t1.log[1:(len-52*k)]
    test.dt = t1.log[(len-52*k+1):(len - 52*(k-1))]
    mod = arima(train.dt, order = order.totry, seasonal = 
                  list(order = seasorder.totry, period = 52))
    fcast = predict(mod, n.ahead = 52)
    MSE[k] = mean((exp(fcast$pred) - exp(test.dt))^2)
  }
  return(MSE)
}

computeCVmse2 = function(order.totry, seasorder.totry){
  MSE=numeric(0)
  for(k in 108:208){
    train.dt = t1.log[1:k]
    test.dt = t1.log[(k+1):len]
    mod = arima(train.dt, order = order.totry, seasonal = list(order=seasorder.totry, period = 52))
    fcast = predict(mod, n.ahead = (len-k))
    MSE[k] = mean((exp(fcast$pred) - exp(test.dt))^2)
  }
  return(MSE)
}


##computeCVmse function output
MSE.m1.1 = computeCVmse(c(0, 1, 1), c(0,1,1)) #19.53087 159.42724
MSE.m2.1 = computeCVmse(c(0, 1, 2), c(0,1,1)) #19.53459 265.75680
MSE.m2.2 = computeCVmse(c(0, 1, 2), c(0,1,2)) #18.57876 265.69014
MSE.m6.1 = computeCVmse(c(0, 1, 6), c(0,1,1)) #21.21828 404.80785

##computeCVmse2 function output
MSE2.m1.1 = computeCVmse2(c(0, 1, 1), c(0,1,1))
MSE2.m2.1 = computeCVmse2(c(0, 1, 2), c(0,1,1))
MSE2.m2.2 = computeCVmse2(c(0, 1, 2), c(0,1,2))
MSE2.m6.1 = computeCVmse2(c(0, 1, 6), c(0,1,1))

##See the tails of MSE2
tail(MSE2.m1.1) #
tail(MSE2.m2.1) #
tail(MSE2.m2.2) #
tail(MSE2.m6.1) #


#Average CV-score using computeCVmse
avg.m1.1 = mean(MSE.m1.1) # 89.47905
avg.m2.1 = mean(MSE.m2.1) # 142.6457
avg.m2.2 = mean(MSE.m2.2) # 142.1345
avg.m6.1 = mean(MSE.m6.1) # 213.0131
##Result : m1.1 model is the best model according to computeCVmse.

#Average CV-score using computeCVmse2
avg2.m1.1 = mean(na.omit(MSE2.m1.1)) #57.09834
avg2.m2.1 = mean(na.omit(MSE2.m2.1)) #62.61343
avg2.m2.2 = mean(na.omit(MSE2.m2.2)) #63.43211
avg2.m6.1 = mean(na.omit(MSE2.m6.1)) #47.34336
##Result : m6.1 models are the best model according to computeCVmse2.

######################################################################################
##Comparing the best candidates using the AIC, BIC, CV results
##m1.1
AIC(m1.1) #-306.2846
BIC(m1.1) #-306.2846
avg.m1.1  # 89.47905
avg2.m1.1 # 57.09834
MSE.m1.1 #19.53087 159.42724
tail(MSE2.m1.1) #36.563576 42.381422 44.023080 62.506252 73.084922  3.907559

##m6.1
AIC(m6.1) #-311.0039
BIC(m6.1) #-286.6051
avg.m6.1  # 213.0131
avg2.m6.1 # 47.34336
MSE.m6.1 #21.21828 404.80785
tail(MSE2.m6.1) #45.12959 54.81341 44.17570 53.27392 74.29521 10.17254

######################################################################################
#Plots
plot(1:length(t1), t1, type = "l", xlab = "Weekly Time", xlim= c(0,length(t1)+52), ylim=c(0,110))
lines((length(t1)+1):(length(t1)+52), exp(predict(m6.1, 52)$pred), col = "red")
lines(x = c(length(t1),length(t1)+1),y = c(t1[length(t1)], exp(predict(m6.1,1)$pred))
      , type = 'l', col="red")

######################################################################################

predictions <- exp(predict(m6.1, n.ahead = 52)$pred)

write.table(predictions,
            sep = "\n",
            col.names = FALSE,
            row.names = FALSE,
            file = "Q2_PHILHOON_OH_23997722.txt")
##########################################################################################################

