library(quantmod)
dev.off()
setwd("C:/Users/vlfgn/Desktop/Stat153")

t1 = as.numeric(unlist(read.csv("1DS.csv", as.is = TRUE, header = FALSE)))
t3 = read.csv("1DS_withWeek.csv", header = FALSE)

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
AIC(a1.0.1) #-369.8849
BIC(a1.0.1) #-360.7353
MSE.a1.0.1 = computeCVmse(c(1,1,0),c(0,1,1))
mean(MSE.a1.0.1) #25.23343
MSE2.a1.0.1 = computeCVmse2(c(1,1,0),c(0,1,1))
mean(na.omit(MSE2.a1.0.1)) #26.74823
tail(MSE2.a1.0.1) #10.74570 20.20786 32.21813 21.75386 19.31622 24.36642

a1.1.1 = arima(t1.log, order=c(1,1,1), seasonal = list(order=c(0,1,1),period=52))
tsdiag(a1.1.1)
AIC(a1.1.1) #-400.6762
BIC(a1.1.1) #-388.4768
MSE.a1.1.1 = computeCVmse(c(1,1,1),c(0,1,1))
mean(MSE.a1.1.1) #17.13972
MSE2.a1.1.1 = computeCVmse2(c(1,1,1),c(0,1,1))
mean(na.omit(MSE2.a1.1.1)) #11.94141
tail(MSE2.a1.1.1) #7.604789  5.493301  5.206586  6.927420 10.482986 19.171434

a2.1.1 = arima(t1.log, order=c(2,1,1), seasonal = list(order=c(0,1,1),period=52))
tsdiag(a2.1.1)
AIC(a2.1.1) #-398.9074
BIC(a2.1.1) #-383.6581
MSE.a2.1.1 = computeCVmse(c(2,1,1),c(0,1,1))
mean(MSE.a2.1.1) #19.381
MSE2.a2.1.1 = computeCVmse2(c(2,1,1),c(0,1,1))
mean(na.omit(MSE2.a2.1.1)) #12.33244
tail(MSE2.a2.1.1) #7.527082  5.635878  5.656881  7.472509 11.172889 19.956677

a2.1.0 = arima(t1.log, order=c(2,1,0), seasonal = list(order=c(0,1,1),period=52))
tsdiag(a2.1.0)
AIC(a2.1.0) #-398.9074
BIC(a2.1.0) #-383.6581
MSE.a2.1.0 = computeCVmse(c(2,1,0),c(0,1,1))
mean(MSE.a2.1.0) #19.381
MSE2.a2.1.0 = computeCVmse2(c(2,1,0),c(0,1,1))
mean(na.omit(MSE2.a2.1.0)) #12.33244
tail(MSE2.a2.1.0) #7.527082  5.635878  5.656881  7.472509 11.172889 19.956677


##If there's a ar componenet a1.1.1 is the best model.
##########################################################################################

m1.1 = arima(t1.log, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1), period = 52))
m1.2 = arima(t1.log, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 2), period = 52))
m1.3 = arima(t1.log, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 3), period = 52))

tsdiag(m1.1)
tsdiag(m1.2)
tsdiag(m1.3)
#Ljung-Box statistic shows that these models are inappropriate

BIC(m1.1) #-374.4202
BIC(m1.2) #-370.1097
BIC(m1.3) #-365.4244
#BIC would choose m1.1.

AIC(m1.1) #-383.5697
AIC(m1.2) #-382.3091
AIC(m1.3) #-380.6737
#AIC would choose m1.1 is the most preferred model 
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
#Ljung Box statistic, m4.1, m5.1, m6.1 are selected.

BIC(m2.1) #-384.9544
BIC(m3.1) #-381.5174
BIC(m4.1) #-382.1302
BIC(m5.1) #-377.4513
BIC(m6.1) #-372.5197
#AIC would choose m2.1.

AIC(m2.1) #-397.1538
AIC(m3.1) #-396.7667
AIC(m4.1) #-400.4293
AIC(m5.1) #-398.8003
AIC(m6.1) #-396.9185
#BIC would choose m4.1.
#Result: Among these models, m4.1, m2.1, m5.1, m6.1 are selected as preferable models.

m2.2 = arima(t1.log, order = c(0, 1, 2), seasonal = list(order = c(0, 1, 2), period = 52))
m2.3 = arima(t1.log, order = c(0, 1, 2), seasonal = list(order = c(0, 1, 3), period = 52))
m2.4 = arima(t1.log, order = c(0, 1, 2), seasonal = list(order = c(0, 1, 4), period = 52))

tsdiag(m2.2)
tsdiag(m2.3)
tsdiag(m2.4)
#According to Ljung Box statistic, above three methods show no difference.

BIC(m2.2)#-382.4255
BIC(m2.3)#-377.3767
BIC(m2.4)#-373.1674
#BIC would choose m2.2

AIC(m2.2)#-397.6748
AIC(m2.3)#-395.6758
AIC(m2.4)#-394.5164
#AIC would choose m2.2
#Result: Among these models, m2.2 is selected as a best model.

#What is the best model among m1.1, m4.1, m2.1, m2.2, m5.1, m6.1 
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
MSE.m1.1 = computeCVmse(c(0, 1, 1), c(0,1,1))
MSE.m4.1 = computeCVmse(c(0, 1, 4), c(0,1,1))
MSE.m2.1 = computeCVmse(c(0, 1, 2), c(0,1,1))
MSE.m2.2 = computeCVmse(c(0, 1, 2), c(0,1,2))
MSE.m5.1 = computeCVmse(c(0, 1, 5), c(0,1,1))
MSE.m6.1 = computeCVmse(c(0, 1, 6), c(0,1,1))

##computeCVmse2 function output
MSE2.m1.1 = computeCVmse2(c(0, 1, 1), c(0,1,1))
MSE2.m4.1 = computeCVmse2(c(0, 1, 4), c(0,1,1))
MSE2.m2.1 = computeCVmse2(c(0, 1, 2), c(0,1,1))
MSE2.m2.2 = computeCVmse2(c(0, 1, 2), c(0,1,2))
MSE2.m5.1 = computeCVmse2(c(0, 1, 5), c(0,1,1))
MSE2.m6.1 = computeCVmse2(c(0, 1, 6), c(0,1,1))

##See the tails of MSE2
tail(MSE2.m1.1) #
tail(MSE2.m4.1) #
tail(MSE2.m2.1) #
tail(MSE2.m2.2) #
tail(MSE2.m5.1) #
tail(MSE2.m6.1) #


#Average CV-score using computeCVmse
avg.m1.1 = mean(MSE.m1.1) # 20.15904
avg.m4.1 = mean(MSE.m4.1) # 17.49653
avg.m2.1 = mean(MSE.m2.1) # 19.11847
avg.m2.2 = mean(MSE.m2.2) # 27.51723
avg.m5.1 = mean(MSE.m5.1) # 17.32303
avg.m6.1 = mean(MSE.m6.1) # 18.79556
##Result : m4.1 model is the best model according to computeCVmse.

#Average CV-score using computeCVmse2
avg2.m1.1 = mean(na.omit(MSE2.m1.1)) #18.37231
avg2.m4.1 = mean(na.omit(MSE2.m4.1)) #12.12896
avg2.m2.1 = mean(na.omit(MSE2.m2.1)) #12.89052
avg2.m2.2 = mean(na.omit(MSE2.m2.2)) #13.03457
avg2.m5.1 = mean(na.omit(MSE2.m5.1)) #12.06189
avg2.m6.1 = mean(na.omit(MSE2.m6.1)) #12.14805
##Result : m4.1, m1.5 models are the best model according to computeCVmse2.

######################################################################################
##Comparing the best candidates using the AIC, BIC, CV results
##m4.1 
AIC(m4.1) #-400.4293
BIC(m4.1) #-382.1302
avg.m4.1  # 17.49653
avg2.m4.1 # 12.12896
MSE.m4.1 #11.95002 23.04304
tail(MSE2.m4.1) #7.263629  5.404993  5.678071  7.438003 11.010460 18.938550

##m5.1 
AIC(m5.1) #-398.8003
BIC(m5.1) #-377.4513
avg.m5.1  # 17.32303
avg2.m5.1 # 12.06189
MSE.m5.1 #11.94706 22.69901
tail(MSE2.m5.1) #7.141472  5.400059  5.169943  6.926618 10.369848 18.460622

##a1.1.1
AIC(a1.1.1) #-400.6762
BIC(a1.1.1) #-388.4768
mean(MSE.a1.1.1) #17.13972
mean(na.omit(MSE2.a1.1.1)) #11.94141
MSE.a1.1.1 #11.89859 22.38086
tail(MSE2.a1.1.1) #7.604789  5.493301  5.206586  6.927420 10.482986 19.1714347

######################################################################################
#Plots
plot(1:length(t1), t1, type = "l", xlab = "Weekly Time", xlim= c(0,length(t1)+52))
lines((length(t1)+1):(length(t1)+52), exp(predict(a1.1.1, 52)$pred), col = "red")
lines(x = c(length(t1),length(t1)+1),y = c(t1[length(t1)], exp(predict(a1.1.1,1)$pred))
      , type = 'l', col="red")

######################################################################################

predictions <- exp(predict(a1.1.1, n.ahead = 52)$pred)

write.table(predictions,
            sep = "\n",
            col.names = FALSE,
            row.names = FALSE,
            file = "Q1_PHILHOON_OH_23997722.txt")
##########################################################################################################

