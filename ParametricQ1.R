library(quantmod)
library(TSA)
library(zoo)
setwd("C:/Users/vlfgn/Desktop/Stat153")
dev.off()
t3 = read.csv("1DS_withWeek.csv", header = FALSE)
t1.1 = ts(t3[,2], start = 1, frequency = 52)

z = time(t1.1)
plot(t1.1)

par(mfrow=c(1,2))

lmod=lm(t1.1~z)
abline(lmod, col = 'blue')
resl = lmod$residuals
plot(resl,type = 'l')
acf(resl,lag.max=100) 
pacf(resl,lag.max=100)

har=harmonic(t1.1,m=20)
smod=lm(resl~har)
plot(resl,type='l')
lines(smod$fitted.values,col='red')

par(mfrow=c(1,2))
acf(smod$residuals,lag.max=100)
pacf(smod$residuals, lag.max=100)

MSE.full=numeric(0)
for(k in 1:26){ 
  har=harmonic(t1.1,k)
  
  smod=lm(resl~har)
  plot(resl,type='l', main = paste0("k value ",k))
  lines(smod$fitted.values,col='red')
  acf(smod$residuals,lag.max=100) 
  #temp=readline(prompt='press key to continue') #allows pause before new graph is displayed
  predict = lmod$fitted.values + smod$fitted.values
  MSE.full[k] = mean((predict - t1.1)^2)
}

dev.off()
tp=5+(1:52)/52 #prediction times, since data ended at time 5
sp=lmod$coefficients[1]+lmod$coefficients[2]*tp+smod$fitted.values[2:53]
#1st part does a+b*t, then the seasonal part, which is the same set of values every 52 values.
#1st value corresponds to 1,2,3..., 2nd value corresponds to 1+1/52,2+1/52,etc
plot(t1.1,xlim=c(1,6)) #adjust x axis to make room for predictions
lines(x=tp,y=sp,col='gray')

dev.off()

qmod = lm(t1.1 ~ z + I(z^2))
summary(qmod)
plot(t1.1)
points(x=z, y= qmod$fitted.values, type = 'l', col = 'red')
q.resl = qmod$residuals
plot(q.resl,type = 'l')
par(mfrow=c(1,2))
acf(q.resl,lag.max=100) 
pacf(q.resl,lag.max=100)


MSE2.full=numeric(0)
for(k in 1:26){ 
  har=harmonic(t1.1,k)
  
  qsmod = lm(q.resl ~ har)
  plot(q.resl,type='l', main = paste0("k value ",k))
  lines(qmod$fitted.values,col='red')
  acf(qsmod$residuals,lag.max=100) 
  #temp=readline(prompt='press key to continue') #allows pause before new graph is displayed
  predict = qmod$fitted.values + qsmod$fitted.values
  MSE2.full[k] = mean((predict - t1.1)^2)
}

##MSE on the full data set might result in overfitting.Let's find the appropriate K

#################################################################################
#Finding Appropriate K using the Cross-validation(CV1)
#Using l.computeCVmse to find out the best k in linear Model
#Using q.computeCVmse to find out the best k in quadratic Model
#So set training equal to 140, and test equal to 69 gives k=20.


k.l=matrix(0, nrow = 2, ncol = 26)
k.q=matrix(0, nrow = 2, ncol = 26)
for(i in 1:26)
{
  k.l[,i] = l.computeCVmse(i)
  k.q[,i] = q.computeCVmse(i)
}
which(k.l[2,] == min(k.l[2,]))
which(k.q[2,] == min(k.q[2,]))
##Both model choose k = 20.
##Quadratic fits better according to CV when k.level = 20.

##################################################################
##Finding Appropriate K using the Cross-validation(CV2)

k.l.2=matrix(0, nrow = 101, ncol = 26)
k.q.2=matrix(0, nrow = 101, ncol = 26)
for(i in 1:26)
{
  k.l.2[,i] = as.numeric(na.omit(l.computeCVmse2(i)))
  k.q.2[,i] = as.numeric(na.omit(q.computeCVmse2(i)))
}

best.k.l = numeric(0)
best.k.q = numeric(0)
for( i in 1:101)
{
  best.k.l[i] = which(k.l.2[i,] == min(k.l.2[i,]))
  best.k.q[i] = which(k.q.2[i,] == min(k.q.2[i,]))
}
table(best.k.l)
table(best.k.q)

mean.k.l = numeric(0)
mean.k.q = numeric(0)
for(i in 1:26)
{
  mean.k.l[i] = mean(k.l.2[,i])
  mean.k.q[i] = mean(k.q.2[,i])
}

which(mean.k.l == min(mean.k.l))
which(mean.k.q == min(mean.k.q))

mean.k.l
mean.k.q

l.computeCVmse = function(level.k){
  MSE = numeric(0)
  for(k in 2:1){
    train.dt = ts(t1.1[1:(70*k)], start = 1, frequency=52)
    start1.test = ceiling(70*k/52)
    start2.test = 70*k - floor(70*k/52)*52 +1
    test.dt = ts(t1.1[(70*k+1):length(t1.1)]
                 , start = c(start1.test,start2.test), frequency=52)
    z = time(t1.1)
    
    z.train.dt = z[1:(70*k)]
    z.test.dt = z[(70*k+1):length(t1.1)]
 
    t.lmod=lm(train.dt ~ z.train.dt)
    t.resl = t.lmod$residuals
      
    har=harmonic(train.dt,m=level.k)
    smod=lm(t.lmod$residuals~har)
    ssmod = rep(smod$fitted.values,3)
    
    s.value=numeric()
    flag = start2.test
    flag2 = 2
    for(i in length(train.dt)+1:length(test.dt))
    {
      s.value[flag2] = ssmod[flag]
      flag = flag+1
      flag2 = flag2+1
    }
    s.value = as.numeric(na.omit(s.value))[-length(s.value)]
      
    sp=t.lmod$coefficients[1]+t.lmod$coefficients[2]*z.test.dt + s.value
    
    MSE[k] = mean((test.dt - sp)^2)
  }
  return(na.omit(MSE))
}

l.computeCVmse2 = function(level.k){
  MSE = numeric(0)
  for(k in 108:208){
    train.dt = ts(t1.1[1:k], start = 1, frequency=52)
    start1.test = ceiling(k/52)
    start2.test = k - floor(k/52)*52 +1
    test.dt = ts(t1.1[(k+1):209]
                 , start = c(start1.test,start2.test), frequency=52)
    z = time(t1.1)
    
    z.train.dt = z[1:k]
    z.test.dt = z[(k+1):209]
    
    t.lmod=lm(train.dt ~ z.train.dt)
    t.resl = t.lmod$residuals
    
    har=harmonic(train.dt,m=level.k)
    smod=lm(t.lmod$residuals~har)
    ssmod = rep(smod$fitted.values,3)
    
    s.value=numeric()
    flag = start2.test
    flag2 = 2
    for(i in 1:length(z.test.dt))
    {
      s.value[flag2] = ssmod[flag]
      flag = flag+1
      flag2 = flag2+1
    }
    s.value = as.numeric(na.omit(s.value))[-length(s.value)]
    
    sp=t.lmod$coefficients[1]+t.lmod$coefficients[2]*z.test.dt + s.value
    
    MSE[k] = mean((test.dt - sp)^2)
  }
  return(MSE)
}


q.computeCVmse = function(level.k){
  MSE = numeric(0)
  for(k in 2:1){
    train.dt = ts(t1.1[1:(70*k)], start = 1, frequency=52)
    start1.test = ceiling(70*k/52)
    start2.test = 70*k - floor(70*k/52)*52 +1
    test.dt = ts(t1.1[(70*k+1):length(t1.1)]
                 , start = c(start1.test,start2.test), frequency=52)
    z = time(t1.1)
    
    z.train.dt = z[1:(70*k)]
    z.test.dt = z[(70*k+1):length(t1.1)]
    
    t.qmod=lm(train.dt ~ z.train.dt + I(z.train.dt^2))
    t.qresl = t.qmod$residuals
    
    har=harmonic(train.dt,m=level.k)
    smod=lm(t.qmod$residuals~har)
    ssmod = rep(smod$fitted.values,3)
    
    s.value=numeric()
    flag = start2.test
    flag2 = 2
    for(i in length(train.dt)+1:length(test.dt))
    {
      s.value[flag2] = ssmod[flag]
      flag = flag+1
      flag2 = flag2+1
    }
    s.value = as.numeric(na.omit(s.value))[-length(s.value)]
    
    sp=t.qmod$coefficients[1]+t.qmod$coefficients[2]*z.test.dt + t.qmod$coefficients[3]*(z.test.dt)^2+ s.value
    
    MSE[k] = mean((test.dt - sp)^2)
  }
  return(MSE)
}

q.computeCVmse2 = function(level.k){
  MSE = numeric(0)
  for(k in 108:208){
    train.dt = ts(t1.1[1:k], start = 1, frequency=52)
    start1.test = ceiling(k/52)
    start2.test = k - floor(k/52)*52 +1
    test.dt = ts(t1.1[(k+1):209]
                 , start = c(start1.test,start2.test), frequency=52)
    z = time(t1.1)
    
    z.train.dt = z[1:k]
    z.test.dt = z[(k+1):209]
    
    t.qmod=lm(train.dt ~ z.train.dt + I(z.train.dt^2))
    t.qresl = t.qmod$residuals
    
    har=harmonic(train.dt,m=level.k)
    smod=lm(t.qmod$residuals~har)
    ssmod = rep(smod$fitted.values,3)
    
    s.value=numeric()
    flag = start2.test
    flag2 = 2
    for(i in 1:length(z.test.dt))
    {
      s.value[flag2] = ssmod[flag]
      flag = flag+1
      flag2 = flag2+1
    }
    s.value = as.numeric(na.omit(s.value))[-length(s.value)]
    
    sp=t.qmod$coefficients[1]+t.qmod$coefficients[2]*z.test.dt + t.qmod$coefficients[3]*(z.test.dt)^2+ s.value
    
    MSE[k] = mean((test.dt - sp)^2)
  }
  return(MSE)
}
########################################################################
#######################################################################
#######################################################################
#K=20 FOR LINEAR MODEL
#K=5 FOR QUADRATIC MODEL
best.l = numeric(0)
best.q = numeric(0)

best.l = as.numeric(na.omit(l.computeCVmse2(20)))
tail(best.l)
mean(best.l)
best.q = as.numeric(na.omit(q.computeCVmse2(5)))
tail(best.q)
mean(best.q)
mean(as.numeric(na.omit(l.computeCVmse2(12))))

mean(l.computeCVmse(20))
#########################################################3
#########################################################
#Best Model : Linear trend with k=20


##########################################################################
##########################################################################
##Below code shows that when using the parametric trend and training set is too small
##the prediction will far off the result.
k=1
level.k=20
train.dt = ts(t1.1[1:(70*k)], start = 1, frequency=52)
start1.test = ceiling(70*k/52)
start2.test = 70*k - floor(70*k/52)*52 +1
test.dt = ts(t1.1[(70*k+1):length(t1.1)]
             , start = c(start1.test,start2.test), frequency=52)
z = time(t1.1)
plot(train.dt, xlim = c(1,5))
z.train.dt = z[1:(70*k)]
z.test.dt = z[(70*k+1):length(t1.1)]
t.qmod=lm(train.dt ~ 1+ z.train.dt + I(z.train.dt^2))
t.qresl = t.qmod$residuals
har=harmonic(train.dt,m=level.k)
smod=lm(t.qmod$residuals~har)
ssmod = rep(smod$fitted.values,3)
s.value=numeric()
flag = start2.test
flag2 = 2
for(i in length(train.dt)+1:length(test.dt))
{
  s.value[flag2] = ssmod[flag]
  flag = flag+1
  flag2 = flag2+1
}
s.value = as.numeric(na.omit(s.value))[-length(s.value)]
sp=t.qmod$coefficients[1]+t.qmod$coefficients[2]*(z.test.dt) + t.qmod$coefficients[3]*((z.test.dt)^2)+ s.value
plot(train.dt, xlim = c(1,5), ylim = c(40,200))
points(t.qmod$fitted.values, x = z.train.dt)
lines(y = test.dt, x= z.test.dt)
points(y = sp, x = z.test.dt, col = "red", ty = 'l')
##########################################################################
##########################################################################


