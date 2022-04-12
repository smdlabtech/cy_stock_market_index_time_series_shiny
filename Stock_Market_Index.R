#------------------------------------------
# Install packages with the commands below
#------------------------------------------
install.packages('PerformanceAnalytics')
install.packages('quantmod')
install.packages('MASS')
install.packages('mvtnorm')
install.packages('mnormt')
install.packages('rugarch')
install.packages('car')
install.packages('FinTS')
install.packages('zoo')
install.packages('carData')
install.packages('xts')
install.packages('lubridate')

#------------------------------------------
# Loading packages
#------------------------------------------
library(PerformanceAnalytics)
library(quantmod)
library(MASS)
library(mvtnorm)
library(mnormt) 
library(rugarch)
library(carData)
library(car)
library(FinTS)
library(xts)
library(lubridate)

options(digits=4)

# download data 
symbol.vec = c("AMZN","FB")
getSymbols(symbol.vec, from ="2012-06-01", to = "2021-12-01")
colnames(AMZN)
start(AMZN)
end(AMZN)

AMZN = AMZN[,"AMZN.Adjusted", drop=F]
FB = FB[,"FB.Adjusted", drop=F]

# plot prices
dataToPlot1 = cbind(AMZN,FB)
colnames(dataToPlot1) = c("AMZN","FB") 
plot.zoo(dataToPlot1, main="Changes in financial rates", col = c("blue", "red"))

# calculate log-returns for GARCH analysis
AMZN.ret = CalculateReturns(AMZN, method="log")
FB.ret = CalculateReturns(FB, method="log")

# remove first NA observation
AMZN.ret = AMZN.ret[-1,]
colnames(AMZN.ret) ="AMZN"
FB.ret = FB.ret[-1,]
colnames(FB.ret) ="FB"

# create combined data series
data.ret = cbind(AMZN.ret,FB.ret)

# plot returns
colnames(data.ret) = c("AMZN.ret","FB.ret") 
plot.zoo(data.ret, main="Evolution des cours", col = c("blue", "red"))

#plot returns^2
dataToPlot = cbind(AMZN.ret^2,FB.ret^2)
colnames(dataToPlot) = c("AMZN.ret^2","FB.ret^2")
plot.zoo(dataToPlot, main = "Evolution of square yields", col = c("blue", "red"))

dataToPlot = cbind(abs(AMZN.ret), abs(FB.ret))
colnames(dataToPlot) = c("abs(AMZN.ret)","abs(FB.ret)")
plot.zoo(dataToPlot, main = "Evolution of absolute yields", col = c("blue", "red"))

# plot autocorrelations of returns, returns^2 and abs(returns)
par(mfrow=c(3,2))
acf(AMZN.ret, main="AMZN Returns")
acf(AMZN.ret, main="FB Returns")

acf(AMZN.ret^2, main="AMZN Returns^2")
acf(FB.ret^2, main="FB Returns^2")

acf(abs(AMZN.ret), main="AMZN abs(Returns)")
acf(abs(FB.ret), main="FB abs(Returns)")

### Histogramme ###
par(mfrow=c(2,1))
hist(AMZN.ret)
hist(FB.ret,col = c("orange"))

# compute summary statistics
A<-table.Stats(AMZN.ret)
B<-table.Stats(FB.ret)
tab_desc<-cbind(A,B)

dat=cbind(AMZN.ret,FB.ret)
summary(dat)

#------------------------------------------
#              Statistics tests 
#------------------------------------------
# use Box.test from stats package
Box.test(coredata(AMZN.ret^2), type="Ljung-Box", lag = 12)
Box.test(coredata(FB.ret^2), type="Ljung-Box", lag = 12)

# use ArchTest() function from FinTS package for Engle's LM test
ArchTest(AMZN.ret)
ArchTest(FB.ret)
ArchTest(dataToPlot)

#-----------------
# Asymmetric arch
arch11.spec = ugarchspec(variance.model = list(garchOrder=c(1,1)), 
                         mean.model = list(armaOrder=c(0,0)))
AMAZN.arch11.fit = ugarchfit(spec=arch11.spec, data=AMZN.ret, 
                             solver.control=list(trace = 1))
FB.arch11.fit = ugarchfit(spec=arch11.spec, data=FB.ret, 
                          solver.control=list(trace = 1))

AMAZN.arch11.fit
FB.arch11.fit

#---------------------
#   Asymmetric garch
#---------------------
garch11.spec = ugarchspec(variance.model = list(garchOrder=c(1,1)),
                          mean.model = list(armaOrder=c(0,0)))
AMZN.garch11.fit = ugarchfit(spec=garch11.spec, data=AMZN.ret, 
                             solver.control=list(trace = 1)) 
FB.garch11.fit = ugarchfit(spec=garch11.spec, data=FB.ret, 
                           solver.control=list(trace = 1)) 

AMZN.garch11.fit
FB.garch11.fit

# estimated coefficients
AMZN_ret<-coef(AMZN.garch11.fit)
FB_ret<-coef(FB.garch11.fit)
tab2<-rbind(AMZN_ret,FB_ret)

# unconditional mean in mean equation
uncmean(FB.garch11.fit)
uncmean(AMZN.garch11.fit)

# unconditional variance: omega/(alpha1 + beta1)
uncvariance(FB.garch11.fit)
uncvariance(AMZN.garch11.fit)

# persistence: alpha1 + beta1
persistence(AMZN.garch11.fit)
persistence(FB.garch11.fit)

# half-life:
halflife(AMZN.garch11.fit)
halflife(FB.garch11.fit)

# residuals: e(t)
plot.ts(residuals(AMZN.garch11.fit), ylab="e(t)", col="blue")
abline(h=0)
plot.ts(residuals(FB.garch11.fit), ylab="e(t)", col="blue")
abline(h=0)

# sigma(t) = conditional volatility
plot.ts(sigma(AMZN.garch11.fit), ylab="sigma(t)", col="blue")
plot.ts(sigma(FB.garch11.fit), ylab="sigma(t)", col="red")

#----------------------------------
# illustrate plot method
plot(AMZN.garch11.fit, which=1)
plot(AMZN.garch11.fit, which="all")
plot(AMZN.garch11.fit, which=9, main="AMZN")
plot(FB.garch11.fit, which=9, main="FB")

# Engle-Ng sign bias test
signbias(AMZN.garch11.fit)
signbias(FB.garch11.fit)

#--------------
# Egarch model
#--------------
egarch11.spec = ugarchspec(variance.model=list(model="eGARCH", 
                                               garchOrder=c(1,1)),
                           mean.model=list(armaOrder=c(0,0)))
AMZN.egarch11.fit = ugarchfit(egarch11.spec, AMZN.ret)
FB.egarch11.fit = ugarchfit(egarch11.spec, FB.ret)
AMZN.egarch11.fit
FB.egarch11.fit

#----------------
# GJR garch model
#----------------
gjrgarch11.spec = ugarchspec(variance.model=list(model="gjrGARCH",
                                                 garchOrder=c(1,1)),
                             mean.model=list(armaOrder=c(0,0)))
AMZN.gjrgarch11.fit = ugarchfit(gjrgarch11.spec, AMZN.ret)
FB.gjrgarch11.fit = ugarchfit(gjrgarch11.spec, FB.ret)
AMZN.gjrgarch11.fit
FB.gjrgarch11.fit

#--------------
# Aparch model
#--------------
aparch11.1.spec = ugarchspec(variance.model=list(model="apARCH",
                                                 garchOrder=c(1,1)), 
                             mean.model=list(armaOrder=c(0,0)), 
                             fixed.pars=list(delta=1))
AMZN.aparch11.1.fit = ugarchfit(aparch11.1.spec, AMZN.ret)
FB.aparch11.1.fit = ugarchfit(aparch11.1.spec, FB.ret)
AMZN.aparch11.1.fit
FB.aparch11.1.fit

#----------------------
# Information criteria
#----------------------
# AMZN
nic.garch11 = newsimpact(AMZN.garch11.fit)
nic.egarch11 = newsimpact(AMZN.egarch11.fit)
nic.gjrgarch11 = newsimpact(AMZN.gjrgarch11.fit)
nic.aparch11.1 = newsimpact(AMZN.aparch11.1.fit)

model.list = list(garch11 = AMZN.garch11.fit,
                  egarch11 = AMZN.egarch11.fit,
                  gjrgarch11 = AMZN.gjrgarch11.fit,
                  aparch11.1 = AMZN.aparch11.1.fit)
info.mat_amzn = sapply(model.list, infocriteria)
rownames(info.mat_amzn) = rownames(infocriteria(AMZN.garch11.fit))
info.mat_amzn<-t(info.mat_amzn)
info.mat_amzn<-abs(info.mat_amzn)

#x=c("Akaike","Bayes","Shibata","Hannan-Quinn")
x=c(1,2,3,4)
plot((info.mat_amzn[1,1:4])~x,xlab="Criteria index",
     ylab="criteria value",type='l',col="blue",
     main="Information criteria",xlim=c(1,4), ylim=c(5.195,5.25))
lines(info.mat_amzn[2,1:4],col="gray")
lines(info.mat_amzn[3,1:4],col="black")
lines(info.mat_amzn[4,1:4],col="orange")
legend("topleft", legend = c("garch(1,1)","egarch(1,1)",
                             "gjrgarch(1,1)","aparch(1,1)"),
       fill = c("blue","gray","black", "orange"))


# FB
nic.garch11 = newsimpact(FB.garch11.fit)
nic.egarch11 = newsimpact(FB.egarch11.fit)
nic.gjrgarch11 = newsimpact(FB.gjrgarch11.fit)
nic.aparch11.1 = newsimpact(FB.aparch11.1.fit)

model.list = list(garch11 = FB.garch11.fit,
                  egarch11 = FB.egarch11.fit,
                  gjrgarch11 = FB.gjrgarch11.fit,
                  aparch11.1 = FB.aparch11.1.fit)
info.mat_fb = sapply(model.list, infocriteria)
rownames(info.mat_fb) = rownames(infocriteria(FB.garch11.fit))
info.mat_fb<-t(info.mat_fb)
info.mat_fb<-abs(info.mat_fb)

#x=c("Akaike","Bayes","Shibata","Hannan-Quinn")
x=c(1,2,3,4)
plot((info.mat_fb[1,1:4])~x,xlab="Criteria index",
     ylab="criteria value",type='l',col="blue",
     main="Information criteria",xlim=c(1,4), ylim=c(4.83,4.9))
lines(info.mat_fb[2,1:4],col="gray")
lines(info.mat_fb[3,1:4],col="black")
lines(info.mat_fb[4,1:4],col="orange")
legend("topleft", legend = c("garch(1,1)","egarch(1,1)",
                             "gjrgarch(1,1)","aparch(1,1)"),
       fill = c("blue","gray","black", "orange"))

#****************************************************************
#       Graphic representation of the model Garch(1,1) 
# recall normal GARCH(1,1), examine standardized residuals
#AMZN
variance_predite <- AMZN.garch11.fit@fit$var
variance_observee <- (AMZN.garch11.fit@fit$residuals)^2
plot(variance_observee, main="pred variance vs actual variance")
lines(variance_predite, col='blue')

#FB
variance_predite <- FB.garch11.fit@fit$var
variance_observee <- (FB.garch11.fit@fit$residuals)^2
plot(variance_observee, main="pred variance vs actual variance")
lines(variance_predite, col='blue')

# show news impact curve from estimated garch(1,1) and egarch(1,1)
par(mfrow=c(2,2))
plot(nic.garch11$zx, type="l", lwd=2, col="blue", main="GARCH(1,1)", 
     nic.garch11$zy, ylab=nic.garch11$yexpr, xlab=nic.garch11$xexpr)
plot(nic.egarch11$zx, type="l", lwd=2, col="blue", main="EGARCH(1,1)",
     nic.egarch11$zy, ylab=nic.egarch11$yexpr, xlab=nic.egarch11$xexpr)
plot(nic.gjrgarch11$zx, type="l", lwd=2, col="blue", main="TGARCH(1,1)",
     nic.gjrgarch11$zy, ylab=nic.gjrgarch11$yexpr, xlab=nic.gjrgarch11$xexpr)
plot(nic.aparch11.1$zx, type="l", lwd=2, col="blue", main="APARCH(1,1,1)", 
     nic.aparch11.1$zy, ylab=nic.aparch11.1$yexpr, xlab=nic.aparch11.1$xexpr)

#**********************************************************************************
#         Estimation of asymmetric models (ApARCH and EGARCH)
#**********************************************************************************
##### garch with non-normal errors
par(mfrow=c(1,1))
# Aparch
AMZN.aparch11.1.zt = residuals(AMZN.aparch11.1.fit)/sigma(AMZN.aparch11.1.fit)
AMZN.aparch11.1.zt = xts(AMZN.aparch11.1.zt, order.by=index(AMZN.ret))
plot(coredata(AMZN.aparch11.1.zt))
par(mfrow=c(1,2))
plot(AMZN.aparch11.1.fit, which=8)
plot(AMZN.aparch11.1.fit, which=9)

FB.aparch11.1.zt = residuals(FB.aparch11.1.fit)/sigma(FB.aparch11.1.fit)
FB.aparch11.1.zt = xts(FB.aparch11.1.zt, order.by=index(FB.ret))
plot(coredata(FB.aparch11.1.zt))
par(mfrow=c(1,2))
plot(FB.aparch11.1.fit, which=8)
plot(FB.aparch11.1.fit, which=9)


# Egarch
par(mfrow=c(1,1))
AMZN.egarch11.zt = residuals(AMZN.egarch11.fit)/sigma(AMZN.egarch11.fit)
AMZN.egarch11.zt = xts(AMZN.egarch11.zt, order.by=index(AMZN.ret))
plot(coredata(AMZN.egarch11.zt))
par(mfrow=c(1,2))
plot(AMZN.egarch11.fit, which=8)
plot(AMZN.egarch11.fit, which=9)

par(mfrow=c(1,1))
FB.egarch11.zt = residuals(FB.egarch11.fit)/sigma(FB.egarch11.fit)
FB.egarch11.zt = xts(FB.egarch11.zt, order.by=index(FB.ret))
plot(coredata(FB.egarch11.zt))
par(mfrow=c(1,2))
plot(FB.egarch11.fit, which=8)
plot(FB.egarch11.fit, which=9)

##### with Student-t errors ####
# Aparch
aparch11.1.t.spec = ugarchspec(variance.model = list(model="apARCH",
                                                     garchOrder=c(1,1)), 
                               mean.model = list(armaOrder=c(0,0)),
                               distribution.model = "std", 
                               fixed.pars=list(delta=1))

par(mfrow=c(1,1))
AMZN.aparch11.1.t.fit = ugarchfit(spec=aparch11.1.t.spec, data=AMZN.ret)                             
AMZN.aparch11.1.t.fit

FB.aparch11.1.t.fit = ugarchfit(spec=aparch11.1.t.spec, data=FB.ret)                             
FB.aparch11.1.t.fit

par(mfrow=c(1,2))
plot(AMZN.aparch11.1.t.fit, which=3)
legend(x="topright", legend=c("AMZN"),lwd=2, lty = "solid")

plot(FB.aparch11.1.t.fit, which=3)
legend(x="topright", legend=c("FB"),lwd=2, lty = "solid")




# Egarch
egarch11.t.spec = ugarchspec(variance.model = list(model="eGARCH", 
                                                   garchOrder=c(1,1)),
                             mean.model = list(armaOrder=c(0,0)),
                             distribution.model = "std",
                             fixed.pars=list(delta=1))

AMZN.egarch11.t.fit = ugarchfit(spec=egarch11.t.spec, data=AMZN.ret)                             
AMZN.egarch11.t.fit

FB.egarch11.t.fit = ugarchfit(spec=egarch11.t.spec, data=FB.ret)                             
FB.egarch11.t.fit

par(mfrow=c(1,2))
plot(AMZN.egarch11.t.fit, which=3)
legend(x="topright", legend=c("AMZN"),lwd=2, lty = "solid")

plot(FB.egarch11.t.fit, which=3)
legend(x="topright", legend=c("FB"),lwd=2, lty = "solid")

#**************************************************************************
# fit skewed t :  
#**************************************************************************
# Aparch
aparch11.1.st.spec = ugarchspec(variance.model = list(model="apARCH",
                                                      garchOrder=c(1,1)), 
                                mean.model = list(armaOrder=c(0,0)),
                                distribution.model = "sstd",
                                fixed.pars=list(delta=1))


AMZN.aparch11.1.st.fit = ugarchfit(spec=aparch11.1.t.spec, data=AMZN.ret)                             
AMZN.aparch11.1.st.fit
# plot(AMZN.aparch11.1.st.fit, which=9) 

FB.aparch11.1.st.fit = ugarchfit(spec=aparch11.1.t.spec, data=FB.ret)                             
FB.aparch11.1.st.fit

par(mfrow=c(1,2))
plot(AMZN.aparch11.1.st.fit, which=9) 
legend(x="topleft", legend=c("fit skewed t AMZN"),lwd=2, lty = "solid")
plot(FB.aparch11.1.st.fit, which=9) 
legend(x="topleft", legend=c("fit skewed t FB"),lwd=2, lty = "solid")


#************************************************************************
## Plot forecasts from competing models : Forecasting for sigma (Monthly)
#************************************************************************
#***AMZN (mensuel) ***
par(mfrow=c(1,2))
AMZN.egarch11.fcst = ugarchforecast(AMZN.egarch11.fit, n.ahead=250)
AMZN.egarch11.t.fcst = ugarchforecast(AMZN.egarch11.t.fit, n.ahead=250)
plot(AMZN.egarch11.fcst, which = 3)   # Which = 1, graphic de la serie ; Which = 3 volatilite
legend(x="topright", legend=c("Egarch(1,1)"),lwd=2, lty = "solid")


AMZN.aparch11.1.fcst = ugarchforecast(AMZN.aparch11.1.fit, n.ahead=250)
AMZN.aparch11.1.t.fcst = ugarchforecast(AMZN.aparch11.1.t.fit, n.ahead=250)
plot(AMZN.aparch11.1.fcst, which = 3) 
legend(x="topright", legend=c("Aparch(1,1)"),lwd=2, lty = "solid")

#***FB***
par(mfrow=c(1,2))
FB.egarch11.fcst = ugarchforecast(FB.egarch11.fit, n.ahead=250)
FB.egarch11.t.fcst = ugarchforecast(FB.egarch11.t.fit, n.ahead=250)
plot(FB.egarch11.fcst, which = 1)
legend(x="topright", legend=c("Egarch(1,1)"),lwd=2, lty = "solid")

FB.aparch11.1.fcst = ugarchforecast(FB.aparch11.1.fit, n.ahead=250)
FB.aparch11.1.t.fcst = ugarchforecast(FB.aparch11.1.t.fit, n.ahead=250)
plot(FB.aparch11.1.fcst, which = 1)
legend(x="topright", legend=c("Aparch(1,1)"),lwd=2, lty = "solid")


#*****************************************************************
# extract volatility forecasts
#*****************************************************************
#***AMZN***
AMZN.egarch11.sigma = as.data.frame(AMZN.egarch11.fcst)$sigma
AMZN.egarch11.t.sigma = as.data.frame(AMZN.egarch11.t.fcst)$sigma
AMZN.aparch11.1.sigma = as.data.frame(AMZN.aparch11.1.fcst)$sigma
AMZN.aparch11.1.t.sigma = as.data.frame(AMZN.aparch11.1.t.fcst)$sigma

ymax = max(AMZN.egarch11.sigma,AMZN.egarch11.t.sigma,AMZN.aparch11.1.sigma, 
           AMZN.aparch11.1.t.sigma)
ymin = min(AMZN.egarch11.sigma,AMZN.egarch11.t.sigma,AMZN.aparch11.1.sigma, 
           AMZN.aparch11.1.t.sigma)

#***FB***
FB.egarch11.sigma = as.data.frame(FB.egarch11.fcst)$sigma
FB.egarch11.t.sigma = as.data.frame(FB.egarch11.t.fcst)$sigma
FB.aparch11.1.sigma = as.data.frame(FB.aparch11.1.fcst)$sigma
FB.aparch11.1.t.sigma = as.data.frame(FB.aparch11.1.t.fcst)$sigma

ymax = max(FB.egarch11.sigma,FB.egarch11.t.sigma,FB.aparch11.1.sigma, 
           FB.aparch11.1.t.sigma)
ymin = min(FB.egarch11.sigma,FB.egarch11.t.sigma,FB.aparch11.1.sigma, 
           FB.aparch11.1.t.sigma)


#------------
#   Plot ts
#------------
#***AMZN***
plot.ts(AMZN.egarch11.sigma, main="Volatility Forecasts",
        ylim=c(ymin,ymax), col="black", 
        lwd=2, ylab="sigma(t+h|t)", xlab="h")
lines(AMZN.egarch11.t.sigma, col="blue", lwd=2)
lines(AMZN.aparch11.1.sigma, col="green", lwd=2)
lines(AMZN.aparch11.1.t.sigma, col="red", lwd=2)
legend(x="topleft", legend=c("GARCH-n", "GARCH-t", "APARCH-n", "APARCH-t"),
       col=c("black", "blue","green","red"), lwd=2, lty = "solid")

#***FB***
plot.ts(FB.egarch11.sigma, main="Volatility Forecasts", ylim=c(ymin,ymax), 
        col="black", 
        lwd=2, ylab="sigma(t+h|t)", xlab="h")
lines(FB.egarch11.t.sigma, col="blue", lwd=2)
lines(FB.aparch11.1.sigma, col="green", lwd=2)
lines(FB.aparch11.1.t.sigma, col="red", lwd=2)
legend(x="topleft", legend=c("GARCH-n", "GARCH-t", "APARCH-n", "APARCH-t"),
       col=c("black", "blue","green","red"), lwd=2, lty = "solid")


#*** AMZN ***
AMZN.egarch11.fit = ugarchfit(spec=egarch11.spec, data=AMZN.ret, 
                              out.sample=100)
AMZN.egarch11.fit = ugarchfit(spec=egarch11.t.spec, data=AMZN.ret, 
                              out.sample=100)

AMZN.aparch11.1.fit = ugarchfit(aparch11.1.spec, AMZN.ret, 
                                out.sample=100)
AMZN.aparch11.1.t.fit = ugarchfit(spec=aparch11.1.t.spec, data=AMZN.ret, 
                                  out.sample=100)
plot(AMZN.aparch11.1.fit)


#*** FB ***
FB.egarch11.fit = ugarchfit(spec=egarch11.spec, data=FB.ret, out.sample=100)
FB.egarch11.fit = ugarchfit(spec=egarch11.t.spec, data=FB.ret, out.sample=100)

FB.aparch11.1.fit = ugarchfit(aparch11.1.spec, FB.ret, out.sample=100)
FB.aparch11.1.t.fit = ugarchfit(spec=aparch11.1.t.spec, data=FB.ret, 
                                out.sample=100)


#****************************************************
#   compare persistence and unconditional variance
#****************************************************
c.mat = matrix(0, 2, 2)
colnames(c.mat) = c("Persistence", "E[sig(t)]")
rownames(c.mat) = c("EGARCH-n", "EGARCH-t", "APARCH-n","APARCH-t")

c.mat["EGARCH-n","Persistence"] = persistence(FB.egarch11.fit)
c.mat["EGARCH-t","Persistence"] = persistence(FB.egarch11.t.fit)
c.mat["APARCH-n","Persistence"] = persistence(FB.aparch11.1.fit)
c.mat["APARCH-t","Persistence"] = persistence(FB.aparch11.1.t.fit)

c.mat["EGARCH-n","Persistence"] = persistence(AMZN.egarch11.fit)
c.mat["EGARCH-t","Persistence"] = persistence(AMZN.egarch11.t.fit)
c.mat["APARCH-n","Persistence"] = persistence(AMZN.aparch11.1.fit)
c.mat["APARCH-t","Persistence"] = persistence(AMZN.aparch11.1.t.fit)

c.mat["EGARCH-n","E[sig(t)]"] = sqrt(uncvariance(FB.egarch11.fit))
c.mat["EGARCH-t","E[sig(t)]"] = sqrt(uncvariance(FB.egarch11.t.fit))
c.mat["APARCH-n","E[sig(t)]"] = sqrt(uncvariance(FB.aparch11.1.fit))
c.mat["APARCH-t","E[sig(t)]"] = sqrt(uncvariance(FB.aparch11.1.fit))

c.mat["EGARCH-n","E[sig(t)]"] = sqrt(uncvariance(AMZN.egarch11.fit))
c.mat["EGARCH-t","E[sig(t)]"] = sqrt(uncvariance(AMZN.egarch11.t.fit))
c.mat["APARCH-n","E[sig(t)]"] = sqrt(uncvariance(AMZN.aparch11.1.fit))
c.mat["APARCH-t","E[sig(t)]"] = sqrt(uncvariance(AMZN.aparch11.1.fit))

c.mat

#*******************************************
# compute 100 1-step ahead rolling forecasts
AMZN.egarch11.fcst = ugarchforecast(AMZN.egarch11.fit, n.roll=100, 
                                    n.ahead=1)
AMZN.egarch11.t.fcst = ugarchforecast(AMZN.egarch11.t.fit, n.roll=100,
                                      n.ahead=1)
AMZN.aparch11.1.fcst = ugarchforecast(AMZN.aparch11.1.fit, n.roll=100, 
                                      n.ahead=1)
AMZN.aparch11.1.t.fcst = ugarchforecast(AMZN.aparch11.1.fit, n.roll=100, 
                                        n.ahead=1)

FB.egarch11.fcst = ugarchforecast(FB.egarch11.fit, n.roll=100, n.ahead=1)
FB.egarch11.t.fcst = ugarchforecast(FB.egarch11.t.fit, n.roll=100, n.ahead=1)
FB.aparch11.1.fcst = ugarchforecast(FB.aparch11.1.fit, n.roll=100, n.ahead=1)
FB.aparch11.1.t.fcst = ugarchforecast(FB.aparch11.1.fit, n.roll=100, n.ahead=1)

#*******************************************
# compute forecast evaluation statistics
fcst.list = list(egarch11=FB.egarch11.fcst,
                 egarch11.t=FB.egarch11.t.fcst,
                 aparch11.1=FB.aparch11.1.fcst,
                 aparch11.t.1=FB.aparch11.1.t.fcst)
fpm.mat = sapply(fcst.list, fpm)
fpm.mat

#*******************************************
# compute forecast evaluation statistics
fcst.list = list(egarch11=AMZN.egarch11.fcst,
                 egarch11.t=AMZN.egarch11.t.fcst,
                 aparch11.1=AMZN.aparch11.1.fcst,
                 aparch11.t.1=AMZN.aparch11.1.t.fcst)
fpm.mat = sapply(fcst.list, fpm)
fpm.mat







