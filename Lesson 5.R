# /Users/Tanner/Documents/KUMC/Linear Regressions
library(dplyr)

# Clean up
# rm(list=ls())

# 3.1
datasetA = read.table("/Users/Tanner/Documents/KUMC/Linear Regressions/APPENC05.txt")
dataset1 = datasetA %>% rename(X=V3, Y=V2)

# V3 must be square feet X, and V2 must be price as Y
c1 = dataset1
# attach(cancer1)
length(c1$Y)
cancer_model = lm(Y ~ X,c1)
# predict(re_model)

# DIAGNOSTICS
summary(cancer_model)
anova(cancer_model)

plot(c1$Y ~ c1$X, xlab="Cancer Volume",ylab="PSA Level",pch=16,
                        main="Figure E: Relationship of Cancer Volume and PSA Level")
abline(cancer_model)

library(Hmisc)
rcorr(c1$X,c1$Y)

# plot summary count of buckets
c2 = c1 %>% mutate(x_b = ifelse(X <= 10,"<= 10",
                          ifelse(X <= 20, "11-20",
                          ifelse(X <= 30, "21-30",
                          ifelse(X <= 40, "31-40",
                          ifelse(X <= 50, "41-50",
                          "> 50")))))) %>%
                group_by(x_b) %>%
                summarise(count = n())
plot(as.factor(c2$x_b),c2$count,pch=2, xlab="Cancer Volume Buckets"
     ,ylab="Count",pch=16
     ,main="Figure B: Cancer Volume Distribution Plot")

# fail
# dotchart(X,groups=re2$x_b)
dev.new()
boxplot(c1$Y, main = "Figure C: PSA Levels Boxplot",xlab = "PSA" , horizontal = T,ylim=c(0,105))
boxplot(c1$X, main = "Figure D: Cancer Volume Boxplot",xlab = "Cancer Volume" , horizontal = T)



# RESIDUALS
Y_hat = predict(cancer_model,data.frame(X))
e = Y - Y_hat
SE.e = summary(cancer_model)$sigma #s.d. of error
e.star = (e-mean(e)) /  SE.e #hand calc semi-studentized

# create a new plot with large dimensions
dev.new()
plot(e.star ~ c1$X,xlim=c(0,50),ylim=c(-10,10),col="red"
     ,main="Figure E: Semi-Studentized Residual Plot",xlab="Cancer Volume")
abline(0,0,lty=2)
# points(e.star ~ re1$X,col="red")
plot(e ~ c1$X,xlim=c(0,50),ylim=c(-30,30),col="red")

# same as manual e.star
c_resid = residuals(cancer_model)
dev.new()
plot(c_resid~c1$X)


# LINEARITY AND CONSTANCY OF VARIANCE
dev.new()
plot(X,c_resid,main="linearity and constant variance",xlab="Cancer Volume"
     ,ylab="Model Residual")
abline(0,0,lty=2)

dev.new()
plot(X,abs(c_resid),main="Figure F: X by Absolute Value of Residual",
     xlab="Cancer Volume",ylab="Absolute Value Model Residual")
abline(0,0,lty=2)


# Breusch-Pagan Test
# Tests H0: sigma^2 is a function of X 
# versus H1: sigma^2 is not a function of X
require(car)
ncvTest(cancer_model)


# OUTLIERS
plot(e.star~X)
points(rstudent(cancer_model)~X,col="red")
points(rstandard(cancer_model)~X,col="blue")
plot((rstudent(cancer_model)-e.star)~X)
plot((rstandard(cancer_model)-e.star)~X,main="Figure G: Leave-one-out Studentization of Residuals",
     xlab="Cancer Volume",ylab="Studentized Residual")
abline(0,0,lty=2)
max(rstandard(cancer_model)-e.star) 
# rstandard good for outlier detection because theyre normed by
# distance from mean(x) to adjust for leverage points

max(rstudent(cancer_model)-e.star)
# rstudent uses leave-one-out method becuase the value of sigma^2
# is derived without using the ith datapoint

plot(rstudent(cancer_model)~X)
abline(0,0,lty=2)
# identify crashed R
# identify(rstandard(re_model)~X)

# redmodel without outliers?
cancer_model2 = lm(Y~X, subset = rstudent(cancer_model)<2)
summary(cancer_model)
length(which(rstudent(cancer_model)<2))
summary(cancer_model2)

# INDEPENDENCE OF ERROR TERMS
plot(residuals(cancer_model),type="l",ylab=expression(e[i]),
     main="Figure I: Sequence Plot of Residuals ")
points(residuals(cancer_model2),pch=16,col="darkgray")
abline(0,0,lty=2)

# not working
# summary(lm(residuals(cancer_model)[-1]~-1+residuals(cancer_model2)[-200]))

# NORMALITY OF ERROR TERMS
qqnorm(residuals(cancer_model),main="Figure J: Normal Probability Plot of Residuals")
qqline(residuals(cancer_model))

hist(residuals(cancer_model))
boxplot(residuals(cancer_model),horizontal = T)

# Linear Relationship?
summary(cancer_model)
plot(Y~X)
abline(cancer_model)

# Confidence Band on Line
newx = seq(min(X),max(X),by=0.05)
plot(Y~X,pch=16,main="Figure K: Cancer Volume by PSA Level",
     xlab="Cancer Volume",ylab="PSA Level")
abline(cancer_model)


W <- sqrt(2*qf(0.95,2,length(X)-2))
SSX <- sum((X - mean(X))^2)
ll <- fitted(cancer_model) - W*(summary(cancer_model)$sigma)*sqrt((1/length(X)+(X-mean(X))^2/SSX))
ul <- fitted(cancer_model) + W*(summary(cancer_model)$sigma)*sqrt((1/length(X)+(X-mean(X))^2/SSX))
lines(X, ll, col="darkgray",lwd=2)
lines(X, ul, col="darkgray",lwd=2)


conf_interval <- predict(cancer_model, newdata=data.frame(X=newx), 
                         interval="confidence",
                         level = 0.95)
lines(newx, conf_interval[,2], col="darkgray", lty=2, lwd=2)
lines(newx, conf_interval[,3], col="darkgray", lty=2, lwd=2)

# Coefficient
coef(cancer_model)


pred_set_B = predict(cancer_model,data.frame(X = c(20),interval="predict"))

