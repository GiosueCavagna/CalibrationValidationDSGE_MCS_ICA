library(vars)
library(tseries)
library(tidyverse)
library(stargazer)
library(fastICA)
rm(list=ls())

#%slicing 1973:Q1 - 2021:Q1 (251, io uso 247 per evitare covid)
Q1_1971=59;
#%Q1_2021=251;
Q1_2019=251;


#%'PCECC96'-3,'GDPC1'-2,'CPIAUCSL'-121,'AWHMAN'-78,'AHETPIx'-132


data=read.csv("2023-09.csv")

time=data$sasdate[Q1_1971:Q1_2019]
  
Y=ts(data$GDPC1[Q1_1971:Q1_2019], start=c(1971,1), frequency=4);
C=ts(data$PCECC96[Q1_1971:Q1_2019], start=c(1971,1), frequency=4);
P=ts(data$CPIAUCSL[Q1_1971:Q1_2019], start=c(1971,1), frequency=4);
N=ts(data$AWHMAN[Q1_1971:Q1_2019], start=c(1971,1), frequency=4);
W=ts(data$AHETPIx[Q1_1971:Q1_2019], start=c(1971,1), frequency=4);

#ADF test in order to check stationary of the series
print(adf.test(Y))
Yd=diff(Y)# since p-value>0.10 I reject null=> I differentiate
print(adf.test(C))
Cd=diff(C)
print(adf.test(P))
Pd=diff(P)
print(adf.test(N))# since p-value>0.10 I accept the null=> I just take one leg less
Nd=diff(N)
#Nd=ts(N[2:185],start=c(1971,1), frequency=4);
print(adf.test(W))
Wd=diff(W)

#Time=time[2:185]
D=data.frame(Yd=as.vector(Yd),Cd=as.vector(Cd),Pd=as.vector(Pd),Nd=as.vector(Nd),Wd=as.vector(Wd));
Dlev=data.frame(Y=as.vector(Y),C=as.vector(C),P=as.vector(P),N=as.vector(N),W=as.vector(W));

#Optimal length
print(VARselect(D)$selection) #it give with AIC 3 
lag=3
print(VARselect(Dlev)$selection) #it give with AIC 2
#laglev
#Estimate
?VAR
Est = VAR(D, p=lag, type= "none");
Estlev = VAR(Dlev, p=lag, type= "none");
summary(Est)

u=Yd







