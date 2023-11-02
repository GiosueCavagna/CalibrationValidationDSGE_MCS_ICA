library(vars)
library(tseries)
library(tidyverse)
library(stargazer)
library(fastICA)
rm(list=ls());

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
summary (Estlev)

#Get the error term
B=Bcoef(Estlev);
B1=B[,1:5]
B2=B[,6:10]
B3=B[,11:15]

mat=matrix(c(Y,C,P,N,W), nrow=965/5, ncol=5)
U = matrix(0, nrow = 5, ncol = 190);

t=0
for (t in 4:193){
  U[1:5,t-3]= t(t(mat[t,]))-B1%*%t(t(mat[t-1,]))+B2%*%t(t(mat[t-2,]))+B3%*%t(t(mat[t-3,]));
}

#Identifying the contemporaneous reletionship matrix A0

library(fastICA)
U=t(U)

X= fastICA(U, 5, alg.typ = "parallel", fun = "logcosh", alpha = 1,
           method = "C", row.norm = FALSE, maxit = 200,
           tol = 0.0001, verbose = TRUE)
MixingMat=X$A

A0=solve(MixingMat);

#getting the SVAR Matrix

A1=A0%*%B1;
A2=A0%*%B2;
A3=A0%*%B3;











