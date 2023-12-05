library(vars)
library(tseries)
library(tidyverse)
library(stargazer)
library(fastICA)
rm(list=ls());

#-------------------------------------REAL DATA-----------------------------------
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
Est = VAR(D, p=lag, type= "none");
Estlev = VAR(Dlev, p=lag, type= "none");
summary (Estlev)

#Get the error term
?Acoef
?Bcoef
B=Bcoef(Estlev);
B1=B[,1:5];
B2=B[,6:10];
B3=B[,11:15];

mat=matrix(c(Y,C,P,N,W), nrow=965/5, ncol=5)
U = matrix(0, nrow = 5, ncol = 190);

t=0
for (t in 4:193){
  U[1:5,t-3]= t(t(mat[t,]))-B1%*%t(t(mat[t-1,]))+B2%*%t(t(mat[t-2,]))+B3%*%t(t(mat[t-3,]));
}

#Identifying the contemporaneous reletionship matrix A0

library(fastICA)
U=t(U)
?fastICA
X= fastICA(U, 5, alg.typ = "parallel", fun = "logcosh", alpha = 1,
           method = "C", row.norm = FALSE, maxit = 200,
           tol = 0.0001, verbose = TRUE)
MixingMat=X$A

A0=solve(MixingMat);

#compute the SVAR Matrix

A1=A0%*%B1;
A2=A0%*%B2;
A3=A0%*%B3;

#-------------------------------------SIMULATED DATA----------------------------

#let replicate all the process for simulated data:
Sdata=read.csv("Simul_data.csv")
SY=ts(Sdata$Y, start=c(2000,1), frequency=4);
SC=ts(Sdata$C, start=c(2000,1), frequency=4);
SP=ts(Sdata$P, start=c(2000,1), frequency=4);
SN=ts(Sdata$N, start=c(2000,1), frequency=4);
SW=ts(Sdata$W, start=c(2000,1), frequency=4);

Stime=1:250;

SDlev=data.frame(SY=as.vector(SY),SP=as.vector(SP),SN=as.vector(SN),SW=as.vector(SW));

#Optimal length
print(VARselect(SDlev)$selection) #it give with AIC 2
Slag=2

#Estimate
SEstlev = VAR(SDlev, p=Slag, type= "none");
summary (SEstlev)

#Get the error term
SB=Bcoef(SEstlev);
SB1=SB[,1:4]
SB2=SB[,5:8]


Smat=matrix(c(SY,SP,SN,SW), nrow=1000/4, ncol=4)
SU = matrix(0, nrow = 4, ncol = 250);

t=0
for (t in Slag+1:248){
  SU[1:4,t-Slag]= t(t(Smat[t,]))-SB1%*%t(t(Smat[t-1,]))+SB2%*%t(t(Smat[t-2,]));
}

#Identifying the contemporaneous reletionship matrix A0

library(fastICA)
SU=t(SU)

SX= fastICA(SU, 5, alg.typ = "parallel", fun = "logcosh", alpha = 1,
           method = "C", row.norm = FALSE, maxit = 200,
           tol = 0.0001, verbose = TRUE)
SMixingMat=SX$A

SA0=solve(SMixingMat);

#getting the SVAR Matrix

SA1=SA0%*%SB1;
SA2=SA0%*%SB2;










