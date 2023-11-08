library(lubridate)
library(dplyr)
rm(list=ls());

#%slicing 1973:Q1 - 2021:Q1 (251, io uso 247 per evitare covid)
Q1_1971=ymd("1971/03/01");
#%Q1_2021=251;
Q1_2019=ymd("2019/03/01");


#-------------------------------------------------------------------------------
#                                 REAL DATA
#-------------------------------------------------------------------------------

Rdata0=read.csv("2023-09.csv")
Rdata0$sasdate=mdy(Rdata0$sasdate)
Rdata0=Rdata0[!is.na(Rdata0$sasdate),]

Rdata=Rdata0;
Rdata = Rdata %>% filter(sasdate >= Q1_1971 & sasdate <= Q1_2019)


N0=read.csv("CLF16OV.csv")
N0$DATE=ymd(N0$DATE)
N= N0 %>% filter( DATE>= Q1_1971 & DATE <= Q1_2019)
N= N %>% filter( month(DATE)==3 |month(DATE)==6 |month(DATE)==9 |month(DATE)==12 )

#GDP per Capita
Y=log((Rdata$GDPC1/N$CLF16OV)*10e6);
Ydiff1=log(Rdata0$GDPC1[Rdata0$sasdate==ymd("1971/03/01")]*10e6/N0$CLF16OV[N0$DATE==ymd("1971/03/01")])-log(Rdata0$GDPC1[Rdata0$sasdate==ymd("1970/12/01")]*10e6/N0$CLF16OV[N0$DATE==ymd("1970/12/01")]);
Ydiff2=diff(Y);
Ydiff=c(Ydiff1,Ydiff2)
#plot(Rdata$sasdate,Y, type="l")

#Inflation
P1=Rdata0$CPIAUCSL[Rdata0$sasdate==ymd("1971/03/01")]-Rdata0$CPIAUCSL[Rdata0$sasdate==ymd("1970/12/01")];
P2=diff(Rdata$CPIAUCSL)
P=c(P1,P2)
rm(P1,P2);
#plot(Rdata$sasdate,P, type="l")

#Nominal Interest Rate
R=(1+Rdata$FEDFUNDS/100)^(1/4)
#plot(Rdata$sasdate,R,type="l")


#-------------------------------------------------------------------------------
#                                 SIMULATED DATA
#-------------------------------------------------------------------------------
Sdata=read.csv("../Matlab_code/Simul_data_NK_NL_DSGE.csv")

SY=Sdata$log_y
SP=Sdata$Pi
SR=Sdata$R

start_date=Rdata$sasdate[1]
Sdata$date=seq(start_date, by = "quarter", length.out = length(Sdata$Y));

#Inserting the trend in Y
#lm=lm(Y~Rdata$sasdate)
#trend=lm$coefficients[2]
#rm(lm)
#SA=Sdata$A
#mu_z=vector("numeric",length = length(SY));
#mu_z[1]=SA[1]#a+min(Y);
#for (i in 2:length(SY)){
#  mu_z[i]=mu_z[i-1]+log(trend)+SA[i]
#}
#SY=SY+mu_z
#plot(Rdata$sasdate,mu_z, type="l")
#plot(Rdata$sasdate,SA, type="l")

par(mfrow = c(2, 3))
plot(Rdata$sasdate,Ydiff, type="l", main="Real GDP per Capita" , xlab="Quarter")
plot(Rdata$sasdate,P, type="l", main="Real Inflation", xlab="Quarter")
plot(Rdata$sasdate,R,type="l", main="Real Nominal Interest Rate",xlab="Quarter")

plot(Sdata$date,SY, type="l", main="Simulated GDP per Capita" , xlab="Quarter")
plot(Sdata$date,SP, type="l", main="Simulated Inflation", xlab="Quarter")
plot(Sdata$date,SR,type="l", main="Simulated Nominal Interest Rate",xlab="Quarter")





