library(lubridate)
library(dplyr)
rm(list=ls());

#%slicing 1973:Q1 - 2021:Q1 (251, io uso 247 per evitare covid)
Q1_1971=ymd("1971/03/01");
#%Q1_2021=251;
Q1_2019=ymd("2019/03/01");


#%'PCECC96'-3,'GDPC1'-2,'CPIAUCSL'-121,'AWHMAN'-78,'AHETPIx'-132

Rdata0=read.csv("2023-09.csv")
Rdata0$sasdate=mdy(Rdata0$sasdate)
Rdata0=Rdata0[!is.na(Rdata0$sasdate),]

Rdata=Rdata0;
Rdata = Rdata %>% filter(sasdate >= Q1_1971 & sasdate <= Q1_2019)
  

N=read.csv("CLF16OV.csv")
N$DATE=ymd(N$DATE)
N= N %>% filter( DATE>= Q1_1971 & DATE <= Q1_2019)
N= N %>% filter( month(DATE)==3 |month(DATE)==6 |month(DATE)==9 |month(DATE)==12 )

#GDP per Capita
Y=(Rdata$GDPC1/N$CLF16OV)*10e6;
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





              