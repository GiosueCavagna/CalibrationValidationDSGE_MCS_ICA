library(lubridate)
library(dplyr)
library(moments)
library(fastICA)
library(clue)
library(steadyICA)
library(data.table)
library(R.matlab)
library(devtools)
install_github("STieleman/svarCAL")
library(svarCAL)

rm(list=ls());


#REAL DATA-----
#-Data import and arrangement----

Rdata0=read.csv("2023-09.csv")
Rdata0$sasdate=mdy(Rdata0$sasdate)
Rdata0=Rdata0[!is.na(Rdata0$sasdate),]


#%slicing 1973:Q1 - 2019:Q1 (247)
Q1_1971=ymd("1971/03/01");
#%Q1_2019=247;
Q1_2019=ymd("2019/03/01");

Rdata=Rdata0;
Rdata = Rdata %>% filter(sasdate >= Q1_1971 & sasdate <= Q1_2019)


N0=read.csv("CLF16OV.csv")
N0$DATE=ymd(N0$DATE)
N= N0 %>% filter( DATE>= Q1_1971 & DATE <= Q1_2019)
N= N %>% filter( month(DATE)==3 |month(DATE)==6 |month(DATE)==9 |month(DATE)==12 )

#GDP per Capita
Y=log((Rdata$GDPC1/N$CLF16OV));
Ydiff1=log(Rdata0$GDPC1[Rdata0$sasdate==ymd("1971/03/01")]/N0$CLF16OV[N0$DATE==ymd("1971/03/01")])-log(Rdata0$GDPC1[Rdata0$sasdate==ymd("1970/12/01")]/N0$CLF16OV[N0$DATE==ymd("1970/12/01")]);
Ydiff2=diff(Y);
Ydiff=c(Ydiff1,Ydiff2)
rm(Y,Ydiff1,Ydiff2);
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

D=data.frame(Column1 = Ydiff, Column2 = P, Column3 = R)
colnames(D) <- c("Ydiff", "P", "R")

#-Data elaboration----

infocrit = VARselect(D, lag.max=6, type="const")
lag = infocrit$selection["AIC(n)"]
varest = VAR(D, p=lag, type="const")
ures = resid(varest)
AA = Acoef(varest) # VAR coefficients
k = ncol(ures)

# Test Gaussianity 
jarque_result = c(jarque.test(ures[,1])$p.value,jarque.test(ures[,2])$p.value,jarque.test(ures[,3])$p.value)
#we refure the null -> not normal at 5%
A_rw = fAp_fastICA(ures,sseed=46) # rw mixing matrix

#SIMULATED DATA----
source(paste("frobICA_mod.R", sep=""))
refM = A_rw

tau= 200 #200 # simulation length
nmr= 500 #500 # number of Monte Carlo runs
ncp= 500 #500 # number of configuration of parameters
k <- nrow(refM) # number of variables

#-Data import----
Simul_logY=readMat("Data/Simul_logY.mat") #readMat("~/Documents/CalibrationValidationDSGE_MCS_ICA/Matlab_code/Simulated_Data/Simul_logY.mat") #
Simul_Pi=readMat("Data/Simul_PI.mat") #readMat("~/Documents/CalibrationValidationDSGE_MCS_ICA/Matlab_code/Simulated_Data/Simul_Pi.mat") #
Simul_R=readMat("Data/Simul_R.mat") #readMat("~/Documents/CalibrationValidationDSGE_MCS_ICA/Matlab_code/Simulated_Data/Simul_R.mat") 
#[tau,nmr,ncp]

#-Creation of simulation dataset----
D = list()
i=1
Error=matrix(0,3,ncp);
for(ii in 1:ncp){
  Error[c(1,2),ii]=c(ii,"Rank Not or Div");
  #Check if the CoP return Nan simulation
  if (!(is.nan(Simul_logY$Simul.logY[1,1,ii])| is.nan(Simul_logY$Simul.logY[200,500,ii]) | abs(Simul_logY$Simul.logY[200,500,ii])>10e+2 | abs(Simul_logY$Simul.logY[200,500,ii])<10e-14)){ #Check if the CoP return divergent simulation 
    D[[i]] = list()
    D[[i]]$list = list()
    for (j in 1:nmr){
      D[[i]]$list[[j]] = as.data.frame(cbind(Simul_logY$Simul.logY[,j,ii],Simul_Pi$Simul.Pi[,j,ii],Simul_R$Simul.R[,j,ii]))
      colnames(D[[i]]$list[[j]]) = c('log_Y','Pi','R')
    }
    Error[,ii]=c(ii,"No",i);
    i=i+1
  }
  print(ii)
}
sum(Error[2,]!="No")
ncp_n=i-1
#-Estimation procedure----
Mres = as.list(1:ncp_n)
i=14
for(i in 1:ncp_n) {
  Mres[[i]] <- list()
  Mres[[i]]$m_ica <- as.list(1:nmr) # it delivers a mixing matrix for each monte carlo run
  Mres[[i]]$f_dist <- rep(NA,nmr)
  Mres[[i]]$ica_perm <- as.list(1:nmr)
  
  #estimate VARs
  m_VAR <- tryCatch(
           { lapply(D[[i]][[1]],VAR,p=la)},# mixing matrices of model data
            error = function(e) {
            m_VAR #In case of problem I put the previous 
           }
           )#lag) # var of model data
  m_res <- lapply(m_VAR,residuals) # residuals of model data
  #VAR coefficients
  AA <- lapply(m_VAR,coef)
  
  #estimate mixing matrix using fastICA
  Mres[[i]]$m_ica <-  tryCatch(
                      { lapply(m_res,fAp_fastICA,sseed=46)},# mixing matrices of model data
                      error = function(e) {
                        Mres[[i-1]]$m_ica #In case of problem I put the previous 
                      }
                      )
  mt_ica <- lapply(Mres[[i]]$m_ica,t)
  fb <- lapply(mt_ica,frobICA_mod,refM)
  perm <- list()
  for(j in 1:nmr) {
    perm[[j]] <- fb[[j]]$perm
    Mres[[i]]$f_dist[[j]] <- fb[[j]]$frob_dist
  }
  Mres[[i]]$ica_perm <- Map(function(x, y)x%*%y,Mres[[i]]$m_ica,perm)
  print(i)
}
counter=length(Mres)
Mres=unique(Mres) #I remove the problematic term
countICA=426-375#counter-length(Mres)

fb_dist <- matrix(0,ncp_n-countICA,nmr)
for (i in 1:(ncp_n-countICA)) {
  fb_dist[i,] <- Mres[[i]]$f_dist # MDI for each CoPs and MC runs
}

save(fb_dist,file="M_dist.Rdata") # change path
save(Mres,file="Mres.Rdata") # change path

#MCS----
M_dist=fb_dist

iMC =500# 200
iCoP =ncp_n-countICA# 200

#load ("~/Documents/CalibrationValidationDSGE_MCS_ICA/R_code/M_dist.Rdata")

vMean.val <- apply(M_dist,1,fnMean)
vVar.val <- apply(M_dist,1,fnVar)
mVal <- cbind(1:iCoP,vMean.val,vVar.val)
mVal <- mVal[order(mVal[,2],decreasing=T),]

vIndex <- 1:iCoP
dA <- 1
iN <- iMC

mPValue.val <- fnMCS(vMean.val,vVar.val,vIndex,dA,verbose=0)
# mPValue.val <- mPValue.val[order(mPValue.val[,2],increasing=T),]
mPValue.val <- cbind(mPValue.val,mVal[,2:3])
vCoP.pass <- mPValue.val[which(mPValue.val[,2]>0.05)]

setwd("...") # change path

MCS <- paste("MCS_calib.csv",sep=",")
write.csv(mPValue.val,MCS,col.names=T)

#Function----
fAp_fastICA<-function(ures, sseed=46){
  set.seed(sseed)
  n<-ncol(ures)
  X<-t(ures)
  #jb.info= c(jarque.test(ures[,1])$p.value,jarque.test(ures[,2])$p.value,jarque.test(ures[,3])$p.value)
  icares <- fastICA(t(X), nrow(X),tol=1e-14, maxit=3000, verbose=FALSE) #Doris: tol=1e-14
  set.seed(NULL)
  K=icares$K
  W <- t((icares$K) %*% (icares$W)) 
  A <-solve(W) # A is the mixing matrix
  # A <- ginv(W)
  eres<- t(W %*% X) # 
  aba<-abs(A)
  CC<-matrix(0,n,2)
  for(i in 1:n){
    cc<-which(aba == max(aba), arr.ind = TRUE) # coordinate of the matrix where is the max entry
    aba[cc[1],]<-0
    aba[,cc[2]]<-0
    CC[i,]<-cc
  }
  sn<-1:n
  for(i in 1:n){sn[i]<-CC[CC[,1]==i,2]}
  Ap<-A[,sn]
  cnames<-paste("s", 1:n, sep="")
  colnames(Ap) <-cnames
  #if(Ap[1,1]>0){Ap[,1]<- -Ap[,1]}
  for (i in 1:n){if(Ap[i,i]<0){Ap[,i]<- -Ap[,i]}}
  Ap
}
## Function mean

fnMean <- function(x) {
  mean(x,na.rm=T)
}

## Function variance

fnVar <- function(x) {
  var(x,na.rm=T)
}

## Function test

fnTest <- function(vMean,vVar,iN,dA=1) {
  iM <- length(vMean)
  if (iM>2) {
    mA <- cbind(rep(1,iM-1),diag(-rep(1,iM-1)))
    mVar <- diag(vVar[-1])
    mVar2 <- mVar+vVar[1]
    dTest <- drop(iN*t(mA%*%vMean)%*%solve(mVar2)%*%(mA%*%vMean))
  } else {
    dTest <- iN*(vMean[1]-vMean[2])^2/(vVar[1]+vVar[2])
  }
  if (iM>1) {
    dQuan <- qchisq(1-dA,iM-1)
    dPVa <- 1-pchisq(dTest,iM-1)
  } else {
    dQuan <- 0
    dPVa <- 1
  }
  dResult <- 1*(dTest>dQuan)
  if (is.na(dResult)) dResult <- 0
  return(list(test=dTest,q=dQuan,p=dPVa,result=dResult))
}

## Function elimination procedure

fnElim	<- function(vMean,vVar,vIndex) {
  iM <-	which.max(vMean)
  return(list(mean=vMean[-iM],var=vVar[-iM],indices=vIndex[-iM],elim=vIndex[iM]))
}

## Function MCS

fnMCS <- function(vMean,vVar,vIndex,dA,verbose=1) {
  iI <- length(vIndex)
  iStop <- fnTest(vMean,vVar,dA)$result
  mPValue <- matrix(0,ncol=2,nrow=iCoP)
  for (i in 1:iI) {
    lTest <- fnTest(vMean,vVar,iN,dA)
    iStop <- lTest$result
    if (iStop!=1) break
    lElim <- fnElim(vMean,vVar,vIndex)
    vMean <- lElim$mean
    vVar <- lElim$var
    vIndex <- lElim$indices
    if (verbose==1) cat("Eliminated COP: ",lElim$elim,"; ","p-value: ", lTest$p,"\n", sep="")
    if (verbose==1) cat("Remaining COPs:",vIndex,"\n")
    mPValue[i,] <- c(lElim$elim,lTest$p)
  }
  if (i==iI) {
    mPValue[iI,1] <- setdiff(vIndex,mPValue[1:(iCoP-1),1])
    mPValue[iI,2] <- 1
  } else {
    mPValue <- mPValue[1:(i-1),]
    mPValue[,2] <- cummax(mPValue[,2])
  }
  return(mPValue)
}



