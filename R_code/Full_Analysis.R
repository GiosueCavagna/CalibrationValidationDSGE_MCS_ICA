rm(list=ls());

#LIBRARIES-----
library(lubridate)
library(dplyr)
library(moments)
library(fastICA)
library(clue)
library(steadyICA)
library(data.table)
library(R.matlab)
library(devtools)
#install_github("STieleman/svarCAL")
library(svarCAL)

#FUNCTION----

#Function to Compute the mininimum distance mixing matrix from the Real data mixing matrix
frobICA_mod<-function (M1 = NULL, M2 = NULL, S1 = NULL, S2 = NULL, standardize = FALSE){
  tfun = function(x) all(x == 0)
  if (is.null(M1) && is.null(M2) && is.null(S1) && is.null(S2))
    stop("need to supply either M1 and M2 or S1 and S2")
  if (!is.null(M1) && !is.null(M2) && !is.null(S1) && !is.null(S2)) {
    stop("provide either (M1 and M2) or (S1 and S2) but not both (M1,M2) and (S1,S2)")
  }
  if (!is.null(M1) && nrow(M1) > ncol(M1))
    stop("The input appears to be S1 and S2, but the arguments were not specified; re-run with S1=<object> and S2=<object>")
  if (is.null(M1)) {
    nS = nrow(S1)
    if (nS != nrow(S2))
      stop("S1 and S2 must have the same number of rows")
    if (sum(apply(S1, 2, tfun)) + sum(apply(S2, 2, tfun)))
      stop("frobICA not defined when S1 or S2 has a column of all zeros")
    if (standardize) {
      S1 = scale(S1)
      S2 = scale(S2)
    }
    
    #START the real frobICA function
    p = ncol(S1)
    q = ncol(S2)
    if (p < q) {
      S1 = cbind(S1, matrix(0, nS, (q - p)))
    }
    if (q < p) {
      S2 = cbind(S2, matrix(0, nS, (p - q)))
    }
    Stemp = matchICA(S = S1, template = S2)
    n.comp = max(q, p)
    indices = c(1:n.comp)[!(apply(Stemp, 2, tfun) | apply(S2,
                                                          2, tfun))]
    return(sqrt(sum((Stemp[, indices] - S2[, indices])^2))/sqrt(nS *
                                                                  min(p, q)))
  }
  else {
    if (sum(apply(M1, 1, tfun)) + sum(apply(M2, 1, tfun)))
      stop("frobICA not defined when M1 or M2 has a row of all zeros")
    if (standardize) {
      temp = diag((diag(M1 %*% t(M1)))^(-1/2))
      M1 = temp %*% M1
      temp = diag((diag(M2 %*% t(M2)))^(-1/2))
      M2 = temp %*% M2
    }
    p = ncol(M1)
    if (p != ncol(M2))
      stop("M1 and M2 must have the same number of columns")
    d = nrow(M1)
    q = nrow(M2)
    n.comp = max(d, q)
    if (n.comp > p)
      warning("M should be d x p")
    if (d < q) {
      M1 = rbind(M1, matrix(0, (q - d), p))
    }
    if (q < d) {
      M2 = rbind(M2, matrix(0, (d - q), p))
    }
    l2.mat1 = l2.mat2 = matrix(NA, nrow = n.comp, ncol = n.comp)
    for (j in 1:n.comp) {
      for (i in 1:n.comp) {
        l2.mat1[i, j] = sum((M2[i, ] - M1[j, ])^2)
        l2.mat2[i, j] = sum((M2[i, ] + M1[j, ])^2)
      }
    }
    l2.mat1 = sqrt(l2.mat1)
    l2.mat2 = sqrt(l2.mat2)
    l2.mat = l2.mat1 * (l2.mat1 <= l2.mat2) + l2.mat2 * (l2.mat2 <
                                                           l2.mat1)
    map = as.vector(solve_LSAP(l2.mat))
    l2.1 = diag(l2.mat1[, map])
    l2.2 = diag(l2.mat2[, map])
    sign.change = -1 * (l2.2 < l2.1) + 1 * (l2.1 <= l2.2)
    perm = diag(n.comp)[, map] %*% diag(sign.change)
    M.perm = t(perm) %*% M1
    indices = c(1:n.comp)[!(apply(M.perm, 1, tfun) | apply(M2,
                                                           1, tfun))]
    return(list(perm=perm,frob_dist=sqrt(sum((M.perm[indices, ] - M2[indices, ])^2))/sqrt(p *
                                                                                            min(d, q))))
  }
}

#Function to compute the mixing matrix
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

#Function to compute the mean of all the MCs
fnMean <- function(x) {
  mean(x,na.rm=T)
}

#Function to compute the variance of all the MCs
fnVar <- function(x) {
  var(x,na.rm=T)
}

#Function to compute the test statistics of the CoP
fnTest <- function(vMean,vVar,nmr,dA=1) {
  iM <- length(vMean)
  if (iM>2) {
    mA <- cbind(rep(1,iM-1),diag(-rep(1,iM-1)))
    mVar <- diag(vVar[-1])
    mVar2 <- mVar+vVar[1]
    dTest <- drop(nmr*t(mA%*%vMean)%*%solve(mVar2)%*%(mA%*%vMean))
  } else {
    dTest <- nmr*(vMean[1]-vMean[2])^2/(vVar[1]+vVar[2])
    cat(dTest)
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

#Function to eliminate the worste CoP
fnElim	<- function(vMean,vVar,vIndex) {
  iM <-	which.max(vMean)
  return(list(mean=vMean[-iM],var=vVar[-iM],indices=vIndex[-iM],elim=vIndex[iM]))
}

#Function to perform the MCS through test and elimination procedure
fnMCS <- function(vMean,vVar,vIndex,dA,verbose=1) {
  iI <- length(vIndex)
  iStop <- fnTest(vMean,vVar,dA)$result
  mPValue <- matrix(0,ncol=2,nrow=ncp)
  for (i in 1:iI) {
    if (verbose==1) cat(i,"-loop \n")
    lTest <- fnTest(vMean,vVar,nmr,dA)
    iStop <- lTest$result
    if (iStop!=1) break
    lElim <- fnElim(vMean,vVar,vIndex)
    vMean <- lElim$mean
    vVar <- lElim$var
    vIndex <- lElim$indices
    if (verbose==1) cat("Eliminated COP: ",lElim$elim,"; ","p-value: ", lTest$p,"\n", sep="")
    if (verbose==1) cat("Remaining COPs:",vIndex,"\n\n\n")
    mPValue[i,] <- c(lElim$elim,lTest$p)
  }
  if (i==iI) {
    mPValue[iI,1] <- setdiff(vIndex,mPValue[1:(ncp-1),1])
    mPValue[iI,2] <- 1
  } else {
    mPValue <- mPValue[1:(i-1),]
    mPValue[,2] <- cummax(mPValue[,2])
  }
  return(mPValue)
}

# Function to compute levels from growth 
ratesfnLev <- function(data,L_0,mu_z) { 
  levels <- rep(0,length(data))#,ncol=ncol(data)) 
  #for (i in 1:ncol(data)) { 
  growth_rates <- diff(data)#,i],) 
  levels[1] <- L_0 
  for (j in 2:(length(growth_rates) + 1)) { 
    levels[j] <- levels[j-1]+growth_rates[j-1] + mu_z 
  }
  return(levels)}



#CODE----
#REAL DATA-----
#-Data import and arrangement----
Rdata0=read.csv("Data/2023-09.csv")
Rdata0$sasdate=mdy(Rdata0$sasdate)
Rdata0=Rdata0[!is.na(Rdata0$sasdate),]

#%slicing 1973:Q1 - 2019:Q1 (247)
Q1_1971=ymd("1971/03/01");
#%Q1_2019=247;
Q1_2019=ymd("2019/03/01");

Rdata=Rdata0;
Rdata = Rdata %>% filter(sasdate >= Q1_1971 & sasdate <= Q1_2019)

#N0=read.csv("Data/CLF16OV.csv") #Workforce population
N0=read.csv("Data/POPTHM.csv") #Population
N0$DATE=ymd(N0$DATE)
N= N0 %>% filter( DATE>= Q1_1971 & DATE <= Q1_2019)
N= N %>% filter( month(DATE)==3 |month(DATE)==6 |month(DATE)==9 |month(DATE)==12 )

#-Preparation of data for VAR on diff----

#GDP per Capita
#Y=log((Rdata$GDPC1*10^6/N$POPTHM))*10;
#Ydiff1=log(Rdata0$GDPC1[Rdata0$sasdate==ymd("1971/03/01")]/N0$CLF16OV[N0$DATE==ymd("1971/03/01")])-log(Rdata0$GDPC1[Rdata0$sasdate==ymd("1970/12/01")]/N0$CLF16OV[N0$DATE==ymd("1970/12/01")]);
#Ydiff2=diff(Y);
#Ydiff=c(Ydiff1,Ydiff2)
#rm(Y,Ydiff1,Ydiff2);

#Inflation
#P1=Rdata0$CPIAUCSL[Rdata0$sasdate==ymd("1971/03/01")]-Rdata0$CPIAUCSL[Rdata0$sasdate==ymd("1970/12/01")];
#P2=diff(Rdata$CPIAUCSL)
#P=c(P1,P2)
#rm(P1,P2);

#Nominal Interest Rate
#R=(1+Rdata$FEDFUNDS/100)^(1/4)

#D=data.frame(Column1 = Ydiff, Column2 = P, Column3 = R)
#colnames(D) <- c("Ydiff", "P", "R")

#-Preparation of data for VAR on levels----
n_periods=length(Rdata$sasdate)
index_periods=1:n_periods

#GDP per Capita
Y=log((Rdata$GDPC1*10^6/N$POPTHM))*10;
lm_Y=lm(Y~index_periods)

#Inflation
P=log(Rdata$CPIAUCSL)*10
lm_P=lm(P~index_periods)

#Nominal Interest Rate
R=(1+Rdata$FEDFUNDS/100)^(1/4)

#Creation DF
D=data.frame(Column1 = Y, Column2 = P, Column3 = R)
colnames(D) <- c("Y", "P", "R")

#-Data elaboration----

infocrit = VARselect(D, lag.max=6, type="const")
lag = infocrit$selection["AIC(n)"]
varest = VAR(D, p=lag, type="const")
ures = resid(varest)
AA = Acoef(varest) # VAR coefficients
k = ncol(ures)

# Test Gaussianity 
jarque_result = c(jarque.test(ures[,1])$p.value,jarque.test(ures[,2])$p.value,jarque.test(ures[,3])$p.value)
#we refuse the null -> not normal at 5%

A_rw = fAp_fastICA(ures,sseed=46) # rw mixing matrix
refM = A_rw

rm(list=setdiff(ls(), c("n_periods","refM","lag","lm_Y","lm_P","frobICA_mod","fAp_fastICA","fnMean","fnVar","fnTest","fnElim","fnMCS","ratesfnLev")))

#SIMULATED DATA----
CoP_store=readMat("~/Documents/CalibrationValidationDSGE_MCS_ICA/Matlab_code/Simulated_Data/Heavy_simul/CoP_store.mat") 
info_simul=readMat("~/Documents/CalibrationValidationDSGE_MCS_ICA/Matlab_code/Simulated_Data/Heavy_simul/info_simul.mat")


nmr= info_simul$info.simul[1] #500 # number of Monte Carlo runs
ncp= info_simul$info.simul[2] #500 # number of configuration of parameters
tau= info_simul$info.simul[3] #200 # simulation length

k = nrow(refM) # number of variables

Simul_Y=readMat("~/Documents/CalibrationValidationDSGE_MCS_ICA/Matlab_code/Simulated_Data/Heavy_simul/Simul_Y.mat")
Simul_P=readMat("~/Documents/CalibrationValidationDSGE_MCS_ICA/Matlab_code/Simulated_Data/Heavy_simul/Simul_P.mat") #readMat("Data/Simul_PI.mat") #
Simul_R=readMat("~/Documents/CalibrationValidationDSGE_MCS_ICA/Matlab_code/Simulated_Data/Heavy_simul/Simul_R.mat") # readMat("Data/Simul_R.mat") #
#[tau,nmr,ncp]

Simul_Y=Simul_Y$Simul.Y
Simul_P=Simul_P$Simul.P
Simul_R=Simul_R$Simul.R

#-Creation of simulation DF and removals of first errors----
if (n_periods>tau) {
  n_periods=tau
}

D = list()
i=1
CoP_err=matrix(NaN,3,ncp);
for(ii in 1:ncp){
  CoP_err[c(1,2),ii]=c(ii,"MatlabErr");
  #Check if the simulation had problems in Matlab
  if (!( is.nan(Simul_Y[tau,nmr,ii])      | 
         abs(Simul_Y[tau,nmr,ii])>10e+2   |
         abs(Simul_Y[tau,nmr,ii])<10e-14  ) ){ #if there are no error
    D[[i]] = list()
    D[[i]]$list = list()
    for (j in 1:nmr){
      D[[i]]$list[[j]] = as.data.frame(cbind(Simul_Y[((tau-n_periods)+1):tau,j,ii],Simul_P[((tau-n_periods)+1):tau,j,ii],Simul_R[((tau-n_periods)+1):tau,j,ii]))
      colnames(D[[i]]$list[[j]]) = c('Y','P','R')
    }
    CoP_err[c(1,2,3),ii]=c(ii,"NoErr",i);
    i=i+1
  }
  print(ii)
}

ncp=i-1 #number of CoP that have passed the first filtration for error
(ncp)
if (ncp==0){
  stop("None configuration pass the test, so the script end here")
}
#-Adding the trend to the simulated data----
for (i in 1:ncp){
  for (j in 1:nmr){
    D[[i]]$list[[j]]$Y=ratesfnLev(c(unlist(D[[i]]$list[[j]]$Y)),lm_Y$coefficients[1],lm_Y$coefficients[2]) 
    D[[i]]$list[[j]]$P=ratesfnLev(c(unlist(D[[i]]$list[[j]]$P)),lm_P$coefficients[1],lm_P$coefficients[2]) 
  }
}

#-Estimation minimum distance mixing matrix procedure----

m_VAR=as.list(1:ncp)#lapply(D[[CoP_err[1,which(CoP_err[3,]==1)]]][[1]],VAR,p=lag)
VarErr_loc1=array(dim=0)
VarErr_loc2=rep(0,ncp)
#-1 I estimate the var and compute the minimum number of mcr

#Function needed otherwise R is not able to handle the variable m_VAR that
#saving everything instead of just residuals and coeff became 21 GB heavy
fn_m_VAR = function(x, p_lag) {
  A=VAR(x, p = p_lag)
  residuals = residuals(A)
  coefficients = coefficients(A)
  results <- list(residuals = residuals, coefficients = coefficients)
  return(results)
}

for (i in 1:ncp){
        print(i)
        err1=0
        m_VAR[[i]] <- tryCatch(
          { lapply(D[[i]][[1]],fn_m_VAR, p_lag=lag)},
          error = function(e) {
            err1<<-1
            return(NaN) #In case of problem I assign a Nan value
          }
        )
        
        if (err1==1){
          VarErr_loc1=c(VarErr_loc1,i)
          VarErr_loc2[i]=1
        }
}

#other solution less elegant----------
# D1=D[1:50];
# D2=D[51:100];
# D3=D[101:150];
# D4=D[151:200];
# D5=D[201:250];
# D6=D[251:300];
# D7=D[301:350];
# D8=D[351:400];
# D9=D[401:ncp];
# m_VAR=as.list(1:50)
# 
# n_n_var=1
# applica_var=function(x) {
#   print(n_n_var)
#   n_n_var<<-n_n_var+1
#   return(lapply(x$list,VAR,p=lag))
# }
# 
# 
# n_quantile=21
# quan=round(quantile(1:ncp,probs =seq(0,1, 1/n_quantile)))
# for(j in 1:n_quantile){
#   if(j==n_quantile){
#     m_VAR=as.list(1:50)
#     m_VAR=lapply(D[quan[j]:(quan[j+1]-1)], applica_var);
#     save(m_VAR, file = paste0("your_variable_q",j,".RData"))
#     # Remove the variable from the environment
#     rm(m_VAR)
#   }
#   else{
#     m_VAR=as.list(1:50)
#     m_VAR=lapply(D[quan[j]:(quan[j+1]-1)], applica_var);
#     save(m_VAR, file = paste0("your_variable_q",j,".RData"))
#     # Remove the variable from the environment
#     rm(m_VAR)
#   }
#   print("one quantile is gone!")
# }

###--------

for (i in VarErr_loc1){
  CoP_err[c(2,3),which(CoP_err[3,]==i)]=c("VarErr",NaN);
}
ncp=ncp-length(VarErr_loc1)
CoP_err[3,which(CoP_err[3,]!=NaN)]=1:ncp
#I remove the problematic Var and update the ncp
m_VAR=m_VAR[!VarErr_loc2]



#2- I Check the normality and save the number of MC how pass JB test

# Function to compute JB test to m_res
fnJB = function(x){
  A=x$residuals
  return(c(jarque.test(A[,1])$p.value,jarque.test(A[,2])$p.value,jarque.test(A[,3])$p.value))
}

fnJBpass=function(x,alppha=NULL){
  nmc=rep(0,length(x))
  if (is.null(alppha)){alppha=0.05}
  for (i in 1:length(x)){
   if(sum(is.nan(x[[i]]))==0){
      if(sum(x[[i]]<alppha)==length(x[[i]])){
        nmc[i]=1
      }
    }
  }
  return(nmc)
}


nmr_pass = as.list(1:ncp)
JBErr_loc1=array(dim=0)
JBErr_loc2=rep(0,ncp)
JB_threshold=70
for(i in 1:ncp) {
  print(i)
  #nmr_pass[[i]] = rep(0,nmr)
  #Test of non normality
  Non_Nomarlity_test=lapply(m_VAR[[i]],fnJB)
  nmr_pass[[i]]=fnJBpass(Non_Nomarlity_test)
  if (sum(nmr_pass[[i]])<JB_threshold) {
    JBErr_loc1=c(JBErr_loc1,i)
    JBErr_loc2[i]=1
  }
}
#CoP_err[3,which(as.numeric(CoP_err[3,])>=i)]=as.numeric(CoP_err[3,which(as.numeric(CoP_err[3,])>=i)])-1


for (i in JBErr_loc1){
  CoP_err[c(2,3),which(CoP_err[3,]==i)]=c("JBErr",NaN);
}

ncp=ncp-length(JBErr_loc1)
CoP_err[3,which(CoP_err[3,]!=NaN)]=1:ncp
m_VAR=m_VAR[!JBErr_loc2]

#4- Estimation minimum distance mixing matrix procedure----
Mres = as.list(1:ncp)
ICAErr_loc1=array(dim=0)
ICAErr_loc2=rep(0,ncp)

for(i in 1:ncp) {
  err2=0
  Mres[[i]] <- list()
  Mres[[i]]$m_ica <- as.list(1:nmr) # it delivers a mixing matrix for each monte carlo run
  Mres[[i]]$f_dist <- rep(NA,nmr)
  Mres[[i]]$ica_perm <- as.list(1:nmr)
  
  #estimate mixing matrix using fastICA
  Mres[[i]]$m_ica <-  tryCatch(
    { A=lapply(m_VAR[[i]],function(x){x$residuals})
      lapply(A,fAp_fastICA,sseed=46)},# mixing matrices of model data
    error = function(e) {
      err2<<-1
      return(NaN) 
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
  
  
  if (err2==1){
    ICAErr_loc1=c(ICAErr_loc1,i)
    ICAErr_loc2[i]=1
  }
  print(i)
  }
 

for (i in ICAErr_loc1){
  CoP_err[c(2,3),which(CoP_err[3,]==i)]=c("ICAErr",NaN);
}

ncp=ncp-length(ICAErr_loc1)
CoP_err[3,which(CoP_err[3,]!=NaN)]=1:ncp


#Computation Frobenius distance
fb_dist = matrix(0,ncp,nmr)
for (i in 1:ncp) {
  fb_dist[i,] = Mres[[i]]$f_dist # MDI for each CoPs and MC runs
}

#MODEL CONFICENSE SET----

vMean.val = apply(fb_dist,1,fnMean)
vVar.val = apply(fb_dist,1,fnVar)
mVal = cbind(1:ncp,vMean.val,vVar.val)
mVal = mVal[order(mVal[,2],decreasing=T),]


dA = 1
alpha=0.05 #significance
#Check if there is enough MC in order to have a variance different from 0
if (sum(vVar.val==0)==0){
  mPValue.val = fnMCS(vMean.val,vVar.val,1:ncp,dA,verbose=1)
  mPValue.val = cbind(mPValue.val,mVal[,2:3])
  vCoP.pass = mPValue.val[which(mPValue.val[,2]>alpha)]  
} else {
  warning("The variance of at least one CoP is equal to 0")
  vCoP.pass = which(vVar.val==0)
}

n_pass=length(vCoP.pass)
CoP_pass=matrix(0,n_pass,size(CoP_store$CoP.store)[2])
for (i in 1:n_pass){
  CoP_pass[i,]= CoP_store$CoP.store[which(CoP_err[3,]==vCoP.pass[i]),]
}
CoP_colnames=c("alppha","betta","rho_a","rho_nu","rho_z","siggma","varphi","phi_pi","phi_y","epsilon","theta","eta")
CoP_pass = data.frame(CoP_pass)
colnames(CoP_pass) = CoP_colnames
#CoP_pass=t(CoP_pass)

#Computation percentage of discarded CoP due to error
err_perc=(info_simul$info.simul[2]-ncp)/info_simul$info.simul[2]


#rm(list=setdiff(ls(), c("CoP_pass","mPValue.val","lag","refM","CoP_store", "CoP_err","err_perc")))

MCS <- paste("MCS_calib.csv",sep=",")
write.csv(mPValue.val,MCS,col.names=T)

COP <- paste("CoP_calib.csv",sep=",")
write.csv(CoP_pass,COP,col.names=T)

