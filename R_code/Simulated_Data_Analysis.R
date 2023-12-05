library(R.matlab)

#-------------------------------------------------------------------------------
#                                 SIMULATED DATA
#-------------------------------------------------------------------------------

Simul_logY=readMat("Data/Simul_logY.mat")
Simul_Pi=readMat("Data/Simul_Pi.mat")
Simul_R=readMat("Data/Simul_R.mat")

load(".../A_rw.Rdata") 
refM <- A_rw

tau= 200 # simulation length
nmr= 500 # number of Monte Carlo runs
ncp= 500 # number of configuration of parameters
#t <- length(106:tau) # time steps (without transient)
#wheredata<-'' # insert path containing the simulated data
k <- nrow(refM) # number of variables

#---------------------Creation of simulation dataset----------------------------
D = list()
for(i in 1:ncp){
  D[[i]] = list()
  D[[i]]$list = list()
  for (j in 1:nmr){
    D[[i]]$list[[j]] = as.data.frame(cbind(Simul_logY$Simul.logY[,j,i],Simul_Pi$Simul.Pi[,j,i],Simul_R$Simul.R[,j,i]))
    colnames(D[[i]]$list[[j]]) = c('log_Y','Pi','R')
  }
  print(i)
}

#---------------------Estimation procedure--------------------------------------
Mres = as.list(1:ncp)
for(i in 1:ncp) {
  Mres[[i]] <- list()
  Mres[[i]]$m_ica <- as.list(1:nmr) # it delivers a mixing matrix for each monte carlo run
  Mres[[i]]$f_dist <- rep(NA,nmr)
  Mres[[i]]$ica_perm <- as.list(1:nmr)
  
  #estimate VARs
  m_VAR <- lapply(D[[i]][[1]],VAR,p=2) # var of model data
  m_res <- lapply(m_VAR,residuals) # residuals of model data
  #VAR coefficients
  AA <- lapply(m_VAR,coef)
  #estimate mixing matrix using fastICA
  Mres[[i]]$m_ica <- lapply(m_res,fAp_fastICA,sseed=46) # mixing matrices of model data
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







t=200;#length(Rdata$sasdate)
i=0;
c=0;
c1=0
store=list();
for (i in 1:(length(Simul_R$Simul.R)/t) ) {
  c1=c1+1;
  SR=Simul_R$Simul.R[(t*(i-1)+1):(t*i)];
  if (is.nan(SR)[1]){
    c=c+1;
    store=c(store,i);
  }
}
e_perc=c/c1
c1-c
#250000 MC of which 28000 Nan




#--------------------------------------------Graphs-----------------------------
timespan=Rdata$sasdate;
SY=Simul_logY$Simul.logY[(t*(i-1)+1):(t*i)];
SP=Simul_Pi$Simul.Pi[(t*(i-1)+1):(t*i)];

par(mfrow = c(2, 3))
plot(timespan,Ydiff, type="l", main="Real GDP per Capita" , xlab="Quarter")
plot(timespan,P, type="l", main="Real Inflation", xlab="Quarter")
plot(timespan,R,type="l", main="Real Nominal Interest Rate",xlab="Quarter")

plot(timespan,SY, type="l", main="Simulated GDP per Capita" , xlab="Quarter")
plot(timespan,SP, type="l", main="Simulated Inflation", xlab="Quarter")
plot(timespan,SR,type="l", main="Simulated Nominal Interest Rate",xlab="Quarter")





#start_date=Rdata$sasdate[1]
#Sdata$date=seq(start_date, by = "quarter", length.out = length(Sdata$Y));

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


