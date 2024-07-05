clc
clear all

n_MC=200;
n_CoP=500;
periods=300;
tic
[CoP_store]=fn_data_simulation(n_MC,n_CoP,periods);
info_simul=[n_MC,n_CoP, periods];
save Simulated_Data/info_simul.mat  info_simul
save Simulated_Data/CoP_store.mat  CoP_store
toc

%time_for_simulation =~ 1.2
%total time=~ time for simulation *n_MC*n_CoP
%Elapsed time is 500,500 200 -> 197343s
%(1.2*(500*500))/(3600*24)=~3,5 days


%load('/Users/giosuecavagna/Documents/CalibrationValidationDSGE_MCS_ICA/Matlab_code/Simulated_Data/Simul_Y.mat')
%Simul_Y(250:end,:,2)