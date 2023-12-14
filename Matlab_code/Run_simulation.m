clear all
n_MC=200;
n_CoP=1250;
periods=300;
tic
[CoP_store]=fn_data_simulation(n_MC,n_CoP,periods);
info_simul=[n_MC,n_CoP, periods];
save Simulated_Data/info_simul2.mat  info_simul
save Simulated_Data/CoP_store2.mat  CoP_store
toc

%time_for_simulation =1.210449
%total time=~ time for simulation *n_MC*n_CoP
%Elapsed time is 500,500 200 -> 197343.640984
%(1.210449*(500*500))/(3600*24)=~3,5 days