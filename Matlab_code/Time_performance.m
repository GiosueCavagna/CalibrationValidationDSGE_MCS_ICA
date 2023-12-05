clear all
n_MC=5;
n_CoP=15;
periods=200;
tic
[c_e1,c_e2,q]=Data_simulation(n_MC,n_CoP,periods);
toc

%time_for_simulation =1.210449
%total time=~ time for simulation *n_MC*n_CoP
%Elapsed time is 500,500 200 -> 197343.640984
%(1.210449*(500*500))/(3600*24)=~3,5 days