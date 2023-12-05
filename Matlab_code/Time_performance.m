clear all
n_MC=5;
n_CoP=50;
periods=200;
tic
[c_e1,c_e2,q]=Data_simulation(n_MC,n_CoP,periods);
toc

%time=1.210449
%Elapsed time is 10,2, 193 20.493271seconds.
%Elapsed time is 500,500 200197343.640984
%(1.210449*(500*500))/(3600*24)=~3,5 days