clear all
clean_dynare_files;



%----------------------------------------------------------------
% 1. Paramethers declaration
%----------------------------------------------------------------
%NK non linear DSGE 

siggma = 1;
varphi=5;
phi_pi = 1.5;
phi_y  = 0.125;
theta=3/4;
rho_nu =0.5;
rho_z  = 0.5;
rho_a  = 0.9;
betta  = 0.99;
eta  =3.77; %footnote 11, p. 115
alppha=1/4;
epsilon=9;
tau=0; %//1/epsilon;

%----------------------------------------------------------------
% 2. Obtain model simulation
%----------------------------------------------------------------
n_MC=2;
n_CoP=2;

periods=193;
params = [alppha,betta,rho_a,rho_nu,rho_z,siggma,varphi,phi_pi,phi_y,eta,epsilon,theta,tau];

Simul_logY=zeros(periods,n_MC,n_CoP);
Simul_Pi=zeros(periods,n_MC,n_CoP);
Simul_R=zeros(periods,n_MC,n_CoP);
for i=1:n_MC
    for j=1:n_CoP
        params_t= params+randn(1,13)*0.01;
        Simul_data_0=Simul_data(params_t,periods);
        Simul_logY(:,i,j)=table2array(Simul_data_0(:,'log_y'));
        Simul_Pi(:,i,j)=table2array(Simul_data_0(:,'Pi'));
        Simul_R(:,i,j)=table2array(Simul_data_0(:,'R'));
    end    
end 

save Simul_logY.mat Simul_logY
save Simul_Pi.mat Simul_Pi
save Simul_R.mat Simul_R
%----------------------------------OLD CODE--------------------------------
%{

%----------------------------------------------------------------
% 1. Paramethers declaration
%----------------------------------------------------------------
%NK non linear DSGE 

siggma = 1;
varphi=5;
phi_pi = 1.5;
phi_y  = 0.125;
theta=3/4;
rho_nu =0.5;
rho_z  = 0.5;
rho_a  = 0.9;
betta  = 0.99;
eta  =3.77; %footnote 11, p. 115
alppha=1/4;
epsilon=9;
tau=0; %//1/epsilon;


%RBC DSGE
%{
alpha = 0.3;  
beta  = 0.99;
sigma = 1;
delta = 0.025;
rhoa = 0;
%}

%NK DSGE
%values to stractural parameters
%{
beta=0.99;      %discount factor
om=0.75;        %firm able to readjust
eta=1;          %log linear utility function
si=1;           %log linear utility function
phi_pi=1.5;
phi_x=0.1;
rho_a=0.9; %persistence technological shock
rho_v=0.5; %persistence monetary shock
rho_u=0.3; %
%}


%----------------------------------------------------------------
% 2. Obtain model simulation
%----------------------------------------------------------------
%params= [alpha,beta,sigma,delta,rhoa];
params = [alppha,betta,rho_a,rho_nu,rho_z,siggma,varphi,phi_pi,phi_y,eta,epsilon,theta,tau];
%params=[beta,om,eta,si,phi_pi,phi_x,rho_a,rho_v,rho_u];
periods=193;
Simul_data_0=Simul_data(params,periods,"NK_NL_DSGE");
%}
