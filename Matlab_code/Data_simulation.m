clear all
clean_dynare_files;



%----------------------------------------------------------------
% 1. Paramethers declaration
%----------------------------------------------------------------
%NK non linear DSGE 
%q=sobolset(13);
alppha  =   1/4;    % capital share [0,1]
betta   =   0.99;   % discount factor (0,1)
rho_a   =   0.9;    % autocorrelation tech schok [0,1)
rho_nu  =   0.5;    % autocorrelation monetary policy [0,1)
rho_z   =   0.5;    % autocorrelation monetary demand [0,1)
siggma  =   1;      % inverse of elasticity of temporal sobstitution of consumption (0,+inf) but we limitate it in [0.5,1,5] 
                    % since u is a CRRA if siggma->0 then the u is linear    %% q(i,6)+0.5
varphi  =   5;      % inverse of elasticity of temporal sobstitution of labor (0,+inf) but we limitate it in [0.5,9.5]
                    % (q(i,7)*9)+0.5
phi_pi  =   1.5;    % inflation feedback taylor rule [0,+inf) but we limitate it in [0,5] q(i,8)*5
phi_y   =   0.125;  % output feedback taylor rule [0,+inf) but we limitate it in [0,1]

theta   =   3/4;    % Calvo parameter [0,1]
epsilon =   9;      % demand elasticity (0,+inf)but we limitate it in [5,15] (q(i,11)*10)+5

%I will not calibrate tau and eta
tau     =   0;      %labor subsidy
eta     =   3.77;   % semielasticity of money demand 

%----------------------------------------------------------------
% 2. Obtain model simulation
%----------------------------------------------------------------
q=sobolset(11);
n_MC=2;
n_CoP=10;

periods=193;
params = [alppha,betta,rho_a,rho_nu,rho_z,siggma,varphi,phi_pi,phi_y,epsilon,theta,tau,eta];

Simul_logY=zeros(periods,n_MC,n_CoP);
Simul_Pi=zeros(periods,n_MC,n_CoP);
Simul_R=zeros(periods,n_MC,n_CoP);

VarNames ={'C';'W_real';'Pi';'A';'N';'R';'realinterest';'Y';'Q';'Z';'S';'Pi_star';'x_aux_1';'x_aux_2';'MC';'M_real';'i_ann';'pi_ann';'r_real_ann';'P';'log_m_nominal';'log_y';'log_W_real';'log_N';'log_P';'log_A';'log_Z';'nu'};

fid=fopen('util.txt','w+');
fprintf(fid,'@#define T = %u\n',periods); 
fclose(fid);

counter_e=0;
for i=1:n_CoP
    %definition of parameters
        
    if i==1
        params_t=params;
    else
        params_t= [q(i+1,1:11),tau,eta];
        params_t(6)= params_t(6)+0.5; %limitation siggma
        params_t(7)= (params_t(7)*9)+0.5; %limitation varphi
        params_t(8)=  params_t(8)*5; %limiation phi_pi
        params_t(11)=  (params_t(11)*10)+5; %limitation epsilon
    end

    for j=1:n_MC
        %Data Simulation
        save myparam_values.mat params_t;

        error=0;
  
        try
            dynare NK_NL_DSGE.mod;
        catch
            %fprintf("Dynare ha dato un errore");
            error=1;
        end
        
        if error==1
            Dd=NaN(periods,length(VarNames));
            Simul_data_t= array2table(Dd, 'VariableNames', VarNames);
            error=0;
            counter_e=counter_e+1;
        else
            load oo_.mat oo_
            Dd=[oo_.endo_simul]';
            Simul_data_t= array2table(Dd, 'VariableNames', VarNames);
            delete oo_.mat
  
        end

        Simul_logY(:,j,i)=table2array(Simul_data_t(:,'log_y'));
        Simul_Pi(:,j,i)=table2array(Simul_data_t(:,'Pi'));
        Simul_R(:,j,i)=table2array(Simul_data_t(:,'R'));
    end    
end 

delete NK_NL_DSGE.log
delete util.txt
delete myparam_values.mat

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
