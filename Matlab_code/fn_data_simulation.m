function [CoP_store]=fn_data_simulation(n_MC,n_CoP,periods)
%This function take as input the number of MC simulation, the number of Cop and the number of periods. It return the Simulation dataset of log_Y, Pi and R.
    %clear all
    clean_dynare_files;
    
    
    rng(1,"twister");
    %----------------------------------------------------------------
    % 1. Paramethers Declaration
    %----------------------------------------------------------------

    %NK non linear DSGE standard parameters 
    alppha  =   1/4;  %1  % capital share [0,1]
    betta   =   0.99; %2  % discount factor (0,1)
    rho_a   =   0.9;  %3  % autocorrelation tech schok [0,1)
    rho_nu  =   0.5;  %4  % autocorrelation monetary policy [0,1)
    rho_z   =   0.5;  %5  % autocorrelation monetary demand [0,1)
    siggma  =   1;    %6  % inverse of elasticity of temporal sobstitution of consumption (0,+inf) 
    varphi  =   5;    %7  % inverse of elasticity of temporal sobstitution of labor (0,+inf) 
    phi_pi  =   1.5;  %8  % inflation feedback taylor rule [0,+inf)
    phi_y   =   0.125;%9  % output feedback taylor rule [0,+inf) 
    epsilon =   9;    %10  %demand elasticity (0,+inf)
    theta   =   3/4;  %11 % Calvo parameter [0,1]
    %I will not calibrate eta
    eta     =    3.77;%13          %params_t(13); %3.77;   % semielasticity of money demand 

    params = [alppha,betta,rho_a,rho_nu,rho_z,siggma,varphi,phi_pi,phi_y,epsilon,theta,eta];
    
    %This CoP do not satisfixes the rank consdition.
    %params =[0.562500000000000,0.437500000000000,0.187500000000000,0.812500000000000,0.687500000000000,1.06250000000000,8.93750000000000,0.312500000000000,0.312500000000000,6.87500000000000,0.562500000000000,0,3.77000000000000];

    q=sobolset(11);
    %----------------------------------------------------------------
    %n_MC=2;
    %n_CoP=2;
    %periods=193;

    %Inizialization of the data set
    Simul_logY=zeros(periods,n_MC,n_CoP);   
    Simul_Pi=zeros(periods,n_MC,n_CoP);
    Simul_R=zeros(periods,n_MC,n_CoP);
    
    VarNames ={'C';'W_real';'Pi';'A';'N';'R';'realinterest';'Y';'Q';'Z';'S';'Pi_star';'x_aux_1';'x_aux_2';'MC';'nu';'P'};%;'M_real';'i_ann';'pi_ann';'r_real_ann';'log_m_nominal';'log_y';'log_W_real';'log_N';'log_P';'log_A';'log_Z'};
    
    %Creation txt file that will be importend during dynare execution in
    %order to insert the number of periods
    fid=fopen('util.txt','w+');
    fprintf(fid,'@#define T = %u\n',periods); 
    fclose(fid);
    
    %inizialization counter
    %c_e1=0; %other error
    %c_e2=0; %rank condition not verified

    %Inizialization store_CoP
    CoP_store=zeros(n_CoP,length(params));
    %Main loop simulation
    for i=1:n_CoP

        %definition of parameters
        if i==1 %The first parametrization is the standard parametrization
            params_t=params; 
        else    %Other Cop with constraint of specific parameters
            params_t= [q(i+1,1:11),eta];
            params_t(1)=(params_t(1)*params(1)*0.40)+params(1)*0.8;     %limitation at -20% + 20%
            params_t(2)=(params_t(2)*params(2)*0.20)+params(2)*0.8;     %limitation at -20% + 0%
            params_t(3)=(params_t(3)*params(3)*0.31)+params(3)*0.8;     %limitation at -20% + 11%
            params_t(4)=(params_t(4)*params(4)*0.40)+params(4)*0.8;     %limitation at -20% + 20%
            params_t(5)=(params_t(5)*params(5)*0.40)+params(5)*0.8;     %limitation at -20% + 20%
            params_t(6)=(params_t(6)*params(6)*0.40)+params(6)*0.8;     %limitation at -20% + 20%
            params_t(7)=(params_t(7)*params(7)*0.40)+params(7)*0.8;     %limitation at -20% + 20%
            params_t(8)=(params_t(8)*params(8)*0.40)+params(8)*0.8;     %limitation at -20% + 20%
            params_t(9)=(params_t(9)*params(9)*0.40)+params(9)*0.8;     %limitation at -20% + 20%
            params_t(10)=(params_t(10)*params(10)*0.40)+params(10)*0.8; %limitation at -20% + 20%
            params_t(11)=(params_t(11)*params(11)*0.40)+params(11)*0.8; %limitation at -20% + 20%
        end

        %MC simulation
        for j=1:n_MC
            my_seed=randi(n_MC*n_CoP*2);
            fprintf('\n\n\n\n\nNUMBER OF CoP = %d',i)
            fprintf('\nNUMBER OF MC = %d\n\n\n\n\n',j)

            save myparam_values.mat params_t;
            save my_seed.mat my_seed;

            error=0;
    
            try
                dynare NK_NL_DSGE.mod;
            catch
                %fprintf("Dynare have an error");
                error=1;
            end
    
            if error==1 %If dynare gives an error due to a weird CoP then Zeros is saved in the dataset
                Dd=zeros(periods,length(VarNames));
                Simul_data_t= array2table(Dd, 'VariableNames', VarNames);
                %c_e1=c_e1+1;
            else
                load oo_.mat oo_
                Dd=[oo_.endo_simul]';
                if isempty (Dd)
                    %c_e2=c_e2+1;

                    Dd=zeros(periods,length(VarNames));
                    Simul_data_t= array2table(Dd, 'VariableNames', VarNames);
                    for k=1:n_MC
                        Simul_Y(:,k,i)=table2array(Simul_data_t(:,'Y'));
                        Simul_Pi(:,k,i)=table2array(Simul_data_t(:,'P'));
                        Simul_R(:,k,i)=table2array(Simul_data_t(:,'R'));
                    end
                    break;
                end
                Simul_data_t= array2table(Dd, 'VariableNames', VarNames);
                delete oo_.mat
            end
    
            %Simul_logY(:,j,i)=table2array(Simul_data_t(:,'log_y'));
            Simul_Y(:,j,i)=table2array(Simul_data_t(:,'Y'));
            %Simul_Pi(:,j,i)=table2array(Simul_data_t(:,'Pi'));
            Simul_P(:,j,i)=table2array(Simul_data_t(:,'P'));
            Simul_R(:,j,i)=table2array(Simul_data_t(:,'R'));
        end  
        CoP_store(i,:)=params_t;
    end 

        
    %it is remove the useless file created during the MC
    delete NK_NL_DSGE.log
    %delete NK_NL_DSGE.jnl
    delete util.txt
    delete myparam_values.mat
    delete my_seed.mat
        
    %it is saved the results of the simulation 
    %save Simulated_Data/Simul_logY.mat Simul_logY
    save Simulated_Data/Simul_Y.mat Simul_Y
    save Simulated_Data/Simul_P.mat Simul_P
    save Simulated_Data/Simul_R.mat Simul_R

    
    
    %----------------------------------MODELS STANDARD PARAMETERS--------------------------------
    %{
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
    %}
    
    