function [Simul_data,counter_e]=Simul_data (params, periods)%, nmodel)
%This function take as imput the parameters, the number of periods that have to be simulated and the name of the model that it is wanted to be simulated, and gives as output the dataset containing all the simulation of the endogenous variable.

myparam= params;
save myparam_values.mat myparam;

VarNames ={'C';'W_real';'Pi';'A';'N';'R';'realinterest';'Y';'Q';'Z';'S';'Pi_star';'x_aux_1';'x_aux_2';'MC';'M_real';'i_ann';'pi_ann';'r_real_ann';'P';'log_m_nominal';'log_y';'log_W_real';'log_N';'log_P';'log_A';'log_Z';'nu'};

fid=fopen('util.txt','w+');
fprintf(fid,'@#define T = %u\n',periods); 
fclose(fid);

error=0;
counter_e=0;

try
    dynare NK_NL_DSGE.mod;
catch
    %fprintf("Dynare ha dato un errore");
    error=1;
end

if error==1
    Dd=NaN(periods,length(VarNames));
    Dd= array2table(Dd, 'VariableNames', VarNames);
    error=0;
    counter_e=counter_e+1;

else
    load oo_.mat oo_
    Dd=[oo_.endo_simul]';
    Dd= array2table(Dd, 'VariableNames', VarNames);
    delete oo_.mat
    delete NK_NL_DSGE.log
    delete util.txt
  
end
  

delete myparam_values.mat

    
%Use the writetable function to save the table as a CSV file
%writetable(Dd, 'Simul_data_.csv');

%Simul_data = readtable('Simul_data_.csv');
Simul_data=Dd;

end

%----------------------------------OLD CODE--------------------------------
%{

myparam= params;
save myparam_values.mat myparam;


fid=fopen('util.txt','w+');
fprintf(fid,'@#define T = %u\n',periods); 
fclose(fid);

%fid1=fopen('utill.m','w+');
%fprintf(fid1,'try\n    dynare '+nmodel+'.mod\ncatch\n    fprintf("Dynare ha dato un errore");\nend\n'+'delete '+nmodel+'.log'); 
fclose(fid1);

run('utill.m')

load oo_.mat oo_
Data=[oo_.endo_simul]';
VarNames=oo_.var_list;
Data = array2table(Data, 'VariableNames', VarNames);

% Use the writetable function to save the table as a CSV file
delete utill.m
delete util.txt
delete myparam_values.mat
delete oo_.mat

writetable(Data, 'Simul_data_'+nmodel+'.csv');


%Simul_data = readtable(Simul_data_'+nmodel+'.csv');
Simul_data=Data;

end

%}
