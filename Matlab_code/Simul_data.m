function Simul_data=Simul_data (params, periods)%, nmodel)
%This function take as imput the parameters, the number of periods that have to be simulated and the name of the model that it is wanted to be simulated, and gives as output the dataset containing all the simulation of the endogenous variable.

myparam= params;
save myparam_values.mat myparam;


fid=fopen('util.txt','w+');
fprintf(fid,'@#define T = %u\n',periods); 
fclose(fid);

try
    dynare NK_NL_DSGE.mod;
catch
    fprintf("Dynare ha dato un errore");
end

load oo_.mat oo_
Dd=[oo_.endo_simul]';
VarNames=oo_.var_list;
Dd= array2table(Dd, 'VariableNames', VarNames);


delete NK_NL_DSGE.log
delete util.txt
delete myparam_values.mat
delete oo_.mat

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
