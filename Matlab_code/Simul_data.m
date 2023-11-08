function Simul_data=Simul_data (params, periods, nmodel)



myparam= params;
save myparam_values.mat myparam;


fid=fopen('util.txt','w+');
fprintf(fid,'@#define T = %u\n',periods); 
fclose(fid);

fid1=fopen('utill.m','w+');
fprintf(fid1,'try\n    dynare '+nmodel+'.mod\ncatch\n    fprintf("Dynare ha dato un errore");\nend\n'+'delete '+nmodel+'.log'); 
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


Simul_data = readtable('Simul_data.csv');

end
