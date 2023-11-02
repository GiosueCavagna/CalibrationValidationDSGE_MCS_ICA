function Simul_data=Simul_data (params, periods, RBC)



myparam= params;
save myparam_values.mat myparam;


fid=fopen('util.txt','w+');
fprintf(fid,'@#define T = %u\n',periods); 
fclose(fid);

if RBC==true
    try
        dynare RBC.mod ;
    catch   
       fprintf('Dynare ha dato un errore');
    end
    delete RBC.log
else 
    try
        dynare DSGE.mod ;
    catch   
        fprintf('Dynare ha dato un errore');
    end
    delete DSGE.log
end

load oo_.mat oo_
Data=[oo_.endo_simul]';
VarNames=oo_.var_list;
Data = array2table(Data, 'VariableNames', VarNames);

% Use the writetable function to save the table as a CSV file



delete util.txt
delete myparam_values.mat
delete oo_.mat

writetable(Data, 'Simul_data.csv');


Simul_data = readtable('Simul_data.csv');

end
