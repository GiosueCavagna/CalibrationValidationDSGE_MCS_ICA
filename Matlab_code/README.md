# Read Me
The purpose of these scripts is running Montecarlo simulation with different configuration of paramether of a DSGE model exploting dynare.

Inside this repository is possible to find:

-	fn_data_simulation.m:
	it is the function that has to be run in order to obtain the simulation, it creates a dataset with the simulations of each variable of interest, 	containing the desired number of MonteCarlo simulation (n_MC) for each of      the desired Configuration of Parameters (n_CoP). This script firstly produce a Sobol sampling for the CoP, then it is run the simulation of the NK_NL_DSGE model (Gali 2015) using the chosen CoP. Finally it store all the      simulation of the varible of interest in a ".mat"   file.

- Run_simulation.m:
	it is the code with the declaration of variables needed for fn_data_simulation, and where this function is called. Moreover the file needed for the R analysis are saved.

-	NK_NL_DSGE.mod:
	it is the dynare file containting the eqaution of the model.

-	simult_NewVersion.m:
	it is the script that has to be substitued inside dynare-#.#/matlab in order to have Laplacian distributed shokcs.

- Clean_dyanre_files.m:
  it is a simple function that clean all the dynare leftovers.
  
