# Simulation Part
The simulation part is all developed in Matlab and Dynare.

fn_data_simulation is the function that have to be run in order to obtain the simulation, it create a a dataset  with the simulation of each variable of interest, containing the desired number of MonteCarlo simulation (MC) for each of the desired Configuration of Parameters (CoP).

This script firstly produce a Sobol sampling for the CoP, then it is run the MC of the NK_NL_DSGE model (Gali 2015) using the chosen CoP. Finally it store all the simulation of the varible of interest in a .mat file. This file will be analyzed in R.

Inside Run_simulation there is the code with the declaration of variables needed for the function which simulate the model. Moreover the file needed for the R analysis are saved. 

## Pay attention:
In order to have non normal distributed shock it is necessary to change the simult.m function inside the Dynare directory with the one present in this repository which contain the easy necessary change that allow to get laplacian shock.
