# Simulation Part
The simulation part is all developed in Matlab and Dynare.

Data_simulation is the script that have to be run in order to obtain the simulation, it create a a dataset  with the simulation of each varibale of interest, containing the desired number of MonteCarlo simulation (MC) for each of the desired Configuration of Parameters (CoP).

This script firstly produce a Sobol sampling for the CoP, then it is run the MC of the NK_NL_DSGE model (Gali 2015) using the chosen CoP. Finally it store all the simulation of the varible of interest in a .mat file. This file will be analyzed in R.

## Pay attention:
In order to have non normal distributed shock it is necessary to change the simult.m function inside the Dynare directory with the one present in this repository which contain the easy necessary change that allow to get laplacian shock.
