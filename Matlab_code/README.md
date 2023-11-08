#Simulation Part
The simulation part is all developed in Matlab and Dynare.

Step1.m is the script that have to be run in order to obtain the simulation.

Inside it is called the function Simul_data.m which is a function that take as imput the parameters, the number of periods that have to be simulated and the name of the model that it is wanted to be simulated, and gives as output the dataset containing all the simulation of the endogenous variable.

## Pay attention:
In order to have non normal distributed shock it is necessary to change the simult.m function inside the Dynare directory with the one present in this repository which contain the easy necessary change that allow to get laplacian shock.