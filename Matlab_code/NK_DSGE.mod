% basic NK model Gali 2008 chapter 3

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

%endogenous variable needed for the model
var x pi int a u v;

%endogenous variable we are interested to plot
var r y yn n;


varexo ea eu ev;


parameters beta om eta si phi_pi phi_x rho_a rho_v rho_u lam kap;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------

load myparam_values.mat

set_param_value('beta',myparam(1));
set_param_value('om',myparam(2));
set_param_value('eta',myparam(3));
set_param_value('si',myparam(4));
set_param_value('phi_pi',myparam(5));
set_param_value('phi_x',myparam(6));

set_param_value('rho_a',myparam(7));
set_param_value('rho_v',myparam(8));
set_param_value('rho_u',myparam(9));

%composite parameters
lam=(1+eta)/(si+eta);
kap=((eta+si)*(1-om)*(1-om*beta))/(om);

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model(linear); %we write linear since the model is already log linearize
%dynamic IS
x=x(+1)-(1/si)*(int-pi(+1))+lam*(a(+1)-a);
%NK Philips curve
pi=beta*pi(+1)+kap*x+u;
%Monethary policy rule
int=phi_x*x+phi_pi*pi+v;

%tree autoregressive process
a=rho_a*a(-1)+ea;
u=rho_u*u(-1)+eu;
v=rho_v*v(-1)+ev;

%additional varibale of interest
r=int-pi(+1);
yn=lam*a;
y=x+yn;
n=y-a;

end;

%since we have already provide log linearized model we can avoid to  declare
%steady state variable since are automatically 0

steady;

check; %check Blanchard-Kan conditions

%shock block
shocks;
var ea; stderr 1;
var eu; stderr 1;
var ev; stderr 1;
end;

@#include "util.txt" 
stoch_simul(periods=@{T},noprint,nofunctions,nomoments,nocorr,irf=0);
save("oo_");