%----------------------------------------------------------------
% 1. Labeling
%----------------------------------------------------------------
var c k a;
varexo ea;
parameters alpha beta sigma delta rhoa;

%----------------------------------------------------------------
% 2. Parameters
%----------------------------------------------------------------

% the following two lines show how to load the value of a parameter (alpha) from another file (myparam_values.mat)
load myparam_values.mat
set_param_value('alpha',myparam(1));
set_param_value('beta',myparam(2));
set_param_value('sigma',myparam(3));
set_param_value('delta',myparam(4));
set_param_value('rhoa',myparam(5));

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------
model;
    c^(-sigma)=beta*c(+1)^(-sigma)*(alpha*exp(a(+1))*k^(alpha-1)+1-delta);
    c+k=exp(a)*k(-1)^alpha+(1-delta)*k(-1); 
    a = rhoa*a(-1)+ ea;    
end;

%model;
%    exp(c)^(-sigma)=beta*exp(c(+1))^(-sigma)*(alpha*exp(a(+1))*exp(k)^(alpha-1)+1-delta);
%    exp(c)+exp(k)=exp(a)*exp(k(-1))^alpha+(1-delta)*exp(k(-1)); 
%    a = rhoa*a(-1)+ ea;
%end;


%----------------------------------------------------------------
% 4. Solution and Simulation
%----------------------------------------------------------------

%steady state
initval; % initial guess for the steay state
    k= exp(log(0.5));
    c = exp(log(-exp(k)+exp(k)^alpha+(1-delta)*exp(k))); 
    a = 1;
end;

steady;


% % solution with rootfinding (no uncertainty)
% shocks;
% var ea; periods 3:5 6; values 0.01 0.005;
% end;  
% 
% perfect_foresight_setup(periods=200);
% perfect_foresight_solver;
% rplot k ;

% solution with perturbation
shocks;
var ea; periods 3:5 6; values 0.01 0.005;
end;

@#include "util.txt" 

stoch_simul(periods=@{T},noprint,nofunctions,nomoments,nocorr,irf=0);
save("oo_");

%stoch_simul(order = 2) c k;
%stoch_simul(order = 2, pruning) c k;