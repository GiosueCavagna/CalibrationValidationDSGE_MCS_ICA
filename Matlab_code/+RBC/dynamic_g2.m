function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
% function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g2
%

if T_flag
    T = RBC.dynamic_g2_tt(T, y, x, params, steady_state, it_);
end
g2_i = zeros(14,1);
g2_j = zeros(14,1);
g2_v = zeros(14,1);

g2_i(1)=1;
g2_i(2)=1;
g2_i(3)=1;
g2_i(4)=1;
g2_i(5)=1;
g2_i(6)=1;
g2_i(7)=1;
g2_i(8)=1;
g2_i(9)=1;
g2_i(10)=1;
g2_i(11)=2;
g2_i(12)=2;
g2_i(13)=2;
g2_i(14)=2;
g2_j(1)=19;
g2_j(2)=46;
g2_j(3)=44;
g2_j(4)=30;
g2_j(5)=47;
g2_j(6)=54;
g2_j(7)=28;
g2_j(8)=31;
g2_j(9)=52;
g2_j(10)=55;
g2_j(11)=1;
g2_j(12)=5;
g2_j(13)=33;
g2_j(14)=37;
g2_v(1)=getPowerDeriv(y(3),(-params(3)),2);
g2_v(2)=(-(T(3)*params(2)*getPowerDeriv(y(6),(-params(3)),2)));
g2_v(3)=(-(T(5)*T(7)));
g2_v(4)=g2_v(3);
g2_v(5)=(-(T(2)*T(5)));
g2_v(6)=g2_v(5);
g2_v(7)=(-(T(1)*params(1)*exp(y(7))*getPowerDeriv(y(4),params(1)-1,2)));
g2_v(8)=(-(T(1)*T(7)));
g2_v(9)=g2_v(8);
g2_v(10)=(-(T(1)*T(2)));
g2_v(11)=(-(exp(y(5))*getPowerDeriv(y(1),params(1),2)));
g2_v(12)=(-T(6));
g2_v(13)=g2_v(12);
g2_v(14)=(-T(4));
g2 = sparse(g2_i,g2_j,g2_v,3,64);
end
