function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double  vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double  vector of endogenous variables in the order stored
%                                                    in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double  matrix of exogenous variables (in declaration order)
%                                                    for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double  vector of steady state values
%   params        [M_.param_nbr by 1]        double  vector of parameter values in declaration order
%   it_           scalar                     double  time period for exogenous variables for which
%                                                    to evaluate the model
%
% Output:
%   T           [#temp variables by 1]       double  vector of temporary terms
%

assert(length(T) >= 25);

T(1) = y(6)^params(6);
T(2) = y(10)^params(7);
T(3) = y(34)/y(6);
T(4) = params(2)*T(3)^(-params(6));
T(5) = y(36)/y(15);
T(6) = T(4)*T(5);
T(7) = y(10)/y(16);
T(8) = T(7)^(1-params(1));
T(9) = 1/params(2);
T(10) = T(9)*y(8)^params(8);
T(11) = y(13)/(steady_state(8));
T(12) = T(11)^params(9);
T(13) = y(13)*(1-params(1))/y(10);
T(14) = y(16)*T(13);
T(15) = (-params(10))/(1-params(1));
T(16) = params(10)/(1-params(1));
T(17) = params(11)*y(8)^T(16);
T(18) = 1+params(10)*params(1)/(1-params(1));
T(19) = y(6)^(-params(6));
T(20) = y(15)*T(19);
T(21) = y(13)*T(20);
T(22) = params(10)+params(1)*params(10)/(1-params(1));
T(23) = params(2)*params(11)*y(35)^T(22);
T(24) = params(2)*params(11)*y(35)^(params(10)-1);
T(25) = y(11)^params(13);

end
