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

assert(length(T) >= 19);

T(1) = y(6)^params(6);
T(2) = y(10)^params(7);
T(3) = params(2)*(y(34)/y(6))^(-params(6));
T(4) = y(36)/y(15);
T(5) = T(3)*T(4);
T(6) = (y(10)/y(16))^(1-params(1));
T(7) = 1/params(2)*y(8)^params(8);
T(8) = (y(13)/(steady_state(8)))^params(9);
T(9) = y(13)*(1-params(1))/y(10);
T(10) = y(16)*T(9);
T(11) = params(12)*y(8)^(params(11)/(1-params(1)));
T(12) = 1+params(11)*params(1)/(1-params(1));
T(13) = y(6)^(-params(6));
T(14) = y(15)*T(13);
T(15) = y(13)*T(14);
T(16) = params(11)+params(1)*params(11)/(1-params(1));
T(17) = params(2)*params(12)*y(35)^T(16);
T(18) = params(2)*params(12)*y(35)^(params(11)-1);
T(19) = y(11)^params(10);

end
