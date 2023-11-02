function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 38);

T = DSGE.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(20) = getPowerDeriv(y(6),params(6),1);
T(21) = (-y(34))/(y(6)*y(6));
T(22) = getPowerDeriv(y(34)/y(6),(-params(6)),1);
T(23) = getPowerDeriv(y(6),(-params(6)),1);
T(24) = y(15)*T(23);
T(25) = y(13)*T(24);
T(26) = 1/params(2)*getPowerDeriv(y(8),params(8),1);
T(27) = params(12)*getPowerDeriv(y(8),params(11)/(1-params(1)),1);
T(28) = params(2)*params(12)*getPowerDeriv(y(35),T(16),1);
T(29) = params(2)*params(12)*getPowerDeriv(y(35),params(11)-1,1);
T(30) = getPowerDeriv(y(10),params(7),1);
T(31) = 1/y(16);
T(32) = getPowerDeriv(y(10)/y(16),1-params(1),1);
T(33) = y(16)*(-(y(13)*(1-params(1))))/(y(10)*y(10));
T(34) = getPowerDeriv(y(11),params(10),1);
T(35) = (-(y(13)*T(34)));
T(36) = 1/(steady_state(8))*getPowerDeriv(y(13)/(steady_state(8)),params(9),1);
T(37) = (-(y(7)*y(16)*(1-params(1))/y(10)));
T(38) = 1/y(15);

end
