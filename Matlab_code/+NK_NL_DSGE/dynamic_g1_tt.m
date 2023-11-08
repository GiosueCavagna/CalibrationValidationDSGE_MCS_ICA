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

assert(length(T) >= 51);

T = NK_NL_DSGE.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(26) = getPowerDeriv(y(6),params(6),1);
T(27) = (-y(34))/(y(6)*y(6));
T(28) = getPowerDeriv(T(3),(-params(6)),1);
T(29) = getPowerDeriv(y(6),(-params(6)),1);
T(30) = y(15)*T(29);
T(31) = y(13)*T(30);
T(32) = 1/y(6);
T(33) = T(9)*getPowerDeriv(y(8),params(8),1);
T(34) = params(12)*getPowerDeriv(y(8),T(16),1);
T(35) = params(2)*params(12)*getPowerDeriv(y(35),T(22),1);
T(36) = params(2)*params(12)*getPowerDeriv(y(35),params(11)-1,1);
T(37) = getPowerDeriv(y(10),params(7),1);
T(38) = 1/y(16);
T(39) = getPowerDeriv(T(7),1-params(1),1);
T(40) = (-(y(13)*(1-params(1))))/(y(10)*y(10));
T(41) = y(16)*T(40);
T(42) = getPowerDeriv(y(11),params(10),1);
T(43) = (-(y(13)*T(42)));
T(44) = 1/(steady_state(8));
T(45) = T(44)*getPowerDeriv(T(11),params(9),1);
T(46) = (1-params(1))/y(10);
T(47) = y(16)*T(46);
T(48) = (-(y(7)*T(47)));
T(49) = (-y(36))/(y(15)*y(15));
T(50) = 1/y(15);
T(51) = (-y(10))/(y(16)*y(16));

end
