function T = dynamic_g2_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g2_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 79);

T = NK_NL_DSGE.dynamic_g1_tt(T, y, x, params, steady_state, it_);

T(52) = getPowerDeriv(y(6),params(6),2);
T(53) = getPowerDeriv(y(10),params(7),2);
T(54) = (-((-y(34))*(y(6)+y(6))))/(y(6)*y(6)*y(6)*y(6));
T(55) = getPowerDeriv(T(3),(-params(6)),2);
T(56) = (-1)/(y(6)*y(6));
T(57) = (-((-y(36))*(y(15)+y(15))))/(y(15)*y(15)*y(15)*y(15));
T(58) = (-1)/(y(15)*y(15));
T(59) = getPowerDeriv(T(7),1-params(1),2);
T(60) = (-1)/(y(16)*y(16));
T(61) = T(39)*T(60)+T(38)*T(51)*T(59);
T(62) = (-((-y(10))*(y(16)+y(16))))/(y(16)*y(16)*y(16)*y(16));
T(63) = T(51)*T(51)*T(59)+T(39)*T(62);
T(64) = T(9)*getPowerDeriv(y(8),params(8),2);
T(65) = T(44)*T(44)*getPowerDeriv(T(11),params(9),2);
T(66) = (-T(47));
T(67) = y(16)*(-((-(y(13)*(1-params(1))))*(y(10)+y(10))))/(y(10)*y(10)*y(10)*y(10));
T(68) = T(14)*T(14)*(-(y(7)*T(67)))-(-(y(7)*T(41)))*(T(14)*T(41)+T(14)*T(41));
T(69) = T(14)*T(14)*T(14)*T(14);
T(70) = y(16)*(-(1-params(1)))/(y(10)*y(10));
T(71) = T(14)*T(47)+T(14)*T(47);
T(72) = params(11)*getPowerDeriv(y(8),T(16),2);
T(73) = getPowerDeriv(y(6),(-params(6)),2);
T(74) = y(15)*T(73);
T(75) = y(13)*T(74);
T(76) = params(2)*params(11)*getPowerDeriv(y(35),T(22),2);
T(77) = params(2)*params(11)*getPowerDeriv(y(35),params(10)-1,2);
T(78) = getPowerDeriv(y(11),params(13),2);
T(79) = (-(y(13)*T(78)));

end
