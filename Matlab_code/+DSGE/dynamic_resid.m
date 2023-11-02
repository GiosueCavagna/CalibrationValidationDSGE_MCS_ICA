function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
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
%   residual
%

if T_flag
    T = DSGE.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(28, 1);
lhs = y(7);
rhs = T(1)*T(2);
residual(1) = lhs - rhs;
lhs = y(14);
rhs = T(5)/y(35);
residual(2) = lhs - rhs;
lhs = y(11);
rhs = 1/y(14);
residual(3) = lhs - rhs;
lhs = y(13);
rhs = y(9)*T(6);
residual(4) = lhs - rhs;
lhs = y(11);
rhs = y(35)*y(12);
residual(5) = lhs - rhs;
lhs = y(11);
rhs = T(7)*T(8)*exp(y(33));
residual(6) = lhs - rhs;
lhs = y(6);
rhs = y(13);
residual(7) = lhs - rhs;
lhs = log(y(9));
rhs = params(3)*log(y(1))+x(it_, 1);
residual(8) = lhs - rhs;
lhs = log(y(15));
rhs = params(5)*log(y(2))-x(it_, 2);
residual(9) = lhs - rhs;
lhs = y(33);
rhs = params(4)*y(5)+x(it_, 3);
residual(10) = lhs - rhs;
lhs = y(20);
rhs = y(7)/T(10);
residual(11) = lhs - rhs;
lhs = 1;
rhs = params(12)*y(8)^(params(11)-1)+(1-params(12))*y(17)^(1-params(11));
residual(12) = lhs - rhs;
lhs = y(16);
rhs = (1-params(12))*y(17)^((-params(11))/(1-params(1)))+T(11)*y(3);
residual(13) = lhs - rhs;
lhs = y(17)^T(12);
rhs = params(11)*y(18)/y(19)*(1-params(13))/(params(11)-1);
residual(14) = lhs - rhs;
lhs = y(18);
rhs = y(20)*T(15)+T(17)*y(37);
residual(15) = lhs - rhs;
lhs = y(19);
rhs = T(15)+T(18)*y(38);
residual(16) = lhs - rhs;
lhs = y(27);
rhs = log(y(13));
residual(17) = lhs - rhs;
lhs = y(28);
rhs = log(y(7));
residual(18) = lhs - rhs;
lhs = y(29);
rhs = log(y(10));
residual(19) = lhs - rhs;
lhs = y(23);
rhs = 4*log(y(8));
residual(20) = lhs - rhs;
lhs = y(22);
rhs = 4*log(y(11));
residual(21) = lhs - rhs;
lhs = y(24);
rhs = 4*log(y(12));
residual(22) = lhs - rhs;
lhs = y(21);
rhs = y(13)/T(19);
residual(23) = lhs - rhs;
lhs = y(26);
rhs = log(y(21)*y(25));
residual(24) = lhs - rhs;
lhs = y(8);
rhs = y(25)/y(4);
residual(25) = lhs - rhs;
lhs = y(30);
rhs = log(y(25));
residual(26) = lhs - rhs;
lhs = y(31);
rhs = log(y(9));
residual(27) = lhs - rhs;
lhs = y(32);
rhs = log(y(15));
residual(28) = lhs - rhs;

end
