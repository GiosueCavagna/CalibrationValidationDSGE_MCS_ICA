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
    T = NK_DSGE.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(10, 1);
lhs = y(4);
rhs = y(14)-1/params(4)*(y(6)-y(15))+params(10)*(y(16)-y(7));
residual(1) = lhs - rhs;
lhs = y(5);
rhs = y(15)*params(1)+y(4)*params(11)+y(8);
residual(2) = lhs - rhs;
lhs = y(6);
rhs = y(4)*params(6)+y(5)*params(5)+y(9);
residual(3) = lhs - rhs;
lhs = y(7);
rhs = params(7)*y(1)+x(it_, 1);
residual(4) = lhs - rhs;
lhs = y(8);
rhs = params(9)*y(2)+x(it_, 2);
residual(5) = lhs - rhs;
lhs = y(9);
rhs = params(8)*y(3)+x(it_, 3);
residual(6) = lhs - rhs;
lhs = y(10);
rhs = y(6)-y(15);
residual(7) = lhs - rhs;
lhs = y(12);
rhs = params(10)*y(7);
residual(8) = lhs - rhs;
lhs = y(11);
rhs = y(4)+y(12);
residual(9) = lhs - rhs;
lhs = y(13);
rhs = y(11)-y(7);
residual(10) = lhs - rhs;

end
