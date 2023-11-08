function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = NK_NL_DSGE.static_resid_tt(T, y, x, params);
end
residual = zeros(28, 1);
lhs = y(2);
rhs = T(1)*T(2);
residual(1) = lhs - rhs;
lhs = y(9);
rhs = params(2)/y(3);
residual(2) = lhs - rhs;
lhs = y(6);
rhs = 1/y(9);
residual(3) = lhs - rhs;
lhs = y(8);
rhs = y(4)*T(3);
residual(4) = lhs - rhs;
lhs = y(6);
rhs = y(3)*y(7);
residual(5) = lhs - rhs;
lhs = y(6);
rhs = T(4)*T(5)*exp(y(28));
residual(6) = lhs - rhs;
lhs = y(1);
rhs = y(8);
residual(7) = lhs - rhs;
lhs = log(y(4));
rhs = log(y(4))*params(3)+x(1);
residual(8) = lhs - rhs;
lhs = log(y(10));
rhs = log(y(10))*params(5)-x(2);
residual(9) = lhs - rhs;
lhs = y(28);
rhs = y(28)*params(4)+x(3);
residual(10) = lhs - rhs;
lhs = y(15);
rhs = y(2)/T(6);
residual(11) = lhs - rhs;
lhs = 1;
rhs = params(12)*T(7)+(1-params(12))*y(12)^(1-params(11));
residual(12) = lhs - rhs;
lhs = y(11);
rhs = (1-params(12))*y(12)^((-params(11))/(1-params(1)))+y(11)*T(8);
residual(13) = lhs - rhs;
lhs = y(12)^(1+params(11)*params(1)/(1-params(1)));
rhs = params(11)*y(13)/y(14)*(1-params(13))/(params(11)-1);
residual(14) = lhs - rhs;
lhs = y(13);
rhs = y(15)*T(11)+y(13)*T(12);
residual(15) = lhs - rhs;
lhs = y(14);
rhs = T(11)+y(14)*T(7)*params(2)*params(12);
residual(16) = lhs - rhs;
lhs = y(22);
rhs = log(y(8));
residual(17) = lhs - rhs;
lhs = y(23);
rhs = log(y(2));
residual(18) = lhs - rhs;
lhs = y(24);
rhs = log(y(5));
residual(19) = lhs - rhs;
lhs = y(18);
rhs = 4*log(y(3));
residual(20) = lhs - rhs;
lhs = y(17);
rhs = 4*log(y(6));
residual(21) = lhs - rhs;
lhs = y(19);
rhs = 4*log(y(7));
residual(22) = lhs - rhs;
lhs = y(16);
rhs = y(8)/T(13);
residual(23) = lhs - rhs;
lhs = y(21);
rhs = log(y(16)*y(20));
residual(24) = lhs - rhs;
lhs = y(3);
rhs = 1;
residual(25) = lhs - rhs;
lhs = y(25);
rhs = log(y(20));
residual(26) = lhs - rhs;
lhs = y(26);
rhs = log(y(4));
residual(27) = lhs - rhs;
lhs = y(27);
rhs = log(y(10));
residual(28) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
