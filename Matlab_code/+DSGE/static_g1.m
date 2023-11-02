function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
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
%   g1
%

if T_flag
    T = DSGE.static_g1_tt(T, y, x, params);
end
g1 = zeros(28, 28);
g1(1,1)=(-(T(2)*getPowerDeriv(y(1),params(6),1)));
g1(1,2)=1;
g1(1,5)=(-(T(1)*getPowerDeriv(y(5),params(7),1)));
g1(2,3)=(-((-params(2))/(y(3)*y(3))));
g1(2,9)=1;
g1(3,6)=1;
g1(3,9)=(-((-1)/(y(9)*y(9))));
g1(4,4)=(-T(3));
g1(4,5)=(-(y(4)*1/y(11)*T(16)));
g1(4,8)=1;
g1(4,11)=(-(y(4)*T(16)*(-y(5))/(y(11)*y(11))));
g1(5,3)=(-y(7));
g1(5,6)=1;
g1(5,7)=(-y(3));
g1(6,3)=(-(exp(y(28))*T(5)*1/params(2)*getPowerDeriv(y(3),params(8),1)));
g1(6,6)=1;
g1(6,8)=(-(exp(y(28))*T(4)*((y(8))-y(8))/((y(8))*(y(8)))*getPowerDeriv(y(8)/(y(8)),params(9),1)));
g1(6,28)=(-(T(4)*T(5)*exp(y(28))));
g1(7,1)=1;
g1(7,8)=(-1);
g1(8,4)=1/y(4)-params(3)*1/y(4);
g1(9,10)=1/y(10)-params(5)*1/y(10);
g1(10,28)=1-params(4);
g1(11,2)=(-(1/T(6)));
g1(11,5)=(-((-(y(2)*y(11)*(-(y(8)*(1-params(1))))/(y(5)*y(5))))/(T(6)*T(6))));
g1(11,8)=(-((-(y(2)*y(11)*(1-params(1))/y(5)))/(T(6)*T(6))));
g1(11,11)=(-((-(y(2)*y(8)*(1-params(1))/y(5)))/(T(6)*T(6))));
g1(11,15)=1;
g1(12,3)=(-(params(12)*T(15)));
g1(12,12)=(-((1-params(12))*getPowerDeriv(y(12),1-params(11),1)));
g1(13,3)=(-(y(11)*params(12)*getPowerDeriv(y(3),params(11)/(1-params(1)),1)));
g1(13,11)=1-T(8);
g1(13,12)=(-((1-params(12))*getPowerDeriv(y(12),(-params(11))/(1-params(1)),1)));
g1(14,12)=getPowerDeriv(y(12),1+params(11)*params(1)/(1-params(1)),1);
g1(14,13)=(-(params(11)*(1-params(13))*1/y(14)/(params(11)-1)));
g1(14,14)=(-(params(11)*(1-params(13))*(-y(13))/(y(14)*y(14))/(params(11)-1)));
g1(15,1)=(-(y(15)*T(14)));
g1(15,3)=(-(y(13)*params(2)*params(12)*getPowerDeriv(y(3),params(11)+params(1)*params(11)/(1-params(1)),1)));
g1(15,8)=(-(y(15)*T(10)));
g1(15,10)=(-(y(15)*y(8)*T(9)));
g1(15,13)=1-T(12);
g1(15,15)=(-T(11));
g1(16,1)=(-T(14));
g1(16,3)=(-(y(14)*params(2)*params(12)*T(15)));
g1(16,8)=(-T(10));
g1(16,10)=(-(y(8)*T(9)));
g1(16,14)=1-T(7)*params(2)*params(12);
g1(17,8)=(-(1/y(8)));
g1(17,22)=1;
g1(18,2)=(-(1/y(2)));
g1(18,23)=1;
g1(19,5)=(-(1/y(5)));
g1(19,24)=1;
g1(20,3)=(-(4*1/y(3)));
g1(20,18)=1;
g1(21,6)=(-(4*1/y(6)));
g1(21,17)=1;
g1(22,7)=(-(4*1/y(7)));
g1(22,19)=1;
g1(23,6)=(-((-(y(8)*getPowerDeriv(y(6),params(10),1)))/(T(13)*T(13))));
g1(23,8)=(-(1/T(13)));
g1(23,16)=1;
g1(24,16)=(-(y(20)/(y(16)*y(20))));
g1(24,20)=(-(y(16)/(y(16)*y(20))));
g1(24,21)=1;
g1(25,3)=1;
g1(26,20)=(-(1/y(20)));
g1(26,25)=1;
g1(27,4)=(-(1/y(4)));
g1(27,26)=1;
g1(28,10)=(-(1/y(10)));
g1(28,27)=1;
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
