function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
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
%   g1
%

if T_flag
    T = DSGE.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(28, 41);
g1(1,6)=(-(T(2)*T(20)));
g1(1,7)=1;
g1(1,10)=(-(T(1)*T(30)));
g1(2,6)=(-(T(4)*params(2)*T(21)*T(22)/y(35)));
g1(2,34)=(-(T(4)*params(2)*T(22)*1/y(6)/y(35)));
g1(2,35)=(-((-T(5))/(y(35)*y(35))));
g1(2,14)=1;
g1(2,15)=(-(T(3)*(-y(36))/(y(15)*y(15))/y(35)));
g1(2,36)=(-(T(3)*T(38)/y(35)));
g1(3,11)=1;
g1(3,14)=(-((-1)/(y(14)*y(14))));
g1(4,9)=(-T(6));
g1(4,10)=(-(y(9)*T(31)*T(32)));
g1(4,13)=1;
g1(4,16)=(-(y(9)*T(32)*(-y(10))/(y(16)*y(16))));
g1(5,35)=(-y(12));
g1(5,11)=1;
g1(5,12)=(-y(35));
g1(6,8)=(-(exp(y(33))*T(8)*T(26)));
g1(6,11)=1;
g1(6,13)=(-(exp(y(33))*T(7)*T(36)));
g1(6,33)=(-(T(7)*T(8)*exp(y(33))));
g1(7,6)=1;
g1(7,13)=(-1);
g1(8,1)=(-(params(3)*1/y(1)));
g1(8,9)=1/y(9);
g1(8,39)=(-1);
g1(9,2)=(-(params(5)*1/y(2)));
g1(9,15)=T(38);
g1(9,40)=1;
g1(10,5)=(-params(4));
g1(10,33)=1;
g1(10,41)=(-1);
g1(11,7)=(-(1/T(10)));
g1(11,10)=(-((-(y(7)*T(33)))/(T(10)*T(10))));
g1(11,13)=(-(T(37)/(T(10)*T(10))));
g1(11,16)=(-((-(y(7)*T(9)))/(T(10)*T(10))));
g1(11,20)=1;
g1(12,8)=(-(params(12)*getPowerDeriv(y(8),params(11)-1,1)));
g1(12,17)=(-((1-params(12))*getPowerDeriv(y(17),1-params(11),1)));
g1(13,8)=(-(y(3)*T(27)));
g1(13,3)=(-T(11));
g1(13,16)=1;
g1(13,17)=(-((1-params(12))*getPowerDeriv(y(17),(-params(11))/(1-params(1)),1)));
g1(14,17)=getPowerDeriv(y(17),T(12),1);
g1(14,18)=(-(params(11)*(1-params(13))*1/y(19)/(params(11)-1)));
g1(14,19)=(-(params(11)*(1-params(13))*(-y(18))/(y(19)*y(19))/(params(11)-1)));
g1(15,6)=(-(y(20)*T(25)));
g1(15,35)=(-(y(37)*T(28)));
g1(15,13)=(-(y(20)*T(14)));
g1(15,15)=(-(y(20)*y(13)*T(13)));
g1(15,18)=1;
g1(15,37)=(-T(17));
g1(15,20)=(-T(15));
g1(16,6)=(-T(25));
g1(16,35)=(-(y(38)*T(29)));
g1(16,13)=(-T(14));
g1(16,15)=(-(y(13)*T(13)));
g1(16,19)=1;
g1(16,38)=(-T(18));
g1(17,13)=(-(1/y(13)));
g1(17,27)=1;
g1(18,7)=(-(1/y(7)));
g1(18,28)=1;
g1(19,10)=(-(1/y(10)));
g1(19,29)=1;
g1(20,8)=(-(4*1/y(8)));
g1(20,23)=1;
g1(21,11)=(-(4*1/y(11)));
g1(21,22)=1;
g1(22,12)=(-(4*1/y(12)));
g1(22,24)=1;
g1(23,11)=(-(T(35)/(T(19)*T(19))));
g1(23,13)=(-(1/T(19)));
g1(23,21)=1;
g1(24,21)=(-(y(25)/(y(21)*y(25))));
g1(24,25)=(-(y(21)/(y(21)*y(25))));
g1(24,26)=1;
g1(25,8)=1;
g1(25,4)=(-((-y(25))/(y(4)*y(4))));
g1(25,25)=(-(1/y(4)));
g1(26,25)=(-(1/y(25)));
g1(26,30)=1;
g1(27,9)=(-(1/y(9)));
g1(27,31)=1;
g1(28,15)=(-T(38));
g1(28,32)=1;

end
