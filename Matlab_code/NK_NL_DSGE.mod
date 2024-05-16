/*
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * It is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * For a copy of the GNU General Public License,
 * see <http://www.gnu.org/licenses/>.
 */
load my_seed.mat
set_dynare_seed(my_seed)
%set_dynare_seed('clock');
%set_dynare_seed(1234);
%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------
var 
    C               ${C}$           (long_name='Consumption')
    W_real          ${\frac{W}{P}}$ (long_name='Real Wage')
    Pi              ${\Pi}$         (long_name='inflation')
    A               ${A}$           (long_name='AR(1) technology process')
    N               ${N}$           (long_name='Hours worked')
    R               ${R^n}$         (long_name='Nominal Interest Rate') 
    realinterest    ${R^{r}}$       (long_name='Real Interest Rate')
    Y               ${Y}$           (long_name='Output') 
    Q               ${Q}$           (long_name='Bond price')
    Z               ${Z}$           (long_name='AR(1) preference shock process')
    S               ${S}$           (long_name='Price dispersion')
    Pi_star         ${\Pi^*}$       (long_name='Optimal reset price')
    x_aux_1         ${x_1}$         (long_name='aux. var. 1 recursive price setting')
    x_aux_2         ${x_2}$         (long_name='aux. var. 2 recursive price setting')
    MC              ${mc}$          (long_name='real marginal costs')
    P               ${P}$           (long_name='price level')
    nu              ${\nu}$         (long_name='AR(1) monetary policy shock process')    
    ;

varexo 
    eps_a        ${\varepsilon_a}$   (long_name='technology shock')
    eps_z        ${\varepsilon_z}$   (long_name='preference shock')
    eps_nu       ${\varepsilon_\nu}$ (long_name='monetary policy shock')
    ;

parameters 
    alppha              ${\alpha}$      (long_name='capital share')
    betta               ${\beta}$       (long_name='discount factor')
    rho_a               ${\rho_a}$      (long_name='autocorrelation technology shock')
    rho_nu              ${\rho_{\nu}}$  (long_name='autocorrelation monetary policy shock')
    rho_z               ${\rho_{z}}$    (long_name='autocorrelation monetary demand shock')
    siggma              ${\sigma}$      (long_name='inverse EIS')
    varphi              ${\varphi}$     (long_name='inverse Frisch elasticity')
    phi_pi              ${\phi_{\pi}}$  (long_name='inflation feedback Taylor Rule')
    phi_y               ${\phi_{y}}$    (long_name='output feedback Taylor Rule')
    epsilon             ${\epsilon}$    (long_name='demand elasticity')
    theta               ${\theta}$      (long_name='Calvo parameter')
    eta                 ${\eta}$        (long_name='semi-elasticity of money demand')
    ;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------
load myparam_values.mat

set_param_value('alppha',params_t(1));
set_param_value('betta',params_t(2));
set_param_value('rho_a',params_t(3));
set_param_value('rho_nu',params_t(4));
set_param_value('rho_z',params_t(5));
set_param_value('siggma',params_t(6));
set_param_value('varphi',params_t(7));
set_param_value('phi_pi',params_t(8));
set_param_value('phi_y',params_t(9));
set_param_value('epsilon',params_t(10));
set_param_value('theta',params_t(11));
set_param_value('eta',params_t(12));

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------
model;
    %[name='Preference Shock, p.54']
    log(Z)=rho_z*log(Z(-1))-eps_z;
    %[name='Definition nominal interest rate), p. 22 top']
    R=1/Q;
    %[name='Definition Real interest rate']
    R=realinterest*Pi(+1);
    %[name='Definition price level']
    Pi=P/P(-1);
    %[name='FOC Wages, eq. (2)']
    W_real=C^siggma*N^varphi;
    %[name='Euler equation eq. (3)']
    Q=betta*(C(+1)/C)^(-siggma)*(Z(+1)/Z)/Pi(+1);
    %[name='Technology Shock, eq. (6)']
    log(A)=rho_a*log(A(-1))+eps_a;
    %[name='LOM prices, eq. (7)']
    1=theta*Pi^(epsilon-1)+(1-theta)*(Pi_star)^(1-epsilon);
    %[name='Market Clearing, eq. (15)']
    C=Y;
    %[name='Aggregate output, above eq. (14)']
    Y=A*(N/S)^(1-alppha);
    %[name='LOM price dispersion']
    S=(1-theta)*Pi_star^(-epsilon/(1-alppha))+theta*Pi^(epsilon/(1-alppha))*S(-1);
    %[name='Definition marginal cost']
    MC=W_real/((1-alppha)*Y/N*S);
    %[name='FOC price setting']
    Pi_star^(1+epsilon*(alppha/(1-alppha)))=x_aux_1/x_aux_2*epsilon/(epsilon-1);
    %[name='Auxiliary price setting recursion 2']
    x_aux_2=Z*C^(-siggma)*Y+betta*theta*Pi(+1)^(epsilon-1)*x_aux_2(+1);
    %[name='Auxiliary price setting recursion 1']
    x_aux_1=Z*C^(-siggma)*Y*MC+betta*theta*Pi(+1)^(epsilon+alppha*epsilon/(1-alppha))*x_aux_1(+1);
    %[name='Monetary Policy Rule, p. 26 bottom/eq. (22)']
    R=1/betta*Pi^phi_pi*(Y/steady_state(Y))^phi_y*exp(nu);
    %[name='Monetary policy shock']
    nu=rho_nu*nu(-1)+eps_nu;
end;

%----------------------------------------------------------------
% 4. Steady state values
%---------------------------------------------------------------

steady_state_model;
    A=1;
    Z=1;
    S=1;
    Pi_star=1;
    P=1;
    MC=(epsilon-1)/epsilon;
    R=1/betta;
    Pi=1;
    Q=1/R;
    realinterest=R;
    N=((1-alppha)*MC)^(1/((1-siggma)*alppha+varphi+siggma));
    C=A*N^(1-alppha);
    W_real=C^siggma*N^varphi;
    Y=C;
    nu=0;
    x_aux_1=C^(-siggma)*Y*MC/(1-betta*theta*Pi^(epsilon/(1-alppha)));
    x_aux_2=C^(-siggma)*Y/(1-betta*theta*Pi^(epsilon-1));
end;



resid;
steady;
check;

%----------------------------------------------------------------
% 5. shocks
%---------------------------------------------------------------
shocks;
    var eps_a  = 0.00025;
    var eps_nu = 0.00025; 
    var eps_z = 0.00025; 

end;

@#include "util.txt" 
stoch_simul(periods=@{T},order=3,pruning, noprint,nofunctions,nomoments,nocorr,irf=0);

save("oo_");