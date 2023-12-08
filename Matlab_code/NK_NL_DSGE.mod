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
 
set_dynare_seed('clock');

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
    M_real          ${M/P}$         (long_name='real money stock')
    i_ann           ${i^{ann}}$     (long_name='annualized nominal interest rate')
    pi_ann          ${\pi^{ann}}$   (long_name='annualized inflation rate')
    r_real_ann      ${r^{r,ann}}$   (long_name='annualized real interest rate')
    P               ${P}$           (long_name='price level')
    log_m_nominal   ${log(M)}$      (long_name='log nominal money stock')
    log_y           ${log(Y)}$      (long_name='log output')
    log_W_real      ${log(W/P)}$    (long_name='log real wage')
    log_N           ${log(N)}$      (long_name='log hours')
    log_P           ${log(P)}$      (long_name='log price level')
    log_A           ${log(A)}$      (long_name='log technology level')
    log_Z           ${log(Z)}$      (long_name='log preference shock')
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
    tau                 ${\tau}$      (long_name='labor subsidy')
    eta                 ${\eta}$        (long_name='semi-elasticity of money demand')
    ;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------


%siggma = 1;
%varphi=5;
%phi_pi = 1.5;
%phi_y  = 0.125;
%theta=3/4;
%rho_nu =0.5;
%rho_z  = 0.5;
%rho_a  = 0.9;
%betta  = 0.99;
%eta  =3.77; %footnote 11, p. 115
%alppha=1/4;
%epsilon=9;
%tau=0; //1/epsilon;

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
set_param_value('tau',params_t(12));
set_param_value('eta',params_t(13));

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------
model;
    [name='FOC Wages, eq. (2)']
    W_real=C^siggma*N^varphi;
    [name='Euler equation eq. (3)']
    Q=betta*(C(+1)/C)^(-siggma)*(Z(+1)/Z)/Pi(+1);
    [name='Definition nominal interest rate), p. 22 top']
    R=1/Q;
    [name='Aggregate output, above eq. (14)']
    Y=A*(N/S)^(1-alppha);
    [name='Definition Real interest rate']
    R=realinterest*Pi(+1);
    [name='Monetary Policy Rule, p. 26 bottom/eq. (22)']
    R=1/betta*Pi^phi_pi*(Y/steady_state(Y))^phi_y*exp(nu);
    [name='Market Clearing, eq. (15)']
    C=Y;
    [name='Technology Shock, eq. (6)']
    log(A)=rho_a*log(A(-1))+eps_a;
    [name='Preference Shock, p.54']
    log(Z)=rho_z*log(Z(-1))-eps_z;
    [name='Monetary policy shock']
    nu=rho_nu*nu(-1)+eps_nu;
    [name='Definition marginal cost']
    MC=W_real/((1-alppha)*Y/N*S);
    [name='LOM prices, eq. (7)']
    1=theta*Pi^(epsilon-1)+(1-theta)*(Pi_star)^(1-epsilon);
    [name='LOM price dispersion']
    S=(1-theta)*Pi_star^(-epsilon/(1-alppha))+theta*Pi^(epsilon/(1-alppha))*S(-1);
    [name='FOC price setting']
    Pi_star^(1+epsilon*(alppha/(1-alppha)))=x_aux_1/x_aux_2*(1-tau)*epsilon/(epsilon-1);
    [name='Auxiliary price setting recursion 1']
    x_aux_1=Z*C^(-siggma)*Y*MC+betta*theta*Pi(+1)^(epsilon+alppha*epsilon/(1-alppha))*x_aux_1(+1);
    [name='Auxiliary price setting recursion 2']
    x_aux_2=Z*C^(-siggma)*Y+betta*theta*Pi(+1)^(epsilon-1)*x_aux_2(+1);

    %obs (maybe)
    [name='Definition log output']
    log_y = log(Y);
    [name='Definition log real wage']
    log_W_real=log(W_real);
    [name='Definition log hours']
    log_N=log(N);
    [name='Annualized inflation']
    pi_ann=4*log(Pi);
    [name='Annualized nominal interest rate']
    i_ann=4*log(R);
    [name='Annualized real interest rate']
    r_real_ann=4*log(realinterest);
    [name='Real money demand, eq. (4)']
    M_real=Y/R^eta;
    [name='definition nominal money stock']
    log_m_nominal=log(M_real*P);
    [name='Definition price level']
    Pi=P/P(-1);
    [name='Definition log price level']
    log_P=log(P);
    [name='Definition log TFP']
    log_A=log(A);
    [name='Definition log preference']
    log_Z=log(Z);

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
MC=(epsilon-1)/epsilon/(1-tau);
R=1/betta;
Pi=1;
Q=1/R;
realinterest=R;
N=((1-alppha)*MC)^(1/((1-siggma)*alppha+varphi+siggma));
C=A*N^(1-alppha);
W_real=C^siggma*N^varphi;
Y=C;
money_growth=0;
money_growth_ann=0;
nu=0;
x_aux_1=C^(-siggma)*Y*MC/(1-betta*theta*Pi^(epsilon/(1-alppha)));
x_aux_2=C^(-siggma)*Y/(1-betta*theta*Pi^(epsilon-1));
log_y = log(Y);
log_W_real=log(W_real);
log_N=log(N);
pi_ann=4*log(Pi);
i_ann=4*log(R);
r_real_ann=4*log(realinterest);
M_real=Y/R^eta;
log_m_nominal=log(M_real*P);
log_P=log(P);
log_A=0;
log_Z=0;
end;

%//write_latex_dynamic_model;

resid;
steady;
check;

%----------------------------------------------------------------
% 5a. Positive technological shock
%---------------------------------------------------------------
shocks;
    var eps_a  = 0.001^2; //unit shock to technology
    var eps_nu = 0.0025^2; %//1 standard deviation shock of 25 basis points, i.e. 1 percentage point annualized
    var eps_z = 0.001^2; 

end;

%----------------------------------------------------------------
% 5b. Monetary policy shock
%---------------------------------------------------------------
%shocks;
%    var eps_nu = 0.0025^2; %//1 standard deviation shock of 25 basis points, i.e. 1 percentage point annualized
%end;

%stoch_simul(order = 1, nograph, nomoments, irf=16);

%add by me
@#include "util.txt" 


stoch_simul(periods=@{T},order=3, noprint,nofunctions,nomoments,nocorr,irf=0);

save("oo_");