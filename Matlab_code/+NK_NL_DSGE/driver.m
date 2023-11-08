%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'NK_NL_DSGE';
M_.dynare_version = '5.4';
oo_.dynare_version = '5.4';
options_.dynare_version = '5.4';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(3,1);
M_.exo_names_tex = cell(3,1);
M_.exo_names_long = cell(3,1);
M_.exo_names(1) = {'eps_a'};
M_.exo_names_tex(1) = {'{\varepsilon_a}'};
M_.exo_names_long(1) = {'technology shock'};
M_.exo_names(2) = {'eps_z'};
M_.exo_names_tex(2) = {'{\varepsilon_z}'};
M_.exo_names_long(2) = {'preference shock'};
M_.exo_names(3) = {'eps_nu'};
M_.exo_names_tex(3) = {'{\varepsilon_\nu}'};
M_.exo_names_long(3) = {'monetary policy shock'};
M_.endo_names = cell(28,1);
M_.endo_names_tex = cell(28,1);
M_.endo_names_long = cell(28,1);
M_.endo_names(1) = {'C'};
M_.endo_names_tex(1) = {'{C}'};
M_.endo_names_long(1) = {'Consumption'};
M_.endo_names(2) = {'W_real'};
M_.endo_names_tex(2) = {'{\frac{W}{P}}'};
M_.endo_names_long(2) = {'Real Wage'};
M_.endo_names(3) = {'Pi'};
M_.endo_names_tex(3) = {'{\Pi}'};
M_.endo_names_long(3) = {'inflation'};
M_.endo_names(4) = {'A'};
M_.endo_names_tex(4) = {'{A}'};
M_.endo_names_long(4) = {'AR(1) technology process'};
M_.endo_names(5) = {'N'};
M_.endo_names_tex(5) = {'{N}'};
M_.endo_names_long(5) = {'Hours worked'};
M_.endo_names(6) = {'R'};
M_.endo_names_tex(6) = {'{R^n}'};
M_.endo_names_long(6) = {'Nominal Interest Rate'};
M_.endo_names(7) = {'realinterest'};
M_.endo_names_tex(7) = {'{R^{r}}'};
M_.endo_names_long(7) = {'Real Interest Rate'};
M_.endo_names(8) = {'Y'};
M_.endo_names_tex(8) = {'{Y}'};
M_.endo_names_long(8) = {'Output'};
M_.endo_names(9) = {'Q'};
M_.endo_names_tex(9) = {'{Q}'};
M_.endo_names_long(9) = {'Bond price'};
M_.endo_names(10) = {'Z'};
M_.endo_names_tex(10) = {'{Z}'};
M_.endo_names_long(10) = {'AR(1) preference shock process'};
M_.endo_names(11) = {'S'};
M_.endo_names_tex(11) = {'{S}'};
M_.endo_names_long(11) = {'Price dispersion'};
M_.endo_names(12) = {'Pi_star'};
M_.endo_names_tex(12) = {'{\Pi^*}'};
M_.endo_names_long(12) = {'Optimal reset price'};
M_.endo_names(13) = {'x_aux_1'};
M_.endo_names_tex(13) = {'{x_1}'};
M_.endo_names_long(13) = {'aux. var. 1 recursive price setting'};
M_.endo_names(14) = {'x_aux_2'};
M_.endo_names_tex(14) = {'{x_2}'};
M_.endo_names_long(14) = {'aux. var. 2 recursive price setting'};
M_.endo_names(15) = {'MC'};
M_.endo_names_tex(15) = {'{mc}'};
M_.endo_names_long(15) = {'real marginal costs'};
M_.endo_names(16) = {'M_real'};
M_.endo_names_tex(16) = {'{M/P}'};
M_.endo_names_long(16) = {'real money stock'};
M_.endo_names(17) = {'i_ann'};
M_.endo_names_tex(17) = {'{i^{ann}}'};
M_.endo_names_long(17) = {'annualized nominal interest rate'};
M_.endo_names(18) = {'pi_ann'};
M_.endo_names_tex(18) = {'{\pi^{ann}}'};
M_.endo_names_long(18) = {'annualized inflation rate'};
M_.endo_names(19) = {'r_real_ann'};
M_.endo_names_tex(19) = {'{r^{r,ann}}'};
M_.endo_names_long(19) = {'annualized real interest rate'};
M_.endo_names(20) = {'P'};
M_.endo_names_tex(20) = {'{P}'};
M_.endo_names_long(20) = {'price level'};
M_.endo_names(21) = {'log_m_nominal'};
M_.endo_names_tex(21) = {'{log(M)}'};
M_.endo_names_long(21) = {'log nominal money stock'};
M_.endo_names(22) = {'log_y'};
M_.endo_names_tex(22) = {'{log(Y)}'};
M_.endo_names_long(22) = {'log output'};
M_.endo_names(23) = {'log_W_real'};
M_.endo_names_tex(23) = {'{log(W/P)}'};
M_.endo_names_long(23) = {'log real wage'};
M_.endo_names(24) = {'log_N'};
M_.endo_names_tex(24) = {'{log(N)}'};
M_.endo_names_long(24) = {'log hours'};
M_.endo_names(25) = {'log_P'};
M_.endo_names_tex(25) = {'{log(P)}'};
M_.endo_names_long(25) = {'log price level'};
M_.endo_names(26) = {'log_A'};
M_.endo_names_tex(26) = {'{log(A)}'};
M_.endo_names_long(26) = {'log technology level'};
M_.endo_names(27) = {'log_Z'};
M_.endo_names_tex(27) = {'{log(Z)}'};
M_.endo_names_long(27) = {'log preference shock'};
M_.endo_names(28) = {'nu'};
M_.endo_names_tex(28) = {'{\nu}'};
M_.endo_names_long(28) = {'AR(1) monetary policy shock process'};
M_.endo_partitions = struct();
M_.param_names = cell(13,1);
M_.param_names_tex = cell(13,1);
M_.param_names_long = cell(13,1);
M_.param_names(1) = {'alppha'};
M_.param_names_tex(1) = {'{\alpha}'};
M_.param_names_long(1) = {'capital share'};
M_.param_names(2) = {'betta'};
M_.param_names_tex(2) = {'{\beta}'};
M_.param_names_long(2) = {'discount factor'};
M_.param_names(3) = {'rho_a'};
M_.param_names_tex(3) = {'{\rho_a}'};
M_.param_names_long(3) = {'autocorrelation technology shock'};
M_.param_names(4) = {'rho_nu'};
M_.param_names_tex(4) = {'{\rho_{\nu}}'};
M_.param_names_long(4) = {'autocorrelation monetary policy shock'};
M_.param_names(5) = {'rho_z'};
M_.param_names_tex(5) = {'{\rho_{z}}'};
M_.param_names_long(5) = {'autocorrelation monetary demand shock'};
M_.param_names(6) = {'siggma'};
M_.param_names_tex(6) = {'{\sigma}'};
M_.param_names_long(6) = {'inverse EIS'};
M_.param_names(7) = {'varphi'};
M_.param_names_tex(7) = {'{\varphi}'};
M_.param_names_long(7) = {'inverse Frisch elasticity'};
M_.param_names(8) = {'phi_pi'};
M_.param_names_tex(8) = {'{\phi_{\pi}}'};
M_.param_names_long(8) = {'inflation feedback Taylor Rule'};
M_.param_names(9) = {'phi_y'};
M_.param_names_tex(9) = {'{\phi_{y}}'};
M_.param_names_long(9) = {'output feedback Taylor Rule'};
M_.param_names(10) = {'eta'};
M_.param_names_tex(10) = {'{\eta}'};
M_.param_names_long(10) = {'semi-elasticity of money demand'};
M_.param_names(11) = {'epsilon'};
M_.param_names_tex(11) = {'{\epsilon}'};
M_.param_names_long(11) = {'demand elasticity'};
M_.param_names(12) = {'theta'};
M_.param_names_tex(12) = {'{\theta}'};
M_.param_names_long(12) = {'Calvo parameter'};
M_.param_names(13) = {'tau'};
M_.param_names_tex(13) = {'{\tau}'};
M_.param_names_long(13) = {'labor subsidy'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 3;
M_.endo_nbr = 28;
M_.param_nbr = 13;
M_.orig_endo_nbr = 28;
M_.aux_vars = [];
M_ = setup_solvers(M_);
M_.Sigma_e = zeros(3, 3);
M_.Correlation_matrix = eye(3, 3);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
M_.surprise_shocks = [];
M_.heteroskedastic_shocks.Qvalue_orig = [];
M_.heteroskedastic_shocks.Qscale_orig = [];
options_.linear = false;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
M_.nonzero_hessian_eqs = [1 2 3 4 5 6 8 9 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28];
M_.hessian_eq_zero = isempty(M_.nonzero_hessian_eqs);
M_.orig_eq_nbr = 28;
M_.eq_nbr = 28;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 0 6 34;
 0 7 0;
 0 8 35;
 1 9 0;
 0 10 0;
 0 11 0;
 0 12 0;
 0 13 0;
 0 14 0;
 2 15 36;
 3 16 0;
 0 17 0;
 0 18 37;
 0 19 38;
 0 20 0;
 0 21 0;
 0 22 0;
 0 23 0;
 0 24 0;
 4 25 0;
 0 26 0;
 0 27 0;
 0 28 0;
 0 29 0;
 0 30 0;
 0 31 0;
 0 32 0;
 5 33 0;]';
M_.nstatic = 19;
M_.nfwrd   = 4;
M_.npred   = 4;
M_.nboth   = 1;
M_.nsfwrd   = 5;
M_.nspred   = 5;
M_.ndynamic   = 9;
M_.dynamic_tmp_nbr = [25; 26; 28; 5; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , 'FOC Wages, eq. (2)' ;
  2 , 'name' , 'Euler equation eq. (3)' ;
  3 , 'name' , 'Definition nominal interest rate), p. 22 top' ;
  4 , 'name' , 'Aggregate output, above eq. (14)' ;
  5 , 'name' , 'Definition Real interest rate' ;
  6 , 'name' , 'Monetary Policy Rule, p. 26 bottom/eq. (22)' ;
  7 , 'name' , 'Market Clearing, eq. (15)' ;
  8 , 'name' , 'Technology Shock, eq. (6)' ;
  9 , 'name' , 'Preference Shock, p.54' ;
  10 , 'name' , 'Monetary policy shock' ;
  11 , 'name' , 'Definition marginal cost' ;
  12 , 'name' , 'LOM prices, eq. (7)' ;
  13 , 'name' , 'LOM price dispersion' ;
  14 , 'name' , 'FOC price setting' ;
  15 , 'name' , 'Auxiliary price setting recursion 1' ;
  16 , 'name' , 'Auxiliary price setting recursion 2' ;
  17 , 'name' , 'Definition log output' ;
  18 , 'name' , 'Definition log real wage' ;
  19 , 'name' , 'Definition log hours' ;
  20 , 'name' , 'Annualized inflation' ;
  21 , 'name' , 'Annualized nominal interest rate' ;
  22 , 'name' , 'Annualized real interest rate' ;
  23 , 'name' , 'Real money demand, eq. (4)' ;
  24 , 'name' , 'definition nominal money stock' ;
  25 , 'name' , 'Definition price level' ;
  26 , 'name' , 'Definition log price level' ;
  27 , 'name' , 'Definition log TFP' ;
  28 , 'name' , 'Definition log preference' ;
};
M_.mapping.C.eqidx = [1 2 7 15 16 ];
M_.mapping.W_real.eqidx = [1 11 18 ];
M_.mapping.Pi.eqidx = [2 5 6 12 13 15 16 20 25 ];
M_.mapping.A.eqidx = [4 8 27 ];
M_.mapping.N.eqidx = [1 4 11 19 ];
M_.mapping.R.eqidx = [3 5 6 21 23 ];
M_.mapping.realinterest.eqidx = [5 22 ];
M_.mapping.Y.eqidx = [4 6 7 11 15 16 17 23 ];
M_.mapping.Q.eqidx = [2 3 ];
M_.mapping.Z.eqidx = [2 9 15 16 28 ];
M_.mapping.S.eqidx = [4 11 13 ];
M_.mapping.Pi_star.eqidx = [12 13 14 ];
M_.mapping.x_aux_1.eqidx = [14 15 ];
M_.mapping.x_aux_2.eqidx = [14 16 ];
M_.mapping.MC.eqidx = [11 15 ];
M_.mapping.M_real.eqidx = [23 24 ];
M_.mapping.i_ann.eqidx = [21 ];
M_.mapping.pi_ann.eqidx = [20 ];
M_.mapping.r_real_ann.eqidx = [22 ];
M_.mapping.P.eqidx = [24 25 26 ];
M_.mapping.log_m_nominal.eqidx = [24 ];
M_.mapping.log_y.eqidx = [17 ];
M_.mapping.log_W_real.eqidx = [18 ];
M_.mapping.log_N.eqidx = [19 ];
M_.mapping.log_P.eqidx = [26 ];
M_.mapping.log_A.eqidx = [27 ];
M_.mapping.log_Z.eqidx = [28 ];
M_.mapping.nu.eqidx = [6 10 ];
M_.mapping.eps_a.eqidx = [8 ];
M_.mapping.eps_z.eqidx = [9 ];
M_.mapping.eps_nu.eqidx = [10 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [4 10 11 20 28 ];
M_.exo_names_orig_ord = [1:3];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(28, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(3, 1);
M_.params = NaN(13, 1);
M_.endo_trends = struct('deflator', cell(28, 1), 'log_deflator', cell(28, 1), 'growth_factor', cell(28, 1), 'log_growth_factor', cell(28, 1));
M_.NNZDerivatives = [87; 120; 112; ];
M_.static_tmp_nbr = [13; 3; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
load myparam_values.mat
set_param_value('alppha',myparam(1));
set_param_value('betta',myparam(2));
set_param_value('rho_a',myparam(3));
set_param_value('rho_nu',myparam(4));
set_param_value('rho_z',myparam(5));
set_param_value('siggma',myparam(6));
set_param_value('varphi',myparam(7));
set_param_value('phi_pi',myparam(8));
set_param_value('phi_y',myparam(9));
set_param_value('eta',myparam(10));
set_param_value('epsilon',myparam(11));
set_param_value('theta',myparam(12));
set_param_value('tau',myparam(13));
resid;
steady;
oo_.dr.eigval = check(M_,options_,oo_);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 1e-06;
options_.k_order_solver = true;
options_.irf = 0;
options_.nocorr = true;
options_.nofunctions = true;
options_.nomoments = true;
options_.noprint = true;
options_.order = 3;
options_.periods = 193;
var_list_ = {};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
save("oo_");


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'NK_NL_DSGE_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'NK_NL_DSGE_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'NK_NL_DSGE_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'NK_NL_DSGE_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'NK_NL_DSGE_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'NK_NL_DSGE_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'NK_NL_DSGE_results.mat'], 'oo_recursive_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
