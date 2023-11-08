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
M_.fname = 'NK_DSGE';
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
M_.exo_names(1) = {'ea'};
M_.exo_names_tex(1) = {'ea'};
M_.exo_names_long(1) = {'ea'};
M_.exo_names(2) = {'eu'};
M_.exo_names_tex(2) = {'eu'};
M_.exo_names_long(2) = {'eu'};
M_.exo_names(3) = {'ev'};
M_.exo_names_tex(3) = {'ev'};
M_.exo_names_long(3) = {'ev'};
M_.endo_names = cell(10,1);
M_.endo_names_tex = cell(10,1);
M_.endo_names_long = cell(10,1);
M_.endo_names(1) = {'x'};
M_.endo_names_tex(1) = {'x'};
M_.endo_names_long(1) = {'x'};
M_.endo_names(2) = {'pi'};
M_.endo_names_tex(2) = {'pi'};
M_.endo_names_long(2) = {'pi'};
M_.endo_names(3) = {'int'};
M_.endo_names_tex(3) = {'int'};
M_.endo_names_long(3) = {'int'};
M_.endo_names(4) = {'a'};
M_.endo_names_tex(4) = {'a'};
M_.endo_names_long(4) = {'a'};
M_.endo_names(5) = {'u'};
M_.endo_names_tex(5) = {'u'};
M_.endo_names_long(5) = {'u'};
M_.endo_names(6) = {'v'};
M_.endo_names_tex(6) = {'v'};
M_.endo_names_long(6) = {'v'};
M_.endo_names(7) = {'r'};
M_.endo_names_tex(7) = {'r'};
M_.endo_names_long(7) = {'r'};
M_.endo_names(8) = {'y'};
M_.endo_names_tex(8) = {'y'};
M_.endo_names_long(8) = {'y'};
M_.endo_names(9) = {'yn'};
M_.endo_names_tex(9) = {'yn'};
M_.endo_names_long(9) = {'yn'};
M_.endo_names(10) = {'n'};
M_.endo_names_tex(10) = {'n'};
M_.endo_names_long(10) = {'n'};
M_.endo_partitions = struct();
M_.param_names = cell(11,1);
M_.param_names_tex = cell(11,1);
M_.param_names_long = cell(11,1);
M_.param_names(1) = {'beta'};
M_.param_names_tex(1) = {'beta'};
M_.param_names_long(1) = {'beta'};
M_.param_names(2) = {'om'};
M_.param_names_tex(2) = {'om'};
M_.param_names_long(2) = {'om'};
M_.param_names(3) = {'eta'};
M_.param_names_tex(3) = {'eta'};
M_.param_names_long(3) = {'eta'};
M_.param_names(4) = {'si'};
M_.param_names_tex(4) = {'si'};
M_.param_names_long(4) = {'si'};
M_.param_names(5) = {'phi_pi'};
M_.param_names_tex(5) = {'phi\_pi'};
M_.param_names_long(5) = {'phi_pi'};
M_.param_names(6) = {'phi_x'};
M_.param_names_tex(6) = {'phi\_x'};
M_.param_names_long(6) = {'phi_x'};
M_.param_names(7) = {'rho_a'};
M_.param_names_tex(7) = {'rho\_a'};
M_.param_names_long(7) = {'rho_a'};
M_.param_names(8) = {'rho_v'};
M_.param_names_tex(8) = {'rho\_v'};
M_.param_names_long(8) = {'rho_v'};
M_.param_names(9) = {'rho_u'};
M_.param_names_tex(9) = {'rho\_u'};
M_.param_names_long(9) = {'rho_u'};
M_.param_names(10) = {'lam'};
M_.param_names_tex(10) = {'lam'};
M_.param_names_long(10) = {'lam'};
M_.param_names(11) = {'kap'};
M_.param_names_tex(11) = {'kap'};
M_.param_names_long(11) = {'kap'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 3;
M_.endo_nbr = 10;
M_.param_nbr = 11;
M_.orig_endo_nbr = 10;
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
options_.linear = true;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
M_.nonzero_hessian_eqs = [];
M_.hessian_eq_zero = isempty(M_.nonzero_hessian_eqs);
M_.orig_eq_nbr = 10;
M_.eq_nbr = 10;
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
 0 4 14;
 0 5 15;
 0 6 0;
 1 7 16;
 2 8 0;
 3 9 0;
 0 10 0;
 0 11 0;
 0 12 0;
 0 13 0;]';
M_.nstatic = 5;
M_.nfwrd   = 2;
M_.npred   = 2;
M_.nboth   = 1;
M_.nsfwrd   = 3;
M_.nspred   = 3;
M_.ndynamic   = 5;
M_.dynamic_tmp_nbr = [0; 0; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , 'x' ;
  2 , 'name' , 'pi' ;
  3 , 'name' , 'int' ;
  4 , 'name' , 'a' ;
  5 , 'name' , 'u' ;
  6 , 'name' , 'v' ;
  7 , 'name' , 'r' ;
  8 , 'name' , 'yn' ;
  9 , 'name' , 'y' ;
  10 , 'name' , 'n' ;
};
M_.mapping.x.eqidx = [1 2 3 9 ];
M_.mapping.pi.eqidx = [1 2 3 7 ];
M_.mapping.int.eqidx = [1 3 7 ];
M_.mapping.a.eqidx = [1 4 8 10 ];
M_.mapping.u.eqidx = [2 5 ];
M_.mapping.v.eqidx = [3 6 ];
M_.mapping.r.eqidx = [7 ];
M_.mapping.y.eqidx = [9 10 ];
M_.mapping.yn.eqidx = [8 9 ];
M_.mapping.n.eqidx = [10 ];
M_.mapping.ea.eqidx = [4 ];
M_.mapping.eu.eqidx = [5 ];
M_.mapping.ev.eqidx = [6 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [4 5 6 ];
M_.exo_names_orig_ord = [1:3];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(10, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(3, 1);
M_.params = NaN(11, 1);
M_.endo_trends = struct('deflator', cell(10, 1), 'log_deflator', cell(10, 1), 'growth_factor', cell(10, 1), 'log_growth_factor', cell(10, 1));
M_.NNZDerivatives = [34; 0; -1; ];
M_.static_tmp_nbr = [0; 0; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
load myparam_values.mat
set_param_value('beta',myparam(1));
set_param_value('om',myparam(2));
set_param_value('eta',myparam(3));
set_param_value('si',myparam(4));
set_param_value('phi_pi',myparam(5));
set_param_value('phi_x',myparam(6));
set_param_value('rho_a',myparam(7));
set_param_value('rho_v',myparam(8));
set_param_value('rho_u',myparam(9));
M_.params(10) = (1+M_.params(3))/(M_.params(3)+M_.params(4));
lam = M_.params(10);
M_.params(11) = (M_.params(3)+M_.params(4))*(1-M_.params(2))*(1-M_.params(2)*M_.params(1))/M_.params(2);
kap = M_.params(11);
steady;
oo_.dr.eigval = check(M_,options_,oo_);
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (1)^2;
M_.Sigma_e(2, 2) = (1)^2;
M_.Sigma_e(3, 3) = (1)^2;
options_.irf = 0;
options_.nocorr = true;
options_.nofunctions = true;
options_.nomoments = true;
options_.noprint = true;
options_.order = 2;
options_.periods = 193;
var_list_ = {};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
save("oo_");


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'NK_DSGE_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'NK_DSGE_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'NK_DSGE_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'NK_DSGE_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'NK_DSGE_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'NK_DSGE_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'NK_DSGE_results.mat'], 'oo_recursive_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
