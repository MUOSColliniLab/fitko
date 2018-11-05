%%% FITKO EXAMPLES %%%
% In this example script the essential functionalities of fitko() will be
% described. The routine is based on the global optimization procedure
% reported in [1].
% [1] Volpato et al., Opt. Express 24 (2016), 24773.

% Generate test data for the examples.
D = generate_test_data();
% Required fields of each dataset D:
%           D.X -- tri-dimensional array with the data (1=nu_3, 2=nu_1, 3=t_2)
%           D.t -- population time axis [fs]
%           D.f -- frequency axis [cm^{-1}]


%% using a parallel model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input parameters of fitko() are passed using a struct variable.
% Type "help fitko" in the command window to get a small description of all
% the inputs. Most of them can be omitted and the routine will proced using
% the internal defaults.

% The simplest model implemented in fitko() is "parallel" which will use a
% sum of exponential decays as fitting model.
% The number of exponential components is controlled via inp.D.
clear inp
inp.model = 'parallel';
inp.D = [200 1000 1e6]; % [fs] decay constants
% The values inserted in inp.D are the initial guesses for the parameters.
% Each exponential component is defined as exp(-t/d), where t is the
% population time and d is the time constant.
% In this case a tri-exponential model will be employed.

% fitting
outp = fitko({D}, inp);
% outp is a struct variable which contains the output of the routine. See
% fitko() help.

% review of the fitting
fitko_check({D}, outp);
% Some figures will be genearted. They contain some relevant information
% that can help you investigate the fitting output.
% In figure 2 are reported the amplitude maps associated to each component.
% In this case we can separate three processes: a relaxation process with a
% time constant of about 500 fs, spectral diffusion of the lower diagonal
% feature with a time constant of about 1 ps, and a slower decay with a
% time constant longer than the investigated spectral window.


%% edit the options available for the fitting routine lsqnonolin()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The convergence of the optimization problem can be followed in the
% command window where the output of the lsqnonlin() routine is displayed.
% Sometimes is neccessary to tweak the default conditions and tollerances
% of the optimization routine in order to achieve a better solution.

% The parameters are passed to the optimization routine using the field
% inp.fit_options. All the defaults can be generated using optimset().
inp.fit_options = optimset('lsqnonlin');
inp.fit_options.TolX = 1e-6;
inp.fit_options.TolFun = 1e-10;
% See help of lsqnonlin() and optimset() for more info on all the possible
% fields controlling the optimization problem.
% Try to write in the command window "inp.fit_options." and then press tab.
% A list of all the parameters should appear.

% fitting
outp = fitko({D}, inp);

% review of the fitting
fitko_check({D}, outp);
% The fitted decay constant for the second exponential component should be
% a little bit different than the previous output.


%% fix a parameter and apply upper and lower bounduaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Some optional parameters of fitko() allow to fix one or more decay
% constant and to impose upper and lower bounduaries.

inp.D_fix = [0 0 1]; % the last time constant is kept fixed to the initial guess
inp.D_lb =  [0    500   0]; % lower bounduaries
inp.D_ub =  [Inf 4000 Inf]; % upper bounduaries
% In this case we fixed the last time constant since it is much longer than
% the investigated time window and we restricted the parameter space for
% the second exponential decay component.

% fitting
outp = fitko({D}, inp);


%% exclude the initial transient from the fit and apply region of interest (ROI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% It is possible to exclude the initial maps from the fitting using
% inp.start_cut.
inp.start_cut = 5;
% In this case the first 5 maps will be discarded.
% It is also possible to select a region of interest using exc_select and
% em_select.
inp.exc_select = 64:100; % excitaton axis selection
inp.em_select = 20:100; % emission axis selection
% The values of inp.exc_selec and inp.em_select are all the indexes of the
% axis D.f which should be considered in the fitting procedure.

% We will now use only two components because the ROI excluded the
% lower diagonal feature where the spectral diffusion is located.
inp.D = [200 1e6];
inp.D_fix = [0 1];
inp.D_lb = [0 0];
inp.D_up = [0 Inf];

% fitting
outp = fitko({D}, inp);

% review of the fitting
fitko_check({D}, outp);
% We obtained a decay constant much closer to the theoretical one of 500
% fs.


%% fitting of exponentially damped oscillations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% It is possible to add damped oscillating components to the fit using 
% inp.w and inp.d.
%
% inp.w = [500 750]; % [cm^{-1}] frequencies
% inp.d = [500 1000]; % [fs] decay constants
%
% As usual you can control fixed parameters and bounduaries using
% inp.w_fix, inp.d_fix, inp.w_lb, inp.d_lb, inp.w_ub and inp.d_ub.


%% multiple datasets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Multiple datasets can be fitted at the same time using the same model.
% For example, this is usefull if you need to fit rephasing and non
% rephasing data simultaneously.

D1 = D; % generate a first copy of the data
D2 = D; % generate a second copy of the data
outp = fitko({D,D1,D2},inp);


%% using a kinetic model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% !!! Warning: this is still a beta !!! Ask Andrea.

% The model "kinetic" can be used to fit the data using a suitable kinetic
% model. The kinetic scheme is defined using matrix inp.T. 
% The values of the matrix represent the time constant associated to the
% kinetic model. The row index represetns the initial state whereas the
% column index represents the arrival state.
% Diagonal time constants define the relaxation to the ground state.

clear inp % clear the variable to reset all the parameters
inp.model = 'kinetic';
inp.T = [0 300; 
         0 1e6];
% With this model we are introducing  a relaxation from stateA to stateB
% with an inital guess of 300 fs and a relaxation to the ground state of
% stateB with a time constant of 1e6 fs.
% Fixing of parameters and upper and lower bounduaries are still possible.
inp.T_fix = [0 0;
             0 1];

% fitting         
outp = fitko({D},inp);

% review of the fitting
fitko_check({D}, outp);

