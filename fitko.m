function output = fitko(data,input)
% FITKO - Fit of kinetics and oscillations of 2DES data
%         The core computation is based on the variable projection
%         algorithm and an optimization via lsqnonlin().
%         Singular value decomposition svd() is used to compressed the data 
%         and to optimize computation time.
%
% USAGE:
% output = fitko(data, input)
%
% data -- cell of 2DES datasets
%       required fields of each dataset D:
%           D.X -- data (1=nu_3, 2=nu_1, 3=t_2);
%           D.t -- population time axis [fs]
%           D.f -- frequency axis
%
% input -- input parameters struct
%           input.model -- 'parallel'
%                          'kinetic'
%                          'kinetic_distributions'
%       Main fields:
%           'parallel'
%               input.D -- time constants of kinetics [fs]
%               input.w -- frequencies of oscillations [cm^{-1}]
%               input.d -- damping times of oscillations [fs]
%
%           'kinetic'
%               input.T -- time constants of kinetics [fs]
%               input.w -- frequencies of oscillations [cm^{-1}]
%               input.d -- damping times of oscillations [fs]
%
%           'kinetic_distributions'
%               input.T -- time constants of kinetics [fs]
%               input.s -- sigmas of time constants distributions [fs]
%               input.realizations -- number of realizations
%               input.distribution -- name of distribution ['gaussian','gamma','lognormal','weibull']
%               input.seed -- seed of the random number generations
%               input.w -- frequencies of oscillations [cm^{-1}]
%               input.d -- damping times of oscillations [fs]
%
%
%       Other relevant input parameters are:
%
%           input.fit_options -- fit options of lsqnonlin(). Use optimset()
%                                for see all the default fields
%                                input.fit_options = optimset('lsqnonlin');
% 
%           input.start_cut -- number of t2 points to remove from the start
%                              (0 is the default)
%           input.end_cut   -- number of t2 points to remove from the end
%                              (0 is the default)
%
%           input.sub_res   -- number of points of the frequency axes nu1 and nu3 in
%                              the subsampled data for core computations
%                              (32 is the default)
%
%           input.exc_select -- index interval to keep in core computation
%                               for the excitation axis
%                               ([] is the default, meaning full axis)
%                               example input.exc_selection = 25:50;
%           input.em_select  -- index interval to keep in core computation
%                               for the emission axis
%                               ([] is the default, meaning full axis)
%                               example input.exc_selection = 25:50;
%           input.svd_cutoff -- cutoff of the singular value of the svd
%                               decomposition of the data (0.0025);
%
% output -- output parameters struct
%       Same fields as input, plus output of the fit:
%           output.T, .D, .w, .d -- parameters updated with fit results
%           output.A -- amplitude maps (cells of datatsets)
%           output.C -- component matrix struct
%           output.regression -- regression statistics
%           output.version -- version of the fitko routine
%
% Andrea Volpato
% andrea.volpato@unipd.it
version='alpha_v0.8_r20180711';

% changelog
% v0.8
% - small bug fixes
%
% v0.7
% - major reshaping of the code, now is much more general
% - implementation of model selection
%
% v0.6
% - major polishing of the code
% - removed population fitting
% - added kinetics projections matrix
%
% v0.5
% - fitting with diagonal basis set
% - fitting of starting populations
%
% v0.4
% - change from population dynamics to population transfer rates
% - refined parallel model selection, now all T T_fix T_lb and T_ub are vectors in input
% - added input.exc_select and input.em_select

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters manipulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of kinetic species
if isfield(input,'T')==1; N_k = size(input.T,2);
elseif isfield(input,'D')==1; N_k = numel(input.D);
elseif isfield(input,'e')==1; N_k = numel(input.e);
else; N_k = 0; end

% number of oscillating components
if isfield(input,'w')==1; N_o = numel(input.w); 
else; N_o = 0; end

% make default parameters struct
par = fitko_make_default(N_k,N_o);

% merge input and default parameters struct
fieldNames = fieldnames(input);
for i=1:numel(fieldNames)
    par.(fieldNames{i}) = input.(fieldNames{i});
end

% starting manipulation of parameters
par = fitko_digest_par(par);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Subsampling and cutting of 2DES data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cut according to start_cut and end_cut
data = fitko_cutting(data,par);

% select portion of maps according to exc_select em_select
data = fitko_exc_em_selection(data,par);

% subsample and reshape
[X, t] = fitko_subsampling(data,par);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Model definions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch par.model
    case 'parallel'
        par.fitted = {'D','w','d'};
        [k0, lb, ub, klabels] = fitko_get_fitted_parameters(par);
        [C, labels, par_k] = fitko_model_parallel(par, t);
        
    case 'kinetic'
        par.fitted = {'T','w','d'};
        [k0, lb, ub, klabels] = fitko_get_fitted_parameters(par);
        [C, labels, par_k] = fitko_model_kinetic(par, t);
        
    case 'kinetic_distributions'
        par.fitted = {'T','s','w','d'};
        [k0, lb, ub, klabels] = fitko_get_fitted_parameters(par);
        [C, labels, par_k] = fitko_model_kinetic_distributions(par, t);
    
    %%% beta models !!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'polynomial'
        par.fitted = {'e'};
        [k0, lb, ub, klabels] = fitko_get_fitted_parameters(par);
        [C, labels, par_k] = fitko_model_polynomial(par, t);
        
    case 'generic'
        par.fitted = {'e'};
        [k0, lb, ub, klabels] = fitko_get_fitted_parameters(par);
        [C, labels, par_k] = fitko_model_generic(par, t);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Variable projection algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% minimization and statistics
% if no parameter is fitted skip the minimization
if numel(k0) == 0
    regression = [];
    k = [];
else
    [k, regression] = fitko_varpro(X,C,k0,lb,ub,...
                                   par.fit_options, par.svd_cutoff);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% update fitted parameters
output = par_k(k);
output.regression = regression; % output regression statistics
output.regression.parlabels = klabels;
output.regression.parameters = k;

% component matrix with only the kinetic basis set
C_struct.X = C(k);
C_struct.t = t(:);
C_struct.labels = labels;
output.C = C_struct;

% amplitude maps of full resolution data
A_basis = fitko_component_projection(data, C_struct);
output.A = A_basis;

% add version lable of fitko routine
output.fitko_version = version;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Subfucntions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% -> 1. parameters manipulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function par = fitko_make_default(N_k,N_o)
% generate all the defaults parameters

% parameters
par.T = zeros(N_k);
par.s = zeros(N_k);
par.D = zeros(1,N_k);
par.w = zeros(1,N_o);
par.d = zeros(1,N_o);
par.e = zeros(1,N_k);

% fix switch
par.T_fix = zeros(N_k);
par.s_fix = zeros(N_k);
par.D_fix = zeros(1,N_k);
par.w_fix = zeros(1,N_o);
par.d_fix = zeros(1,N_o);
par.e_fix = zeros(1,N_k);

% lower bounduaries
par.T_lb = 0*ones(N_k);
par.s_lb = 0*ones(N_k);
par.D_lb = 0*ones(1,N_k);
par.w_lb = 0*ones(1,N_o);
par.d_lb = 0*ones(1,N_o);
par.e_lb = 0*ones(1,N_k);

% upper bounduaries
par.T_ub = Inf*ones(N_k);
par.s_ub = Inf*ones(N_k);
par.D_ub = Inf*ones(1,N_k);
par.w_ub = Inf*ones(1,N_o);
par.d_ub = Inf*ones(1,N_o);
par.e_ub = Inf*ones(1,N_k);

% fitting fit_options
par.fit_options = optimset('lsqnonlin'); 

% additional fit_options
par.start_cut = 0;
par.end_cut = 0;
par.sub_res = 32;
par.em_select = [];
par.exc_select = [];
par.svd_cutoff = 0.0025;
par.model = 'kinetic';
par.seed = 1;

end

function par = fitko_digest_par(par)
% make w and d row vectors
par.w = par.w(:)';
par.d = par.d(:)';

% fix T parameters when they are zero
par.T_fix(par.T == 0) = 1;
par.s_fix(par.s == 0) = 1;

% enable display iteration of the minimization algorithm
par.fit_options.Display = 'iter';
end



%% -> 3 subsampling and cutting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = fitko_cutting(data,par)

% extract time axis
t = data{1}.t;

% cut the time axis and the subsampled using par
start_cut = par.start_cut;
end_cut = par.end_cut;
t = t(1+start_cut:end-end_cut);
for i = 1:numel(data)
    data{i}.X = data{i}.X(:,:,1+start_cut:end-end_cut);
    data{i}.t = t;
end

end

function data = fitko_exc_em_selection(data,par)

% extract data
f = data{1}.f;
exc_select = par.exc_select;
em_select = par.em_select;

% default
if numel(exc_select)==0
    exc_select=1:numel(f);
end
if numel(em_select)==0
    em_select=1:numel(f);
end

% boulean selection
exc_selbou = zeros(numel(f),1);
exc_selbou(exc_select) = 1;
exc_selbou = logical(exc_selbou);
em_selbou = zeros(numel(f),1);
em_selbou(em_select) = 1;
em_selbou = logical(em_selbou);

% select area of interest, zeroing the rest
for i = 1:numel(data)
    data{i}.X(~em_selbou,:,:)=0;
    data{i}.X(:,~exc_selbou,:)=0;
end

end

function [X, t] = fitko_subsampling(data,par)

% extract time axis
t = data{1}.t;

% extract subsampling resolution
sub_res = par.sub_res;

% data extraction and subsampling
X = [];
for i = 1:numel(data)
    D = data{i};
    Y = imresize3(D.X, [sub_res sub_res numel(t)],'nearest');
    Z = reshape(Y, [sub_res*sub_res numel(t)]);
    Z = transpose(Z);
    X = [X Z];
end

% remove null columns
x = sum(abs(X));
X = X(:,x~=0);

end



%% -> 4. variable projection algorithm minimization %%%%%%%%%%%%%%%%%%%%%%%

function  [k, regression] = fitko_varpro(X,C,k0,lb,ub,fit_options,svd_cutoff);

% svd compression of data
[U,S,V] = svd_compression(X,svd_cutoff); % cutoff 0.25%

% residues function
f_lsq = @(k) fitko_varpro_residues(U,S,V,C,k);

% lsqnonlin minimization and statistics
[k,wresid_norm2,~,exitflag,output,~,J] = lsqnonlin(f_lsq,k0,lb,ub,fit_options);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Statistical calculation (adapted from varpro.m)

% evaluate components matrix
Ce=C(k);

% Compute the degree of freedom df.
[~,n_C]=size(Ce);
m = numel(X); % number of observation
q = numel(k0); % number of nonlinear parameters
n = size(X,2)*n_C; % number of linear parameters
df = m-n-q;
df_reduced = m-q;

% Insert default output of lsqnonlin in Regression.report.
regression.report = output;
wresid_norm = sqrt(wresid_norm2);
regression.report.exitflag = exitflag;

% Calculate sample variance,  the norm-squared of the residual
%    divided by the number of degrees of freedom.
sigma2 = wresid_norm2 / df;

% Compute  Regression.sigma:        
%                square-root of weighted residual norm squared divided 
%                by number of degrees of freedom.
regression.sigma = sqrt(sigma2);

% Compute  Regression.R_squared:
%                The coeficient of determination for the fit,
%                also called the square of the multiple correlation
%                coefficient, or R^2.
%                It is computed as 1 - wresid_norm^2/CTSS,
%                where the "corrected total sum of squares"
%                CTSS is the norm-squared of W*(y-y_bar),
%                and the entries of y_bar are all equal to
%                (the sum of W_i y_i) divided by (the sum of W_i).
X_bar = mean(X(:));
CTTS = norm(X(:)-X_bar)^2;
regression.R_squared = 1 - wresid_norm^2/CTTS;

% Compute  Regression.RMS = sigma^2:
%                the weighted residual norm divided by the number of
%                degrees of freedom.
%                RMS = wresid_norm / sqrt(m-n+q)
regression.RMS = sigma2;

% Compute  Regression.Cov
%               the covariance matrix. Then, with corr2cov() function compute
%               standard deviation (Regression.Std) and correlation matrix
%               (Regression.Corr). The confidence interval (Regression.Conf)
%               is computed from standard devation.
% see ref. Mullen VanStokkum, Numer. Algor. 2009 51 pag 326
Cov = sigma2.*inv(transpose(J)*J);
[Std,Corr] = cov2corr(full(Cov));
Std = Std(:);
Conf = 1.96*Std;
Conf_relative = abs(Conf./k);

% save into regression field
regression.Cov = Cov;
regression.Std = Std;
regression.Corr = Corr;
regression.Conf = Conf;
regression.Conf_relative = Conf_relative;
regression.DF=df;
regression.DF_reduced=df_reduced;

% Evaluate Rank of components matrix and warn if singular
myrank=rank(Ce);
mycomp=size(Ce,2);
if myrank<mycomp
    fprintf('*** Warning: The linear parameters are not well-determined.\n')
    fprintf('*** The number of conponents is %d.\n',mycomp)
    fprintf('*** The rank of the component matrix is %d.\n',myrank)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

function [U,S,V] = svd_compression(X,cutoff)
% svd compression of data data retaining only singular values greater than
% cutoff

% svd
[U,S,V] = svd(X,'econ');

% singolar values to keep
keep = diag(S)>cutoff*max(S(:));

% discard small singular values and vectors
U = U(:,keep);
S = S(keep,keep);
V = V(:,keep);
end

function f = fitko_varpro_residues(U,S,V,C,k)
% reference: KM Mullen IHM van Stokkum Numer Alogr (2009) 51:319-340

% evaluate component matrix
Ce=C(k);

% QR decomposition of Ce
[Q,~]=qr(Ce);
Q2 = Q(:,size(Ce,2)+1:end);

% evaluate residues
% f=Q2'*X; % attivare questa per accellerare il calcolo
f=Q2*transpose(Q2)*U*S*V';

f=[real(f(:)); imag(f(:))];

end



%% -> 5. output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = fitko_component_projection(data, C_struct)

% component matrix
Ce = C_struct.X;

% pseudo inverse
Cp = pinv(Ce);

% project the matrix of component on data
for i = 1:numel(data)
    Xi = data{i}.X;
    AiX = fitko_etprod('ijk',Cp,'km',Xi,'ijm');
    Ai.X = AiX;
    Ai.f = data{i}.f;
    Ai.labels = C_struct.labels;
    A{i} = Ai;
end

end
