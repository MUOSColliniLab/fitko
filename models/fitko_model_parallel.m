function [C, labels, par_k] = fitko_model_parallel(par, t)

par_k = @(k) update_par(par,k);
C = @(k) [C_k_fun(par_k, k, t) C_o_fun(par_k, k, t)];

% labels
N = numel(par.D);
M = numel(par.w);
alphabet = {'A','B','C','D','E','F','H','I','J','k','M','N','O',...
          'P','Q','R','S','T','U','V','W','X','Y','Z'};
labels_k = alphabet(1:N);
numbers_plus = {'1+','2+','3+','4+','5+','6+','7+','8+','9+','10+','11+','12+','13+',...
           '14+','15+','16+','17+','18+','19+','20+','21+'};
numbers_minus = {'1-','2-','3-','4-','5-','6-','7-','8-','9-','10-','11-','12-','13-',...
           '14-','15-','16-','17-','18-','19-','20-','21-'};
labels_o = [numbers_plus(1:M); numbers_minus(1:M)];
labels = [labels_k, labels_o(:)'];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function par = update_par(par, k)
    fix = [];
    lab_ = [];
    for i = 1:numel(par.fitted)
        fix = [fix; vec(par.([par.fitted{i},'_fix']))];
        labi = char( par.fitted{i} * vec(ones(size(par.(par.fitted{i})))) );
        lab_ = [lab_; labi];
    end
    
    clear labi
    lab = lab_(~fix); % label dei parametri k

    for i = 1:numel(par.fitted)
        labi = par.fitted{i}; % label of parameter i
        fixi = fix(lab_==labi); % fix of parameter i
        ki = k(lab==labi); % values of parameter i

        par.(labi)(~fixi) = ki; % update par struct
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kinetics
function C_k = C_k_fun(par_k, k, t)
% kinetic matrix

par = par_k(k);

% extract time constant matrix
D = par.D;

% number of species
N = numel(D);

if N == 0
    C_k = [];
else
    C_k = exp(- t(:)*(1./D));
end

end

% oscillations
function C_o = C_o_fun(par_k, k, t)
% oscillation matrix

par = par_k(k);

w = par.w(:)';
d = par.d(:)';

% apply fitted parameters
N = numel(w);

% time matrix
T = diag(t)*ones(numel(t),2*N);

% extract and convert parameters
freq = 2 * pi * w / 33356; % [cm^{-1]] -> [PHz rad]
rates = - 1 ./d; % [fs] -> [PHz rad]

% exponential factors matrix
b = [rates + 1i*freq; rates - 1i*freq];
b = b(:)';
B = ones(numel(t),2*N) * diag(b);

% oscillating components matrix
C_o = exp(B .* T);

end