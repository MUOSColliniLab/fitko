function [C, labels, par_k] = fitko_model_kinetic(par, t)

par_k = @(k) update_par(par,k);
C = @(k) [C_k_fun(par_k, k, t) C_o_fun(par_k, k, t)];

% labels
N = size(par.T,1);
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
T = par.T;

% number of species
N = size(T,1);

if N == 0
    C_k = [];
else
    % generate exponential constant matrix
    M=1./T;
    M(T==0)=0;

    % matrix of the kinetic model
    R = transpose(M);
    R(linspace(1,numel(R),N)) = - sum(M,2);

    % solve kinetics
    [X,L]=eig(R);
    E=zeros(N,N,numel(t));
    for i=1:numel(t)
        S=L;
        S(linspace(1,numel(L),N))=exp(diag(L)*t(i));
        E(:,:,i) = X*S/X;
    end

    % rate transfer matrix
    % Butkus et al. Lithuanian Journal of Physics 50:3 (2010) page 285
    E = reshape(permute(E,[3 2 1]), [numel(t) numel(T)]);
    C_k = E(:, vec(eye(size(T))==1) ); % diagonal kinetics for fitting (basis set)
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