function [C, labels, par_k] = fitko_model_kinetic_distributions(par, t)

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
function [C_k, C_k_all] = C_k_fun(par_k, k, t)
% kinetic matrix

par = par_k(k);

% extract
T = par.T; % time constants
s = par.s; % sigmas
realizations = par.realizations; % realizations
seed = par.seed;
distribution = par.distribution; % distribution

% number of species
N = size(T,1);

%
E = zeros(numel(t),N*N);

if N == 0
    C_k = [];
else
    rng(seed);
    for i = 1:realizations
       
        % gaussian distribution of times
        Tr = parameter_distribution(T,s,distribution);
        
        % generate exponential constant matrix
        M=1./Tr;
        M(Tr==0)=0;

        % matrix of the kinetic model
        R = transpose(M);
        R(linspace(1,numel(R),N)) = - sum(M,2);

        % solve kinetics
        [X,L]=eig(R);
        Ei=zeros(N,N,numel(t));
        for k=1:numel(t)
            S=L;
            S(linspace(1,numel(L),N))=exp(diag(L)*t(k));
            Ei(:,:,k) = X*S/X;
        end

        % rate transfer matrix
        % Butkus et al. Lithuanian Journal of Physics 50:3 (2010) page 285
        Ei = reshape(permute(Ei,[3 2 1]), [numel(t) numel(T)]);

        % ad to average
        E = E + Ei;
    end
    E = E / realizations;
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Gaussian distribution of times
% Gaussian distribution of times
function T = parameter_distribution(T,s,distribution)
    for i = 1:numel(T)
        if s(i) == 0
        else
            switch distribution
                case 'lognormal'
                    T(i) = lognrnd(log(T(i)),s(i)); % lognormal
                case 'gamma'
                    T(i) = gamrnd(s(i),T(i)/s(i)); % gamma
                case 'weibull'
                    T(i) = wblrnd(T(i),s(i)); % weibull
                case 'gaussian'
                    T(i) = T(i) + randn()*s(i); % gaussian
                    if T(i) < 0; T(i) = 1e-4; end % se i tempi sono negativi usa un valore piccolo positivo
            end
        end
    end
end