function [C, labels, par_k] = fitko_model_generic(par, t)

par_k = @(k) update_par(par,k);
C = @(k) [C_k_fun(par_k, k, t)];

% labels
N = numel(par.e);
alphabet = {'A','B','C','D','E','F','H','I','J','k','M','N','O',...
          'P','Q','R','S','T','U','V','W','X','Y','Z'};
labels = alphabet(1:N);

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

% esponenti
e = par.e;
f = par.functions;

C_k = zeros(numel(t),numel(e)); 
for i = 1:numel(e)
    C_k(:,i) = f{i}(t,e(i));
end

end