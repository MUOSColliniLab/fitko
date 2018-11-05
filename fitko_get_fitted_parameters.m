function [k0, lb, ub, klabels, par_k] = fitko_get_fitted_parameters(par)

[k_, lb_, ub_, fix, ~, klabels_] = vec_par(par);

k0 = k_(~fix); 
lb = lb_(~fix);
ub = ub_(~fix);
klabels = klabels_(~fix);

par_k = @(k) update_par(par, k);

end

function [k, lb, ub, fix, lab, klabels] = vec_par(par)
    k = [];
    lb = [];
    ub = [];
    fix = [];
    lab = [];
    klabels = char(zeros(0,3));
    for i = 1:numel(par.fitted)
        k =   [k;   vec(par.(par.fitted{i}         ))];
        lb =  [lb;  vec(par.([par.fitted{i},'_lb'] ))];
        ub =  [ub;  vec(par.([par.fitted{i},'_ub'] ))];
        fix = [fix; vec(par.([par.fitted{i},'_fix']))];
        
        labi = char( par.fitted{i} * vec(ones(size(par.(par.fitted{i})))) );
        lab = [lab; labi];
       
        klabelsi = [labi num2str((1:numel(labi))','%02i')];
        klabels = [klabels;  klabelsi];
    end
    klabels = cellstr(klabels);
end

function par = update_par(par, k)
    [~ , ~, ~, fix, lab_] = vec_par(par); % tutti i parametri
    lab = lab_(~fix); % label dei parametri k

    for i = 1:numel(par.fitted)
        labi = par.fitted{i}; % label of parameter i
        fixi = fix(lab_==labi); % fix of parameter i
        ki = k(lab==labi); % values of parameter i

        par.(labi)(~fixi) = ki; % update par struct
    end
end



