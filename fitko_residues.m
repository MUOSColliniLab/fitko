function [R, M] = fitko_residues(data, output)
% FITKO_RESIDUES evaluates residues and model function.
%
% USAGE:
% [R, M] = fitko_residues(data, output)
%
% data -- cell of 2DES datasets
% output -- output of fitko
%
% R -- cell of residues with 2DES dataset structure
% M -- cell of model functions with 2DES dataset structure



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Preliminary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cutting
data = fitko_cutting(data,output);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%c% Evaluate model function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% project data on the component matrix
C = output.C;
A = fitko_component_projection(data, C);

% evaluate model
M = fitko_evaluate_model(A, C);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute residues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clone datasets into R
R = data;

% subtract model from data
for i = 1:numel(data)
    R{i}.X = R{i}.X - M{i}.X;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Subfucntions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% preliminary
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

%% evaluate model function

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

function M = fitko_evaluate_model(A, C)

for i = 1:numel(A)
    % extract A and C matrices
    Ai = A{i}.X;
    Ci = C.X;
    
    % matrix multiplication
    X = fitko_etprod('ijk',Ai,'ijn',Ci,'kn');
    
    % saving
    M{i}.X = X;
    M{i}.t = C.t;
    M{i}.f = A{i}.f;
end
end
