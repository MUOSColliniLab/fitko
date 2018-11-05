%% TEST

% generate sample data
D.t=0:5:1000;
D.f=0:100;
gauss=gaussmf(D.f,[10 50]);
A=diag(gauss)*ones(numel(gauss))*diag(gauss);
dinamics=exp((-1/500+1i*2*pi*100/33356)*D.t) ...
    + exp((-1/100+1i*2*pi*400/33356)*D.t) ...
    + 2*exp((-1/250)*D.t);
D.X = fitko_etprod('ijk',A,'ijm',dinamics,'mk');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test all
clear input
input.model = 'parallel';
input.D=[300 50];
input.em_select=1:50;
input.exc_select=1:50;
input.w=[100 200]; 
input.d=[100 300];
input.sub_res = 10;

output = fitko({D,D},input);
[R, M] = fitko_residues({D,D},output);
fitko_check({D,D},output)
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test no dynamics
clear input
input.T=[];
input.em_select=1:50;
input.exc_select=1:50;
input.w=[100 300]; 
input.d=[100 520];
input.sub_res = 10;

output = fitko({D,D},input);
[R, M] = fitko_residues({D,D},output);
fitko_check({D,D},output)
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test no oscillations
clear input
input.T=[0 100 200 0; 0 500 1000 0; 0 0 10000 0; 0 0 0 0];
input.em_select=1:50;
input.exc_select=1:50;
input.w=[]; 
input.d=[];
input.sub_res = 10;

output = fitko({D},input);
[R, M] = fitko_residues({D},output);
fitko_check({D},output)
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test kinetics distributions
clear input
input.model = 'kinetic_distributions';
input.distribution = 'gamma';
input.realizations = 10;
input.T=[0 100 200; 0 500 1000; 0 0 10000];
input.s=[0 0 10; 0 0 0; 0 0 0];
input.em_select=1:50;
input.exc_select=1:50;
input.w=[]; 
input.d=[];
input.sub_res = 10;

output = fitko({D},input);
[R, M] = fitko_residues({D},output);
fitko_check({D},output)
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test oscillations and kinetics
clear input
input.T=[0 100 200; 0 500 1000; 0 0 10000];
input.em_select=1:50;
input.exc_select=1:50;
input.w=[100 200 300]; 
input.d=[100 200 300];
input.sub_res = 10;

output = fitko({D},input);
[R, M] = fitko_residues({D},output);
fitko_check({D},output)
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test parallel
clear input
input.model = 'parallel';
input.D=[100 200 1000];
input.em_select=1:50;
input.exc_select=1:50;
input.w=[]; 
input.d=[];
input.sub_res = 10;

output = fitko({D},input);
[R, M] = fitko_residues({D},output);
fitko_check({D},output)
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test parallel and populations
clear input
input.model = 'parallel';
input.D=[100 200];
% input.em_select=1:50;
% input.exc_select=1:50;
input.w=[]; 
input.d=[];
input.sub_res = 10;

output = fitko({D},input);
[R, M] = fitko_residues({D},output);
fitko_check({D},output)
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test fixed kinetics
clear input
input.T=[100 200 1000; 0 0 100; 0 10 0];
input.T_fix=ones(size(input.T));
% input.em_select=1:50;
% input.exc_select=1:50;
input.w=[]; 
input.d=[];
input.sub_res = 10;

output = fitko({D},input);
[R, M] = fitko_residues({D},output);
fitko_check({D},output)
close all


%% PROTOTIPING

D.t=0:5:1000;
D.f=0:100;
gauss=gaussmf(D.f,[10 50]);
A=diag(gauss)*ones(numel(gauss))*diag(gauss);
dinamics=exp((-1/500+1i*2*pi*100/33356)*D.t) + 2*exp((-1/250)*D.t);
% dinamics=2*exp((-1/250)*D.t);
D.X = fitko_etprod('ijk',A,'ijm',dinamics,'mk');

clear input
input.model = 'kinetic';
input.realizations = 5;
input.T=[0  60 50; 0 0 500; 0 0 0];
input.s=[0  10  0; 0 0   0; 0 0 0];
% input.P=[1 0.3 0.1];
% input.T_fix=ones(size(input.T));
% input.em_select=1:50;
% input.exc_select=1:50;
input.w=[10];
input.d=[10];
% input.sub_res = 10;

output = fitko({D},input);
[R, M] = fitko_residues({D},output);
fitko_check({D},output)
