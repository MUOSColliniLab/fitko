function D = generate_test_data()
% axes
f = linspace(15000, 19000, 128);
t = 0:10:2000;
X = zeros(numel(f), numel(f), numel(t));

% parameters
e = [16250 17500]; % energies
so = [150 150]; % omogeneous broadening sigma
si = [250 250]; % inhomogeneous broadening sigma
ss = [50 50]; % stokes shift
t_sd = [1000 1000]; % spectral diffusion time
t_tr = [1e6 0
        500  0]; % tranfer times i=initial_state, j=final_state
    
% create data
X = X + feature(f,e(1),e(1),      so(1),si(1),        exp(-t./t_sd(1))     .* exp(-t/t_tr(1,1))     );
X = X + feature(f,e(1),e(1)-ss(1),si(1),si(1),        (1-exp(-t./t_sd(1))) .* exp(-t/t_tr(1,1))     ); 
X = X + 0.5*feature(f,e(2),e(2),      so(2),si(2),    exp(-t./t_sd(2))     .* exp(-t/t_tr(2,1))     );
X = X + 0.5*feature(f,e(2),e(2)-ss(2),si(2),si(2),    (1-exp(-t./t_sd(2))) .* exp(-t/t_tr(2,1))     ); 
X = X + 0.5*feature(f,e(2),e(1)-ss(1),si(1),si(1),    exp(-t./t_tr(1,1))   .* (1-exp(-t/t_tr(2,1))) );

% output
D.X = X + 0.03*randn(size(X));
D.t = t;
D.f = f;
D.type = 'T';
end

function X = feature(f,ex,ey,so,si,time_evolution)
    [X,Y] = meshgrid(f,f);
    gauss =@(f,sig) exp(-(f).^2./(2*sig^2));
    M = gauss((X+Y+(ex-ey))/2-ex, si) .* gauss((X-Y-(ex-ey))/2, so);
    X = fitko_etprod('ijk', M, 'ijm', time_evolution,'mk');
end