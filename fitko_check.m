function fitko_check(data, output)
% FITKO_CHECK generates a set of figure to check the output of the fit.
%
% USAGE:
% fitko_check(data, output)
%
% data -- cell of 2DES datasets
% output -- output of fitko()

f1=figure(1);
f2=figure(2);
f3=figure(3);
f4=figure(4);
set(f1, 'WindowStyle', 'docked');
set(f2, 'WindowStyle', 'docked');
set(f3, 'WindowStyle', 'docked');
set(f4, 'WindowStyle', 'docked');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Evaluate model and residues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[R, M] = fitko_residues(data,output);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Correlation Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

correlation_matrix(f1,output);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Component Report
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

component_report(f2,output);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Data and Fit Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

explore_fit(f3,data,M,R);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Kinetics projections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch output.model
    case 'kinetic'
        kinetic_projections(f4,output.T,output.C);
end

%%%
f1=figure(1);

end

%% 1.

function correlation_matrix(h_figure,output)

% extract regression statistics
regression = output.regression;
labels = regression.parlabels;

% select figure
figure(h_figure); clf;
set(h_figure,'name','Correlation Matrix')

% imagesc of regression
if isfield(regression,'Corr')==1
subplot 122
imagesc(abs(regression.Corr))
axis square
title('Correlation matrix')
colorbar
caxis([0 1]);
colormap(gray);

set(gca,'xtick',1:numel(labels))
set(gca,'ytick',1:numel(labels))
set(gca,'xticklabel',labels)
set(gca,'yticklabel',labels)
end

% table with parameters
if numel(output.regression.parameters)==0; output.regression.Conf_relative = []; end
t = uitable(h_figure);
d = [ output.regression.parlabels(:),...
      num2cell(output.regression.parameters(:)),...
      num2cell(output.regression.Conf_relative(:)*100)  ];
t.Data = d;
t.Units = 'normalized';
t.Position = [0.1 0.1 0.3 0.8];
t.ColumnName = {'Label','Value','StdError %'};

uicontrol('Style','text','unit','normalized','position',[0 0 1 0.075],...
    'fontsize',10,...
    'string',['Legend:   T - time constant of kinetic model [fs]    ',...
                        'D - time constant of parallel model[fs]     ',...
                        'w - frequency of a damped oscillation [cm^{-1}]    ',...
                        'd - decay time of a damped oscillation [fs]']);
end

%% removed.

function plot_populations(h_figure,output)

% extract regression statistics
X = output.Pt.X;
t = output.Pt.t;
labels = output.Pt.labels;

% select figure
figure(h_figure); clf;
set(h_figure,'name','Populations')

% imagesc of regression
plot(t,X)
legend(labels)

end


%% 2.

function component_report(h_figure,output)

% select figure 3
figure(h_figure); clf;
set(h_figure,'name','Associated Spectra')

% number of components
n=numel(output.C.labels);
if n==1; n=2; end % slider compatibility with one dataset

% number of datasets
N=numel(output.A);
if N==1; N=2; end % slider compatibility with one dataset

% dataset index slider
uicontrol('style','text',...
          'string','Select Dataset',...
          'unit','normalized','position',[0.7,0.06,0.2,0.03])
hindex = uicontrol('style','slider','unit','normalized','position',[0.7,0.03,0.2,0.03],...
                'value',1, 'min',1, 'max',N, 'sliderstep', [1 1]/(N-1));            

% components slider
uicontrol('style','text',...
          'string','Select Component and Refresh Figure',...
          'unit','normalized','position',[0.4,0.06,0.2,0.03]);
hcomp = uicontrol('Parent',h_figure,'Style','slider',...
                'unit','normalized','Position',[0.4,0.03,0.2,0.03],...
                'value',1, 'min',1, 'max',n, 'sliderstep', [1 1]/(n-1),...
                'callback',{@callback_component_plot,output,hindex});
            
% table with parameters
t = uitable(h_figure);
d = [ num2cell(output.w(:)),...
      num2cell(output.d(:))  ];
t.Data = d;
t.Units = 'normalized';
t.Position = [0.6 0.1 0.3 0.3];
t.ColumnName = {'w / cm^{-1}','d / fs'};

end             

function callback_component_plot(h,evt,output,hindex)

% get indexes for extract data to plot
component_index = round(h.Value);
dataset_index = round(hindex.Value);
h.Value = component_index;
hindex.Value = dataset_index;

% extract relevant data according to dataset_index
A = output.A{dataset_index};
C = output.C;
T = output.T;

% number of stuffs
N_0 = sum(T(:)~=0 | vec(eye(size(T))==1));
N_f = length(output.w);

% extract more
Ai = A.X(:,:,component_index);
Ai_phase = angle(Ai);
f = A.f;
t = C.t;

% maximum amplitude
M=max(abs(Ai(:)));

% plot real map
subplot 231
cla
imagesc(f,f,real(Ai))
set(gca,'ydir','normal')
line(f,f,'color','k')
axis square
colorbar
caxis([-M M])
xlabel('Excitation Frequency / cm^{-1}')
ylabel('Emission Frequency / cm^{-1}')
title('Real Amplitude Map')

% plot abs map
subplot 232
cla
imagesc(f,f,abs(Ai))
set(gca,'ydir','normal')
line(f,f,'color','k')
axis square
colorbar
caxis([-M M])
xlabel('Excitation Frequency / cm^{-1}')
ylabel('Emission Frequency / cm^{-1}')
title('Absolute Amplitude Map')

% plot phase map
subplot 233
cla
imagesc(f,f,Ai_phase)
set(gca,'ydir','normal')
line(f,f,'color','k')
axis square
colorbar
colormap(hsv)
caxis([-pi pi])
xlabel('Excitation Frequency / cm^{-1}')
ylabel('Emission Frequency / cm^{-1}')
title('Phase Map')

% plot component
subplot 223
cla
plot(t, real(C.X(:,component_index)), 'k')
hold on
plot(t, imag(C.X(:,component_index)), 'k--')
hold off
xlabel('Time / fs')
title([output.C.labels{component_index},' Component'])
% ylim(1.05*[-1 1])
legend({'Real','Imag'})
legend('boxoff')

colormap(parula)

end

%% 3.

function explore_fit(h_figure,data,M,R)

% select figure
figure(h_figure); clf;
set(h_figure,'name','Explore fit')

% number of datasets
N=numel(data);
if N==1; N=2; end % slider compatibility with one dataset

uicontrol('style','text','string','Select Dataset -- click on slider','unit','normalized','position',[0.1,0.73,0.2,0.03])

if N==1; N=N+0.01; end % per non incorrere errori con lo slider
hs4 = uicontrol('Parent',h_figure,'Style','slider',...
                'value',1, 'min',1, 'max',N, 'sliderstep', [1 1]/(N-1),...
                'unit','normalized','Position',[0.1,0.7,0.2,0.03],...
                'callback',{@update_mean_map,data});

hbtn = uicontrol(h_figure,...
    'Style','pushbutton',...
    'String','Select a point',...
    'unit','normalized','Position',[0.33 0.7 0.05 0.05],...
    'callback',{@callback_explore_fit,data,M,R,hs4}); 

end

function update_mean_map(h,evt,data)

% get the slider value
dataset_index = round(h.Value);
h.Value = dataset_index;

% extract data
X_data = data{dataset_index}.X;
f = data{dataset_index}.f;

% Mean map
Dm=mean(X_data,3);
M=max(abs(Dm(:)));
ax=subplot(2,3,2); cla;
imagesc(f,f,real(Dm));
set(gca,'ydir','normal')
line(f,f,'color','k')
axis square
colormap parula
caxis(M*[-1 1])
xlabel('Excitation Frequency / cm^{-1}')
ylabel('Emission Frequency / cm^{-1}')
title('Mean Real Map')

end

function callback_explore_fit(h,evt,data,M,R,hs4)

% get the slider value
dataset_index = round(hs4.Value);
hs4.Value = dataset_index;

% extract data
X_fit = M{dataset_index}.X;
X_data = data{dataset_index}.X;
X_res = R{dataset_index}.X;
t_fit = M{dataset_index}.t;
t_data = data{dataset_index}.t;
t_res = R{dataset_index}.t;
f = data{dataset_index}.f;

% Mean map
Dm=mean(X_data,3);
M=max(abs(Dm(:)));
ax=subplot(2,3,2); cla;
imagesc(f,f,real(Dm));
set(gca,'ydir','normal')
line(f,f,'color','k')
axis square
colormap parula
caxis(M*[-1 1])
xlabel('Excitation Frequency / cm^{-1}')
ylabel('Emission Frequency / cm^{-1}')
title('Mean Real Map')

[x,y]=getpts(ax);
x=x(1); y=y(1);
[~,xind]= min(abs(f-x));    
[~,yind]= min(abs(f-y));

signal = permute(X_data(yind,xind,:),[3 2 1]);
fit = permute(X_fit(yind,xind,:),[3 2 1]);
res = permute(X_res(yind,xind,:),[3 2 1]);

% update title of mean map figure
subplot 232
title(['Mean Map (indexes:',num2str(yind),', ',num2str(xind),')'])
hold on
scatter(x,y,'k')
hold off

% real signal and fit
subplot 223
plot(t_data,real(signal),'k');
hold on
plot(t_fit,real(fit),'r');
plot(t_res,real(res),'b')
hold off
legend({'Signal','Fit','Residues'})
legend('boxoff')
title(['Real'])
xlabel('Time / fs')

% imag signal and fit
subplot 224
plot(t_data,imag(signal),'k');
hold on
plot(t_fit,imag(fit),'r');
plot(t_res,imag(res),'b')
hold off
legend({'Signal','Fit','Residues'})
legend('boxoff')
title(['Imaginary'])
xlabel('Time / fs')


end

%% 4.

function kinetic_projections(h_figure, T, C)

t=C.t;

% SOLVE KINETICS
% generate exponential constant matrix
M=1./T;
M(T==0)=0;

% matrix of the kinetic model
R = transpose(M);
R(linspace(1,numel(R),size(T,1))) = - sum(M,2);

% solve kinetics
[X,L]=eig(R);
E=zeros(size(T,1),size(T,2),numel(t));
for i=1:numel(t)
    S=L;
    S(linspace(1,numel(L),size(T,1)))=exp(diag(L)*t(i));
    E(:,:,i) = X*S/X;
end

% diagonal components
C=reshape(E,[size(E,1)*size(E,2) size(E,3)]); C=C'; A=C;
C=C(:,linspace(1,size(C,2),size(T,1)));

% calcolo proiezioni
projC = pinv(C)*A;

% figura
N = size(C,2);
M = size(A,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(h_figure); clf;
set(h_figure,'name','Kinetics projections')

if N>0

% diagonal components
subplot(2,4,[1])
plot(t,ones(size(t))*[0:-1: min([-1,-(N-2)]) ],'k--')
hold on
plot(t,C+[0 : -1 : min([-1,-(N-1)]) ],'k','linewidth',2)
hold off
ylim([-(N-1) 1])
yticks((-(N-1):1:0)+0.5)
yticklabels(M:-(N+1):1)
xlabel('Time / fs')
ylabel('Component index')

% projection
subplot(2,4,[2 3 4])
imagesc(projC)
xticks({})
yticks({})
% colormap(spectral(256))
caxis([-1 1])
colorbar('location','manual','position',[0.925 0.585 0.025 0.34])
title('Projection Matrix')

% all the components
subplot(2,4,[6 7 8])
plot(ones(size(t))*(1 : max([1 (size(A,2)-1)]) ),t,'k--')
hold on
plot(transpose(flip(A)+[0:1:(size(A,2)-1)]) , flip(t) ,'k','linewidth',2)
hold off
set(gca,'ydir','reverse')
xticks((1:M)-0.5)
xticklabels({1:M})
ylabel('Time / fs')
xlabel('Component index')

uicontrol('style','text','unit','normalized',...
    'string','Ask Andrea if you don''t know what to do with this plot.',...
    'position',[0.05 0.01 0.9 0.025]);

% subplot(2,4,[5])
% L = reshape(1:M,[N N]); L = L';
% text(0.5,0.5,num2str(L))
% xticks({})
% yticks({})

end

end