% Compare theory to exp't
reformat_on = 0; %reformat numerics into matrix; only needs to be done once per trial
path(path,'/Users/appm_admin/Documents/MATLAB/conduit_numerics/solver');
                
%% Numerical Parameters
tmax     = 175;    % Solver will run from t=0 to t=tmax
zmax     = 900;     % Solver will solve on domain z=0 to z=zmax
numout   = round(tmax) ;           % Number of output times
t        = linspace(0,tmax,numout);  % Desired output times
dzinit   = 1/10; % Set to 1/500 for optimum
Nz       = round(zmax/dzinit);
h        = 4;                
periodic = 0;
bc_type = 'time_dependent';


%% IC Parameters
Am     = 2; uplus = Am;
asoli = 5;
hstretch = 5; hstretches = hstretch;
zjump = zmax/2;
z0    = zjump - 50;
zminus = -(zjump - z0);
wave_type = 'r'; % r for RW, d for DSW

% load numerics files
load('sample_numerics.mat','A_full','t','zplot');
tm = 110; % found based on data set
t = t(1:tm);
A_full = A_full(1:tm,:);

% Calculate theoretical soliton position
[zs] = soliton_position_fun(1,asoli,zminus,Am);

[foo,zm] = min(abs(max(zplot)-(zs(t)-zminus+z0)));
figure(1); clf;
    contourf(zplot(1:25:end),t,A_full(:,1:25:end),100,'edgecolor','none');
        cmap = load('CoolWarmFloat257.csv');
        colormap(cmap); 
        xlabel('z'); ylabel('t'); colorbar;
    hold on;
    plot(zs(t(1:zm))-zminus+z0,t(1:zm),'k','LineWidth',2);
    hold off;
    if numel(hstretches)>1
        input('r');
    end
    
% Find soliton position over time and compare to theory
locs = zeros(size(t));
for ti = 1:numel(t)
    [pk,loc] = findpeaks(A_full(ti,:),zplot,'SortStr','descend','NPeaks',1);
    locs(ti) = loc;
end
figure(2); clf;
    plot(zs(t(1:zm))-zminus+z0,t(1:zm),locs,t);


% Check RW
    % Pick random time step in the middle
        tchosen = floor(tm/2);
        Achosen = A_full(tchosen,:);
        figure(3); clf;
            plot(zplot,Achosen);
    % Find expected RW at that timestep
    RW = @(z,t) ones(size(z)) .* (z<=(2*t)) + ...
                z/(2*t)       .* (z>2*t) .* (z<2*uplus*t) + ...
                uplus*ones(size(z)) .* (z>=2*uplus*t);
        figure(3); hold on;
            plot(zplot,RW(zplot-zjump,tchosen),'LineWidth',2);
    
