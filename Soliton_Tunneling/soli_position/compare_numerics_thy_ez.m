% Compare theory to exp't
reformat_on = 1; %reformat numerics into matrix; only needs to be done once per trial
path(path,'/Users/appm_admin/Documents/MATLAB/conduit_numerics/solver');

% for various hstretch
hstretches = 5;%[5 2 1 0.1 0.01];                
for hstretch = hstretches
                
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
Am     = 2;
asoli = 5;
% hstretch = defined in loop;
zjump = zmax/2;
z0    = zjump - 50;
zminus = -(zjump - z0);
wave_type = 'r'; % r for RW, d for DSW

% load numerics files
load('sample_numerics.mat','A_full','t','zplot');

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
end
    
