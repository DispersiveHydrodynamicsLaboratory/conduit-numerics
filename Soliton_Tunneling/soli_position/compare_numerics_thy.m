% Compare theory to exp't
reformat_on = 1; %reformat numerics into matrix; only needs to be done once per trial

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

if periodic
    dz       = zmax/Nz;    % Spatial  discretization
else
    dz       = zmax/(Nz+1);    % Spatial  discretization
end




data_dir = ['/Volumes/APPM-DHL/data/conduit_eqtn/',...
            '_tmax_',  num2str(round(tmax)),...
            '_zmax_', num2str(round(zmax)),...
            '_Nz_',   num2str(Nz),...
            '_order_',num2str(h),...
            '_init_condns_','soli_tunneling',...
                '_Amax_',num2str(Am),...
                '_asoli_',num2str(asoli),...
                '_hstretch_',num2str(hstretch),...
                '_wave_type_',wave_type,...
            '_bndry_condns_',bc_type,...
            '/'];
load([data_dir,'parameters.mat'],'t','dz','Nz');

% Find maximum time in numerics files
tind = length(t)-1;
for ii=1:length(t)-1
  [fid,foo] = fopen(strcat(data_dir,num2str(ii,'%05d'),'.mat'),'r');
  if fid == -1 % File does not exist
    tind = ii-1;
    disp(['Maximum time = ',num2str(t(tind))]);
    break;
  end
  fclose(fid);
end
tm = tind;

% t- & z-axis
t = t(1:tm);
zplot  = dz*[1:Nz]';

if reformat_on
    A_full = zeros(tind,length(zplot));
    for ti = 1:tm
        load(strcat(data_dir,num2str(ti,'%05d'),'.mat'),'A','tnow','inc');
        A_full(ti,:) = A;
    end
    save([data_dir,'matrix.mat'],'A_full','t','zplot');
else
    load([data_dir,'matrix.mat'],'A_full','t','zplot');
end

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
    
