save_on  = 1;  % Set to nonzero if you want to run the solver, set
                % to 0 if you want to plot
periodic = 0; % set to nonzero to run periodic solver (no BCs need)
              % set to 0 to run solver with time-dependent BCs      
check_IC = 0; % Set to 1 to ONLY plot ICs (and BCs, if applicable)
              % Set to 0 to run the solver
plot_on  = 1;  % Set to 1 if you want to plot just before and just
                % after (possibly) calling the solver
%% Numerical Parameters
tmax     = 50;    % Solver will run from t=0 to t=tmax
zmax     = 100;    % Solver will solve on domain z=0 to z=zmax
numout   = 200;     % Number of output times
t        = linspace(0,tmax,numout);  % Desired output times
dzinit   = .05;    % Set to 1/500 for optimum
Nz       = round(zmax/dzinit); % Works on its own
h        = 4;      % Order of the scheme; can use 2 or 4

if periodic
    dz       = zmax/Nz;    % Spatial  discretization
else
    dz       = zmax/(Nz+1);    % Spatial  discretization
end
% load('/Volumes/Data Storage/Experiments/Photos/2016_05_17/quantities.mat');

%% PDE Initial and Boundary Conditions
% ICs and BCs consistent with DSW generation
z0 = 50;
Aplus = 2; Aminus = 1;
t0 = z0/(2*Aplus);
% deltat = dzinit/100;
z = 0; %0:1/100:z0;
    g0      = @(t) BC_DSW( t, t0 , z, z0, Aplus , Aminus);         % BCs at z=0
    dg0     = @(t) BC_DSW_deriv( t, t0 , z, z0, Aplus , Aminus);   % derivative of BC
    g1      = @(t) ones(size(t)) ;              % BCs at z=zmax
    dg1     = @(t) zeros(size(t));              % derivative of BC
    f       = @(z) ones(size(z));               % ICs at t=0
    ic_type = 'constant_fun';

    if periodic
        bc_type = 'periodic';
    else
        bc_type = ['time_dependent_DSW_Aplus_',num2str(Aplus),'_Aminus_',num2str(Aminus)];
    end
    

%% Create directory run will be saved to
data_dir = ['./data/conduit_eqtn',...
            '_tmax_',  num2str(round(tmax)),...
            '_zmax_', num2str(round(zmax)),...
            '_Nz_',   num2str(Nz),...
            '_order_',num2str(h),...
            '_init_condns_',ic_type,...
            '_bndry_condns_',bc_type,...
            '/'];
zplot = dz:dz:zmax-dz;

for ii=40:length(t)-1
    load(strcat(data_dir,num2str(ii,'%05d.mat')));
    [pks,locs] = findpeaks(A,'MinPeakProminence',3);
    if isempty(pks)==0
        Acut = A(1:locs);
        if sum(Acut(:)<2)>0
            break;
        end
    end
    figure(19);
    plot(zplot,A);
    axis([0 zmax 0.9*Aminus 3]);
    if ~isempty(pks)
        hold on;
        plot(zplot(locs),pks,'*');
        hold off;
    end
    title(['t=',num2str(tnow)]);
    drawnow; pause(0.05);
end

%     plot(zplot,A);
%     title(['t=',num2str(tnow)]);
%     pause(0.5);

