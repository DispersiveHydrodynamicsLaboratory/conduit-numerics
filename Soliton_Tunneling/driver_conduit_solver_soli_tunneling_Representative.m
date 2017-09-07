% Driver script for running the conduit equation solver
% conduit_solver.m

%% MAJOR EDIT: parameters file includes ALL parameters
%% makes processing easier

save_on  = 1;  % Set to nonzero if you want to run the solver, set
                % to 0 if you want to plot
periodic = 0; % set to nonzero to run periodic solver (no BCs need)
              % set to 0 to run solver with time-dependent BCs      
check_IC = 0; % Set to 1 to only plot ICs
plot_on  = 0;  % Set to 1 if you want to plot just before and just
               % after (possibly) calling the solver
                                
%% Numerical Parameters
tmax     = 400;    % Solver will run from t=0 to t=tmax
zmax     = 2000;     % Solver will solve on domain z=0 to z=zmax
numout   = round(tmax) ;           % Number of output times
t        = linspace(0,tmax,numout);  % Desired output times
dzinit   = 1/4; % Set to 1/500 for optimum
Nz       = round(zmax/dzinit);
h        = 4;                

%% IC Parameters
Am     = 2;
asoli  = 7;
hstretch = 5;
zjump = 200;
z0    = zjump - 75;
wave_type = 'r'; % r for RW, d for DSW

if periodic
    dz       = zmax/Nz;    % Spatial  discretization
else
    dz       = zmax/(Nz+1);    % Spatial  discretization
end


%% PDE Initial and Boundary Conditions
% ICs and BCs consistent with soli-wave generation
[f] = soli_tunneling_IC_varied( Am, asoli, zmax, hstretch, zjump, z0, wave_type);
    if strcmp(wave_type,'d')
        g0      = @(t) Am*ones(size(t));
        g1      = @(t) ones(size(t)) ;  % BCs at z=zmax
    else % wave_type = 'r'
        g0      = @(t) ones(size(t));
        g1      = @(t) Am*ones(size(t)) ;  % BCs at z=zmax
    end
    
    dg0     = @(t) zeros(size(t));  % derivative of BCs
    dg1     = @(t) zeros(size(t));  % derivative of BCs
    ic_type = ['soli_tunneling_Amax_',num2str(Am),...
               '_asoli_',num2str(asoli),...
               '_hstretch_',num2str(hstretch),...
               '_wave_type_',wave_type,'_trapping'];

    if periodic
        bc_type = 'periodic';
    else
        bc_type = 'time_dependent';
    end
   

%% Create directory run will be saved to
data_dir = ['/Volumes/APPM-DHL/data/conduit_eqtn/',...
            '_tmax_',  num2str(round(tmax)),...
            '_zmax_', num2str(round(zmax)),...
            '_Nz_',   num2str(Nz),...
            '_order_',num2str(h),...
            '_init_condns_',ic_type,...
            '_bndry_condns_',bc_type,...
            '/'];
% Create the data directory if necessary
if ~exist(data_dir,'dir')
    mkdir(data_dir);
else
    disp(['Warning, directory ',data_dir]);
    disp('already exists, possibly overwriting data');
end

savefile = sprintf('%sparameters.mat',data_dir);

%% If chosen, run the solver using the parameters and conditions above
if save_on
    % Load initial data
      zplot  = dz*[1:Nz]';
      tplot  = linspace(0,tmax,floor(tmax*10));
      A_init = f(zplot);
    if plot_on
        % Plot initial conditions and boundary conditions
        fontsize = 12;
        figure(1); clf;
        if ~periodic
            subplot(3,1,1);
        end
            plot(zplot,A_init); xlabel('z'); ylabel('A(z)'); 
            title('Initial Conditions');
        if ~periodic
        subplot(3,1,2);
            plot(tplot,g0(tplot),tplot,g1(tplot));
            xlabel('t'); ylabel('A(t)'); title('Boundary Conditions');
            legend('At z=0','At z=zmax');
        subplot(3,1,3);
            plot(tplot,dg0(tplot),tplot,dg1(tplot));
            xlabel('t'); ylabel('A(t)'); title('Derivs of Boundary Conditions');
            legend('At z=0','At z=zmax');
        end
        set(gca,'fontsize',fontsize,'fontname','times');
        pause(0.25);
        if check_IC
            return;
        end
    end
    
    if periodic
    % Save parameters
%         save(savefile,'t','Nz','dz','zmax','f','periodic');
    save(savefile);
    % Run timestepper
        conduit_solver_periodic( t, zmax, Nz, h, f, data_dir );      
    else    
    % Save parameters
%         save(savefile,'t','Nz','dz','zmax','g0','dg0','g1','dg1','f','periodic');
    save(savefile);
    % Run timestepper
        conduit_solver( t, zmax, Nz, h, g0, dg0, g1, dg1, f, data_dir );
    end
else
    load(savefile);
end
% end

% If chosen, plot data associated with the parameters and conditions above
if plot_on
    plot_data_fun(data_dir);
end
end