% Driver script for running the conduit equation solver
% conduit_solver.m

save_on  = 1;  % Set to nonzero if you want to run the solver, set
                % to 0 if you want to plot
periodic = 0; % set to nonzero to run periodic solver (no BCs need)
              % set to 0 to run solver with time-dependent BCs      
check_IC = 0; % Set to 1 to ONLY plot ICs (and BCs, if applicable)
              % Set to 0 to run the solver
plot_on  = 1;  % Set to 1 if you want to plot just before and just
                % after (possibly) calling the solver

% Directory where data is presently saved                
main_dir = '/Users/appm_admin/Documents/MATLAB/conduit_numerics_backup/DSW_runs';

Aplus_opts = [4];
z0_opts = [93.4079 125.3071 164.6100];
DSW_opts = 1:(length(Aplus_opts)*length(z0_opts));
for ii = DSW_opts
    Aplus = Aplus_opts(ceil(ii/5));
    zind = mod(ii,3);
    if zind==0
        zind = 3;
    end
    z0 = z0_opts(zind);
        t0 = z0/(2*Aplus);
        if Aplus>=4
            toffset = 10;
        else
            toffset = 5;
        end
    %% Numerical Parameters
    tmax     = round(t0*(toffset+1)+10);    % Solver will run from t=0 to t=tmax
    zmax     = z0+100;    % Solver will solve on domain z=0 to z=zmax
    numout   = round(tmax)*4+1;     % Number of output times
    t        = linspace(0,tmax,numout);  % Desired output times
    dzinit   = 1/250;    % Set to 1/500 for optimum
    Nz       = round(zmax/dzinit); % Works on its own
    h        = 4;      % Order of the scheme; can use 2 or 4

    if periodic
        dz       = zmax/Nz;    % Spatial  discretization
    else
        dz       = zmax/(Nz+1);    % Spatial  discretization
    end


%% PDE Initial and Boundary Conditions
% ICs and BCs consistent with DSW generation
Aminus = 1;

% deltat = dzinit/100;
z = 0; %0:1/100:z0;
tvec        = linspace(0,tmax,floor(tmax*100));

[A,Adot]    = BC_DSW(tvec, t0 , z, z0, Aplus , Aminus, toffset);   
    g0      = @(t) interp1(tvec,A   ,t,'spline','extrap');  % BCs at z=0
    dg0     = @(t) interp1(tvec,Adot,t,'spline','extrap');   % derivative of BC
    g1      = @(t) ones(size(t)) ;              % BCs at z=zmax
    dg1     = @(t) zeros(size(t));              % derivative of BC
    f       = @(z) ones(size(z));               % ICs at t=0
    ic_type = 'constant_fun';

    if periodic
        bc_type = 'periodic';
    else
        bc_type = ['time_dependent_DSW_Aplus_',num2str(Aplus),'_Aminus_',num2str(Aminus),'_z0_',num2str(z0)];
    end

%% Create directory run will be saved to
data_dir = [main_dir,'/data_for_Dalton/conduit_eqtn',...
            '_tmax_',  num2str(round(tmax)),...
            '_zmax_', num2str(round(zmax)),...
            '_Nz_',   num2str(Nz),...
            '_order_',num2str(h),...
            '_init_condns_',ic_type,...
            '_bndry_condns_',bc_type,...
            '/'];
% Create the data directory if necessary
disp(data_dir);
if ~exist(data_dir,'dir')
    mkdir(data_dir);
else
    disp(['Warning, directory ',data_dir]);
    disp('already exists, possibly overwriting data');
%     input('Return to continue. Cancel to not');
end

savefile = sprintf('%sparameters.mat',data_dir);

%% If chosen, run the solver using the parameters and conditions above
if save_on
    % Load initial data
      zplot  = dz*[1:Nz];
      tplot  = linspace(0,tmax,floor(tmax*100));
      A_init = f(zplot);
    if plot_on
        % Plot initial conditions and boundary conditions
        fontsize = 12;
        figure(1); clf;
        subplot(3,1,1);
            plot(zplot,A_init); xlabel('z'); ylabel('A(z)'); 
            title('Initial Conditions');
        subplot(3,1,2);
            plot(tplot,g0(tplot),tplot,g1(tplot));
            xlabel('t'); ylabel('A(t)'); title('Boundary Conditions');
            legend('At z=0','At z=zmax');
        subplot(3,1,3);
            plot(tplot,dg0(tplot),tplot,dg1(tplot));
            xlabel('t'); ylabel('A(t)'); title('Derivative of Boundary Conditions');
            legend('At z=0','At z=zmax');
        set(gca,'fontsize',fontsize,'fontname','times');
        pause(1); drawnow;
        if check_IC
            continue;
        end
    end
    
    if periodic
    % Save parameters
        save(savefile,'t','Nz','dz','zmax','f','periodic');
    % Run timestepper
        conduit_solver_periodic( t, zmax, Nz, h, f, data_dir );      
    else    
    % Save parameters
        save(savefile,'t','Nz','dz','zmax','g0','dg0','g1','dg1','f','periodic');
    % Run timestepper
        conduit_solver( t, zmax, Nz, h, g0, dg0, g1, dg1, f, data_dir );
    end
else
    load(savefile);
end

if plot_on
    plot_data_fun(data_dir,4,6,floor((toffset)*5*t0),floor((toffset-2)*5*t0));
end

end

