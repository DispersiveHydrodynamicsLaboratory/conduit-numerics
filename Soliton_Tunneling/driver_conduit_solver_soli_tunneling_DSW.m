% Driver script for running the conduit equation solver
% conduit_solver.m


% save_on  = 0;  % Set to nonzero if you want to run the solver, set
%                 % to 0 if you want to plot
% periodic = 0; % set to nonzero to run periodic solver (no BCs need)
%               % set to 0 to run solver with time-dependent BCs      
% check_IC = 0; % Set to 1 to only plot ICs
% plot_on  = 1;  % Set to 1 if you want to plot just before and just
%                 % after (possibly) calling the solver

% for ind = 1%:length(aDSW)            
%     disp(['Soliton pair: ',num2str(aDSW(ind)),' ',num2str(aRW(ind))]);

%% PDE Initial and Boundary Conditions
[fRW, f] = soli_tunneling_IC( ADSW, aDSW(ind), aRW(ind), zmax);
% ICs and BCs consistent with DSW generation
    g0      = @(t) ADSW*ones(size(t));
    dg0     = @(t) zeros(size(t));  % derivative of BCs
    g1      = @(t) ones(size(t)) ;  % BCs at z=zmax
    dg1     = @(t) zeros(size(t));  % derivative of BCs
    ic_type = ['soli_tunneling_DSW_ADSW_',num2str(ADSW),...
               '_aDSW_',num2str(aDSW(ind)),...
               '_aRW_' ,num2str(aRW(ind))];
           
    if periodic
        bc_type = 'periodic';
    else
        bc_type = 'time_dependent';
    end
   

%% Create directory run will be saved to
data_dir = ['/Volumes/Data Storage/Numerics/conduit_eqtn/',...
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
% end

% If chosen, plot data associated with the parameters and conditions above
if plot_on
    plot_data_fun(data_dir);
end