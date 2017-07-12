% Driver script for running the conduit equation solver
% conduit_solver.m

%% Numerical Parameters (not defined in wrapped script)
numout   = round(tmax) ;           % Number of output times
t        = linspace(0,tmax,numout);  % Desired output times
dzinit   = 1/100; % Set to 1/500 for optimum
Nz       = round(zmax/dzinit);
h        = 4;
    if periodic
        dz       = zmax/Nz;    % Spatial  discretization
    else
        dz       = zmax/(Nz+1);    % Spatial  discretization
    end


%% PDE Initial and Boundary Conditions
[f] = soli_tunneling_IC_exm( m_plus, aDSW, zmax);
% ICs and BCs consistent with DSW generation
    g0      = @(t) m_plus*ones(size(t));
    dg0     = @(t) zeros(size(t));  % derivative of BCs
    g1      = @(t) ones(size(t)) ;  % BCs at z=zmax
    dg1     = @(t) zeros(size(t));  % derivative of BCs
    ic_type = ['soli_tunneling_DSW_ADSW_',num2str(m_plus),...
               '_aDSW_',num2str(aDSW)];
           
    if periodic
        bc_type = 'periodic';
    else
        bc_type = 'time_dependent';
    end
   

%% Create directory run will be saved to
data_dir = ['./conduit_eqtn/',...
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