% Driver script for running the conduit equation solver
% conduit_solver.m
save_on  = 1;  % Set to nonzero if you want to run the solver, set
                % to 0 if you want to plot
periodic = 0; % set to nonzero to run periodic solver (no BCs need)
              % set to 0 to run solver with time-dependent BCs                
plot_on  = 1;  % Set to 1 if you want to plot just before and just
                % after (possibly) calling the solver          
check_IC = 0; % Set to nonzero to plot the ICs and BCs without running the solver

%% Numerical Parameters
% for zmax = [1000 2000]
tmax     = 200;    % Solver will run from t=0 to t=tmax
zmax     = 250;     % Solver will solve on domain z=0 to z=zmax
numout   = round(tmax);           % Number of output times
t        = linspace(0,tmax,numout);  % Desired output times
dzinit =  1/100; % Spatial Discretization: for most accurate runs
                  % With O(h^4), 0.1 gives 10^{-3} max error over t= [0,53]
Nz       = round(zmax/dzinit);
    h        = 4   ;           % Order of method used     

if periodic
    dz       = zmax/Nz;    % Spatial  discretization
else
    dz       = zmax/(Nz+1);    % Spatial  discretization
end

%% PDE Initial and Boundary Conditions
% load('output.mat');
% f = @(z) interp1(domain,area,z,'spline',1);
%     ic_type = 'soligas';
f = @(z) ones(size(z)) + 0.5*sin((2*pi)*z/(25)).*(z<125);
 g0 = @(t) ones(size(t));
dg0 = @(t) zeros(size(t));
 g1 = @(t) ones(size(t));
dg1 = @(t) zeros(size(t));
ic_type = '_sin_example_';

    if periodic
        bc_type = 'periodic';
    else
        bc_type = 'time_dependent';
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
      zplot  = dz*[1:Nz];
      tplot  = linspace(0,tmax,floor(tmax*10));
      A_init = f(zplot);
    if plot_on
        % Plot initial conditions and boundary conditions
        fontsize = 12;
        figure(1); clf;
        subplot(2,1,1);
            plot(zplot,A_init); xlabel('z'); ylabel('A(z)'); 
            title('Initial Conditions');
        subplot(2,1,2);
            plot(tplot,g0(tplot),tplot,g1(tplot));
            xlabel('t'); ylabel('A(t)'); title('Boundary Conditions');
            legend('At z=0','At z=zmax');
        set(gca,'fontsize',fontsize,'fontname','times');
        
        if check_IC
            return;
        else
            input('Return to continue');
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

%% If chosen, plot data associated with the parameters and conditions above
if plot_on
    disp('Calculating maximum time increment in saved data files...');
    for ii=1:length(t)+1
        [fid,foo] = fopen(strcat(data_dir,num2str(ii,'%05d.mat')),'r');
        if fid == -1 % File does not exist
            tm = ii-1;
            disp(['Maximum time = ',num2str(t(tm))]);
            break;
        end
        fclose(fid);
    end
    if ii == length(t)
        disp(['Maximum time = ',num2str(t(tm))]);
    end
    % Get rid of larger t values
    t = t(1:tm);
    if periodic
        zplot  = dz:dz:zmax;
    else
        zplot  = dz:dz:zmax-dz;
    end
    A_full = zeros(length(t)-1,length(zplot));
    % Load first time step
    load(strcat(data_dir,num2str(0,'%05d')),'A_init');
        fontsize = 12;
        fig=figure(2); clf;
        plot(zplot,A_init);
        qaxis = axis;
        title(['Time: ',num2str(t(1))]);
        set(gca,'fontsize',fontsize,'fontname','times');
        drawnow
%         input('Return');
    %Plot subsequent time steps
    for tind=2:tm
        load(strcat(data_dir,num2str(tind,'%05d')),'A','tnow');
        fontsize = 12;
        fig=figure(2); clf;
        plot(zplot,A);
        A_full(tind-1,:) = A;
        axis(qaxis)
        hold off;
        title(['Time: ',num2str(tnow)]);
        set(gca,'fontsize',fontsize,'fontname','times');
        drawnow;
    end

end
