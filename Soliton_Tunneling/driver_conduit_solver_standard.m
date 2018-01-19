% Driver script for running the conduit equation solver
% conduit_solver.m
save_on  = 1;  % Set to nonzero if you want to run the solver, set
                % to 0 if you want to plot
periodic = 0; % set to nonzero to run periodic solver (no BCs need)
              % set to 0 to run solver with time-dependent BCs                
plot_on  = 0;  % Set to 1 if you want to plot just before and just
                % after (possibly) calling the solver          
check_IC = 0; % Set to nonzero to plot the ICs and BCs without running the solver

%% Numerical Parameters
tmax     = 700;    % Solver will run from t=0 to t=tmax
zmax     = 6500;     % Solver will solve on domain z=0 to z=zmax
numout   = round(tmax);           % Number of output times
t        = linspace(0,tmax,numout);  % Desired output times
dzinit =  1/3; % Spatial Discretization for most accurate runs
                  % With O(h^4), 0.1 gives 10^{-3} max error over t= [0,53]
Nz       = round(zmax/dzinit);
if periodic
    dz       = zmax/Nz;    % Spatial  discretization
else
    dz       = zmax/(Nz+1);    % Spatial  discretization
end
    h        = 2   ;           % Order of method used     

%% PDE Initial and Boundary Conditions
z0soli = 50; z0RW = 275; z0bump = z0RW;
ARW = 1.25; vs = (ARW-1)/2; hsRW = 50;
Abump = ARW + 1; hsbump = 300;
asoli  = 13;
[ zsoli,psoli ] = conduit_soliton_newton_cg( asoli, zmax, 0 );
zsoli = zsoli + z0soli;
fsoli = @(z) interp1(zsoli,psoli,z,'spline',0);
RW       = @(z) (1+vs) + vs*tanh(1/hsRW*(z-z0RW));
% hump     = @(z) (ARW+1)*sech(1/(1*hs)*(z-z0sech));
zbump   = ((-hsbump/2+1/25):1/25:hsbump/2) + z0bump;
sinbump =  Abump*(0.5+0.5*(sin(linspace(-pi/2,pi*3/2,25*hsbump))));
bump = @(z) interp1(zbump,sinbump,z,'spline',0);
f = @(z) fsoli(z) + RW(z) + bump(z);

 g0 = @(t) ones(size(t));
dg0 = @(t) zeros(size(t));
 g1 = @(t) ARW*ones(size(t));
dg1 = @(t) zeros(size(t));
    ic_type = '';
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
      zplot  = dz*[1:Nz];
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
            legend(ic_type);
            drawnow; 
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
plot_data_fun(data_dir);
print('exm','-dpng','-r600')
send_mail_message('mdmaide2','Matlab is finished!','Your sims are done!','exm.png')