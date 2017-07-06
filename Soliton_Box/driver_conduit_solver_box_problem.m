% Driver script for running the conduit equation solver
% conduit_solver.m
save_on  = 1;  % Set to nonzero if you want to run the solver, set
                % to 0 if you want to plot
periodic = 1; % set to nonzero to run periodic solver (no BCs need)
              % set to 0 to run solver with time-dependent BCs                
plot_on  = 0;  % Set to 1 if you want to plot just before and just
                % after (possibly) calling the solver          
check_IC = 0; % Set to nonzero to plot the ICs and BCs without running the solver
% '/Users/dhl/Documents/MATLAB\conduit\data\_tmax_400_zmax_1000_Nz_10000_order_4_init_condns_ampl_1_w_100_asoli_7_bndry_condns_periodic/'
%% PDE Initial and Boundary Conditions
for ampl = [1]
    for w = [100 200]
        for asoli = [3 5 7]
        vstretch = 5;
        f1 = make_smooth_box(vstretch,ampl,w,500);
        Nexp = 12;
        L    = 200;
        if ismember(w,[100 ]) && ismember(asoli,[3 5])
            continue;
        end
        
        if ampl==1 && ismember(w,[100 200]) && asoli==3
            zsoli = 500-w;
            [ zphi,phi ] = conduit_soliton_newton_cg( asoli, Nexp, L, 0 );
        else
            [ zphi,phi ] = conduit_soliton_newton_cg( asoli, Nexp, L, 0 );
            zsoli = 500-2*w;
        end
        
        f2 = @(zq) interp1(zphi+zsoli,phi,zq,'spline',0);
        f = @(z) f1(z) + f2(z) + ones(size(z));
            ic_type = ['ampl_',num2str(ampl),...
                       '_w_',num2str(w),...
                       '_asoli_',num2str(asoli)];
            if periodic
                bc_type = 'periodic';
            else
                bc_type = 'time_dependent';
            end

%% Numerical Parameters
tmax     = 400;    % Solver will run from t=0 to t=tmax
zmax     = 1000;     % Solver will solve on domain z=0 to z=zmax
numout   = round(tmax);           % Number of output times
t        = linspace(0,tmax,numout);  % Desired output times
dzinit   = 0.1;
Nz       = round(zmax/dzinit);
if periodic
    dz       = zmax/Nz;    % Spatial  discretization
else
    dz       = zmax/(Nz+1);    % Spatial  discretization
end
    h        = 4   ;           % Order of method used   
    
    if periodic
        bc_type = 'periodic';
    else
        bc_type = 'time_dependent';
    end

%% Create directory run will be saved to
data_dir = ['/Users/dhl/Documents/MATLAB/conduit/data/',...
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
        end
    end
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
