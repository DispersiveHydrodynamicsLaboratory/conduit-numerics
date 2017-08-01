% Driver script for running the conduit equation solver
% conduit_solver.m
make_matrix_on = 0;
find_peaks_on  = 0;
% Edited on 2017/07/12 by MM & NF
save_on  = 1;  % Set to nonzero if you want to run the solver, set
                % to 0 if you want to plot
periodic = 1; % set to nonzero to run periodic solver (no BCs need)
              % set to 0 to run solver with time-dependent BCs                
plot_on  = 0;  % Set to 1 if you want to plot just before and just
                % after (possibly) calling the solver          
check_IC = 0; % Set to nonzero to plot the ICs and BCs without running the solver

amaxes = 1:1:4;
ws = 100:100:400;

[AMAXES WS] = meshgrid(amaxes,ws);
NS = zeros(size(AMAXES));
AMPLS = struct();

for ii = 1:numel(AMAXES);
    disp(['Run: ',num2str(ii)]);
    % Numerical Parameters
    Amax     = AMAXES(ii);     % Height of the box
    w        = WS(ii);   % Width of the box
    v        = 5*Amax;     % Edges of the box (how smooth the edges are) Possible Radiation Reducer
    tmax     = 500;    % Solver will run from t=0 to t=tmax
%     zmax     = 1000*Amax;     % Solver will solve on domain z=0 to z=zmax
    
    numout   = round(tmax);           % Number of output times
    t        = linspace(0,tmax,numout);  % Desired output times
    if ii >=13
        zmax = 500*Amax;
        dzinit = 1/3;
    else
        zmax     = 1000*Amax;     % Solver will solve on domain z=0 to z=zmax
        dzinit =  1/10; % Spatial Discretization for most accurate runs
                      % With O(h^4), 0.1 gives 10^{-3} max error over t= [0,53]
    end
    z0       = zmax/2;   % Center of the box
    Nz       = round(zmax/dzinit);
    if periodic
        dz       = zmax/Nz;    % Spatial  discretization
    else
        dz       = zmax/(Nz+1);    % Spatial  discretization
    end
        h        = 4   ;           % Order of method used     

    %% PDE Initial and Boundary Conditions
    f = @(z) Amax/2*(tanh((z-z0+w/2)/v)-tanh((z-z0-w/2)/v));
    f = @(z) ones(size(z)) + f(z) .* (f(z)>10^(-4));
    % The initial condition 
        ic_type = '';
        if periodic
            bc_type = 'periodic';
        else
            bc_type = 'time_dependent';
        end

    %% directory run was saved to
    data_dir = ['/Volumes/APPM-DHL/data/conduit_eqtn/',...
                '_tmax_',  num2str(round(tmax)),...
                '_zmax_', num2str(round(zmax)),...
                '_Nz_',   num2str(Nz),...
                '_Amax_', num2str(Amax),...
                '_z0_',   num2str(z0),...
                '_w_',    num2str(w),...
                '_v_',    num2str(v),...
                '_order_',num2str(h),...
                '_init_condns_',ic_type,...
                '_bndry_condns_',bc_type,...
                '/'];
    savefile = sprintf('%sparameters.mat',data_dir);
    load(savefile);
    % Combine all timesteps into one matrix
    if make_matrix_on
        z  = dz*[1:Nz];
        disp('Calculating maximum time increment in saved data files...');
        for ti=1:length(t)+1
            [fid,foo] = fopen(strcat(data_dir,num2str(ii,'%05d.mat')),'r');
            if fid == -1 % File does not exist
                tm = ti-1;
                disp(['Maximum time = ',num2str(t(tm))]);
                break;
            end
            fclose(fid);
        end
        if ti == length(t)+1
            tm = length(t);
            disp(['Maximum time = ',num2str(t(tm))]);
        end
        % Get rid of larger t values
        t = t(1:tm);
        if periodic
            z  = dz:dz:zmax;
        else
            z  = dz:dz:zmax-dz;
        end
        A_full = zeros(length(t)-1,length(z));
        % Load first time step
        load(strcat(data_dir,num2str(0,'%05d')),'A_init');
        for tind=2:tm
            load(strcat(data_dir,num2str(tind,'%05d')),'A','tnow');
            A_full(tind-1,:) = A;
        end
        save([data_dir,'matrix.mat'],'A_full','z','t');
    else
        load([data_dir,'matrix.mat'],'A_full','z','t');
    end
    if find_peaks_on
        allpeaks = struct();
        Nests = [];
        for ti = 1:length(t)-1
            [pks,lcs] = findpeaks(A_full(ti,:),z,'SortStr','descend','MinPeakProminence',5*10^-3);
            allpeaks.(['ind',num2str(ti,'%05d')]).pks = pks;
            allpeaks.(['ind',num2str(ti,'%05d')]).lcs = lcs;
            Nests = [Nests, length(pks)];
        end
        save([data_dir,'peaks.mat'],'allpeaks','Nests');
    else
        load([data_dir,'peaks.mat'],'allpeaks','Nests');
    end
    
%     figure(2); clf;
%         plot(t(1:length(Nests)),Nests,'*'); drawnow; pause(0.2);
        keepsize = 10;
    NS(ii) = mean(Nests(end-keepsize:end));
    amplitudes = zeros(keepsize,round(NS(ii)));
    for jj = 1:keepsize
        try
            amplitudes(jj,:) = allpeaks.(['ind',num2str(round(length(t))-jj-1,'%05d')]).pks(1:length(amplitudes));
        catch
            M = length(allpeaks.(['ind',num2str(round(length(t))-jj-1,'%05d')]).pks);
            amplitudes(jj,:) = [allpeaks.(['ind',num2str(round(length(t))-jj-1,'%05d')]).pks(1:M) nan(1,round(NS(ii))-M)];
        end
                
    end
    AMPLS.(['run',num2str(ii,'%05d')]) =  nanmean(amplitudes);
%     input('r');
end

save([pwd,'/boxdata.mat'],'AMAXES','WS','AMPLS','NS');

	
