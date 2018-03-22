save_on = 1; % 1 = run data processing, 0 = just plot
reformat = 0; % set to 1 to reformat numerics files as .mat files 
              % only needs to be done once (when new trials are added)
plot_debug_on = 0; % Doesn't plot or show output during processing if off
main_dir = '/Volumes/Data Storage/Numerics/conduit_eqtn/DSW_Breaking/';
% Inputs: jump ratios and break heights (non-dimensional)
load('fig_quants.mat');
Abacks = 2:0.5:4;
z0_opts = [in([1 6 11]).zb_fitND];

Afront  = ones(size(Abacks));
if save_on
    % set_quantities;
    % Preallocate data matrices
    all_z0s_in       = zeros(length(z0_opts),length(Abacks));
    all_t0s_in       = zeros(length(z0_opts),length(Abacks));
    all_Abacks       = zeros(length(z0_opts),length(Abacks));
    all_Afronts      = zeros(length(z0_opts),length(Abacks));
    all_t0s_out      = zeros(length(z0_opts),length(Abacks));
    all_z0s_out      = zeros(length(z0_opts),length(Abacks));

    % Process data; find breaking height (based on experimental code)
    for wi = 1:length(z0_opts)
        z0   = z0_opts(wi);
        trials = 1:length(Abacks);
        breaktime = zeros(size(trials));
        breakheight_ND = zeros(size(trials));

    for ii = trials
        disp(['Jump Height: ',num2str(Abacks(ii)),' z0: ',num2str(z0)]);
        t0 = z0/(2*Abacks(ii)); % Expected breaking time 
        if Abacks(ii)>=4
            toffset = 10;
        else
            toffset = 5;
        end
        %% Numerical Parameters
            tmax     = round(t0*(toffset+1)+10);    % Solver will run from t=0 to t=tmax
            zmax     = z0+100;    % Solver will solve on domain z=0 to z=zmax
            numout   = round(tmax)*4+1;     % Number of output times
            t        = linspace(0,tmax,numout);  % Desired output times
            dzinit   = 1/100;    % Set to 1/500 for optimum
            Nz       = round(zmax/dzinit); % Works on its own

            data_dir = [main_dir,'/data_for_Dalton/conduit_eqtn',...
                       '_tmax_',  num2str(round(tmax)),...
                       '_zmax_', num2str(round(zmax)),...
                       '_Nz_',   num2str(Nz),...
                       '_order_4_init_condns_constant_fun_',...
                       'bndry_condns_time_dependent_DSW_',...
                       'Aplus_',num2str(Abacks(ii)),...
                       '_Aminus_1_z0_',num2str(z0),'/'];

        loadfile = sprintf('%sparameters.mat',data_dir);
        load(loadfile);

        if reformat
            % Reformat numerics files into one file of areas
            numfile = reformat_numerics_files(data_dir);
        else
            % Name of reformatted file
            numfile = [data_dir, 'conduit_edges.mat'];
        end
        load(numfile,'area_mat');
        A = area_mat; clear('area_mat');
        Amax = max(A(:));
        [m,n] = size(A);
        numpics = m;
        tm      = m;
        t       = t(1:tm);
        z = dz:dz:zmax-dz;
        
        %% Calculate breaking height/time
        [Abreak] = find_breaking(A,z,t,Abacks(ii),Afront(ii),z0,t0,toffset,0, plot_debug_on, 1);
        disp(['Found Breaking Height: ',num2str(Abreak.height)]);
        % Save data 
        all_z0s_in(wi,ii)   = z0;
        all_z0s_out(wi,ii)  = Abreak.height;
        all_t0s_in(wi,ii)   = t0;             
        all_t0s_out(wi,ii)  = Abreak.time/toffset;
        all_Abacks(wi,ii)    = Abacks(ii);
        all_Afronts(wi,ii)   = Afront(ii);
        
    end
    end
    save([main_dir,'processed_data.mat'],'all_z0s_in','all_z0s_out','all_t0s_in','all_t0s_out','all_Abacks','all_Afronts','z0_opts','Abacks');
else
    load([main_dir,'processed_data.mat'],'all_z0s_in','all_z0s_out','all_t0s_in','all_t0s_out','all_Abacks','all_Afronts','z0_opts','Abacks');
end

% if length(Aback)==1
%     return;
% end

% Compare to Expected Break Heights
% Initialize break height lines to compare to
    nn = length(z0_opts);
    A0 = repmat([min(Abacks) max(Abacks)],nn,1);
    % Commands to generate legend
        zb = [];
        legzb = [];
        for ii = 1:nn
            zb = [zb; [z0_opts(ii) z0_opts(ii)] ];
            legzb = [legzb;{['z0: ',num2str(z0_opts(ii))]}];
        end
% Plot Comparison
    figure(3); clf; 
        h = plot(all_Abacks',all_z0s_out','-*',all_Abacks',all_z0s_in','--');
        set(gca,'FontSize',14,'FontName','times');
        for ii = 1:nn
            set(h(ii+nn),'Color',get(h(ii),'Color'));
        end
    axis([min(Abacks) max(Abacks),0, 1.05*max(z0_opts)]);
    xlabel('Jump Height (ND)'); ylabel('Breaking Height'); title('Numerical DSWs')
    legend(legzb);
    
% Plot Comparison
    t0fun = @(Aback,z0_opts) z0_opts./(2*Aback);
    Abackvec = 2:0.5:4;
    plot_t0_vec = [];
    for ii = 1:length(z0_opts)
        plot_t0_vec = [plot_t0_vec; t0fun(Abackvec,z0_opts(ii))];
    end
    
    figure(4); clf; 
        h = plot(all_Abacks',all_t0s_out','-*',repmat(Abackvec,length(z0_opts),1)',plot_t0_vec','--');
        set(gca,'FontSize',14,'FontName','times');
        for ii = 1:nn
            set(h(ii+nn),'Color',get(h(ii),'Color'));
        end
    axis([min(Abacks) max(Abacks),0, 1.05*max([all_t0s_in(:); plot_t0_vec(:)])]);
    xlabel('Jump Height (ND)'); ylabel('Breaking Time'); title('Numerical DSWs')
    legend(legzb);

