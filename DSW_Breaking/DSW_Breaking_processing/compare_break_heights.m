%% Load Numerics Results
main_dir = '/Users/appm_admin/Documents/MATLAB/conduit_numerics_backup/DSW_runs/data_for_Dalton/';
load([main_dir,'processed_data.mat'],'all_z0s_in','all_z0s_out','all_t0s_in','all_t0s_out','z0_opts','Abacks');%'all_Aminus','all_breaktimes','all_breakheights',

Aback = Abacks;
all_Aminus = ones(size(all_z0s_in));
all_breakheights = all_z0s_in;

%% NONDIMENSIONAL THINGS
% Compare to Expected Break Heights
% Initialize break height lines to compare to
    nn = length(z0_opts);
    A0 = repmat([Aback],nn,1);
    % Commands to generate legend
        zb = [];
        leg = [];
        for ii = 1:nn
            zb = [zb; [z0_opts(ii) z0_opts(ii)] ];
            leg = [leg;{['z0: ',num2str(z0_opts(ii))]}];
        end
% Plot Comparison (ND)
    figure(3); clf; 
        h = plot(A0',all_z0s_out','-*',A0',all_z0s_in','--');
        set(gca,'FontSize',14,'FontName','times');
        for ii = 1:nn
            set(h(ii+nn),'Color',get(h(ii),'Color'));
        end
    axis([min(Aback) max(Aback),0, 1.05*max(z0_opts)]);
    xlabel('Jump Height (ND)'); ylabel('Breaking Height'); title('Numerical DSWs')
    legend(leg);
    
%% Send data to 'dimensional' space
load('/Users/appm_admin/Documents/MATLAB/Dalton/2017_04_01_Dalton_WaveBreaking/quantities.mat','T0fac','U0fac','L0fac','alpha','epsilon');
% load('/Volumes/Data Storage/Experiments/data_processing/2016_05_17_Dalton_Glycerol/RESULTS.mat','alpham');
alpham = alpha;
Q0 = [0.25]; % mL/min
D0 = alpham*Q0.^(1/4);
% %    R = L0fac*Diam * R_nd
% %    Z = L0fac*Diam/epsilon^0.5 * Z_nd
    L0 =  0.157363848283770; %found experimentally
% %    vel = U0fac*Diam^2 * vel_nd
    U0 =  0.080776104456255; %found experimentally
% %     t = T0fac/Diam/epsilon^0.5 * t_nd
    T0 = 1.950756251929821; %found experimentally

% Convert to set_quantitites values
    %% Check 'dimensionality'
    a = load('/Users/appm_admin/Documents/MATLAB/conduit_numerics_backup/DSW_runs/data_for_Dalton/conduit_eqtn_tmax_150_zmax_193_Nz_48352_order_4_init_condns_constant_fun_bndry_condns_time_dependent_DSW_Aplus_2_Aminus_1_z0_93.4079/conduit_edges.mat');
    b = load('/Users/appm_admin/Documents/MATLAB/conduit_numerics_backup/DSW_runs/data_for_Dalton/conduit_eqtn_tmax_150_zmax_193_Nz_48352_order_4_init_condns_constant_fun_bndry_condns_time_dependent_DSW_Aplus_2_Aminus_1_z0_93.4079/parameters.mat');
    zold = b.dz:b.dz:b.zmax-b.dz;
    znew = linspace(b.dz,b.zmax-1,300);
    [Zold,Told] = meshgrid(zold,b.t);
    [Znew,Tnew] = meshgrid(znew,b.t);
    Zold = Zold*L0;
    Znew = Znew*L0;
    Told = Told*T0;
    Tnew = Tnew*T0;
    A    = interp2(Zold,Told,a.area_mat,Znew,Tnew);
%     figure(1); clf;
%     contourf(Znew,Tnew,A,'edgecolor','none');
%             cmap = load('CoolWarmFloat257.csv');
%             colormap(cmap); colorbar;
%             xlabel('z'); ylabel('t');
%             hold on;
%             plot([20 20]*L0,[0,max(b.t)]/T0);
            
    % Convert to Dimensional Quantities
        z0_opts_25        = z0_opts*L0;
        z0_expects_dim_25 = all_z0s_in*L0;
        z0_actuals_dim_25 = all_z0s_out*L0;
%         t0_opts_25        = t0_opts*T0;
        t0_expects_dim_25 = all_t0s_in*T0;
        t0_actuals_dim_25 = all_t0s_out*T0;

%% FLOW RATE 0.25 
% Breaking Height Comparison
% Initialize break height lines to compare to
    nn = length(z0_opts_25);
    A0 = repmat([Abacks],nn,1);
    % Commands to generate legend
        zb = [];
        leg = [];
        for ii = 1:nn
            zb = [zb; [z0_opts_25(ii) z0_opts_25(ii)] ];
            leg = [leg;{['z0: ',num2str(round(z0_opts_25(ii)))]}];
        end
    figure(4); clf; 
        h = plot(A0',z0_actuals_dim_25','-*',A0',z0_expects_dim_25','--');
        set(gca,'FontSize',14,'FontName','times');
        for ii = 1:nn
            set(h(ii+nn),'Color',get(h(ii),'Color'));
        end
    axis([min(Aback) max(Aback),0, 1.05*max(z0_opts_25)]);
    xlabel('Jump Height (ND)'); ylabel('Breaking Height (cm)'); title('Numerical DSWs (Dim, Q0=0.25mL/min)')
    legend(leg);
% Breaking time Comparison
% Initialize break height lines to compare to
    nn = length(z0_opts_25);
    A0 = repmat([Abacks],nn,1);
    % Commands to generate legend
        zb = [];
        leg = [];
        for ii = 1:nn
            zb = [zb; [z0_opts_25(ii) z0_opts_25(ii)] ];
            leg = [leg;{['z0: ',num2str(round(z0_opts_25(ii)))]}];
        end
    figure(5); clf; 
        h = plot(A0',t0_actuals_dim_25','-*',A0',t0_expects_dim_25','--');
        set(gca,'FontSize',14,'FontName','times');
        for ii = 1:nn
            set(h(ii+nn),'Color',get(h(ii),'Color'));
        end
    axis([min(Aback) max(Aback),0, 1.05*max(t0_expects_dim_25(:))]);
    xlabel('Jump Height (ND)'); ylabel('Breaking Time (s)'); title('Numerical DSWs (Dim, Q0=0.25mL/min)')
    legend(leg);
    
% Save to compare with experiments
save('numerical_results.mat','z0_expects_dim_25','z0_actuals_dim_25',...
                             't0_expects_dim_25','t0_actuals_dim_25',...
                              'A0');
                          
save('/Users/appm_admin/Documents/MATLAB/Dalton/2017_04_01_Dalton_WaveBreaking/Processing/numerical_results.mat',...
            'z0_expects_dim_25','z0_actuals_dim_25',...
            't0_expects_dim_25','t0_actuals_dim_25',...
            'A0');
                          
%     % Convert to Dimensional Quantities
%         z0_opts_25        = z0_opts*L0;
%         z0_expects_dim_25 = all_z0s_in*L0;
%         z0_actuals_dim_25 = all_z0s_out*L0;
% %         t0_opts_25        = t0_opts*T0;
%         t0_expects_dim_25 = all_t0s_in*T0;
%         t0_actuals_dim_25 = all_t0s_out*T0;
    