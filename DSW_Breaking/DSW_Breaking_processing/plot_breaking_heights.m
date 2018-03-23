num_dir = '/Users/appm_admin/Documents/MATLAB/conduit_numerics/DSW_Breaking/DSW_Breaking_processing/';
load([num_dir,'fig_quants.mat']);
load([num_dir,'fit_params.mat']);
load([num_dir,'numerical_results.mat']);

%% Rename variables
zbin  = [dsw.zb];
zbfit = [dsw.zb_fit];
zbND  = [dsw.zb_fitND];
zbout = [out.zb];

tbin  = [dsw.tb];
tbfit = [dsw.tb_fit];
tbND  = [dsw.tb_fitND];
tbout = [out.tb];

%% Some sanity checks
efitzb = (mean(zbin)/mean(zbfit)*Lfit/L)^2*epsilon;
efittb = (mean(tbin)/mean(tbfit)*Tfit/T)^2*epsilon;
disp(['Fitted eps: ',num2str(epsilon_fit),'    zb eps: ',num2str(efitzb),'    tb eps: ',num2str(efittb)]);

% return;
%% Plotting
figure(1); clf;
    plot(...A0',t0_actuals_dim_25','o',...
        [dsw.A0],tbout,'x',...
        [dsw.A0],tbfit,'p',...
        ...[dsw.A0],tbin,'<',...
        'Markersize',8);
    
    
figure(2); clf;
    plot(...A0',z0_actuals_dim_25','o',...
        [dsw.A0],zbout,'x',...
        [dsw.A0],zbfit,'p',...
        ...[dsw.A0],zbin,'<',...
        'Markersize',8);
    
