% Compare theory to exp't
reformat_on = 1; %reformat numerics into matrix; only needs to be done once per trial
wave_type = 'd';
trapping  = 0; trap_lines = 0;
% Plotting parameters
wfac  = 2;
wplot = 8.5*wfac;
hplot = 5.25*wfac;

% data_dir = '/Volumes/APPM-DHL/data/conduit_eqtn/_tmax_200_zmax_600_Nz_6000_order_4_init_condns_soli_tunneling_Amax_1.5_asoli_3_hstretch_2_wave_type_r_bndry_condns_time_dependent/';
% data_dir = '/Volumes/APPM-DHL/data/conduit_eqtn/_tmax_400_zmax_1200_Nz_12000_order_4_init_condns_soli_tunneling_Amax_1.5_asoli_3_hstretch_2_wave_type_r_bndry_condns_time_dependent/';
% data_dir = '/Volumes/APPM-DHL/data/conduit_eqtn/_tmax_400_zmax_1200_Nz_12000_order_4_init_condns_soli_tunneling_Amax_1.5_asoli_3_hstretch_2_wave_type_r_bndry_condns_time_dependent/';
% data_dirRW = '/Volumes/APPM-DHL/data/conduit_eqtn/_tmax_400_zmax_1200_Nz_12000_order_4_init_condns_RW_Amax_1.5_hstretch_1_wave_type_r_bndry_condns_time_dependent/';
% data_dir = '/Volumes/APPM-DHL/data/conduit_eqtn/_tmax_700_zmax_2500_Nz_10000_order_4_init_condns_soli_tunneling_Amax_1.5_asoli_3_hstretch_1_wave_type_r_bndry_condns_time_dependent/';
data_dir = '/Volumes/APPM-DHL/data/conduit_eqtn/_tmax_300_zmax_1500_Nz_6000_order_4_init_condns_soli_tunneling_Amax_2_asoli_2_hstretch_5_wave_type_d_trapping_bndry_condns_time_dependent/';

if trapping
    data_dir = '/Volumes/APPM-DHL/data/conduit_eqtn/_tmax_500_zmax_1500_Nz_6000_order_4_init_condns_soli_tunneling_Amax_2_asoli_1.5_hstretch_5_wave_type_d_trapping_bndry_condns_time_dependent/';
end

load([data_dir,'parameters.mat'],'zjump','z0','Am','hstretch','t','zmax','dz','Nz','asoli');
hstretches = hstretch;
zminus = -(zjump - z0);
uplus = Am;

% Find maximum time in numerics files
tind = length(t)-1;
for ii=1:length(t)-1
  [fid,foo] = fopen(strcat(data_dir,num2str(ii,'%05d'),'.mat'),'r');
  if fid == -1 % File does not exist
    tind = ii-1;
    disp(['Maximum time = ',num2str(t(tind))]);
    break;
  end
  fclose(fid);
end
tm = tind;

% t- & z-axis
t = t(1:tm);
zplot  = dz*[1:Nz]';

if reformat_on
    A_full = zeros(tind,length(zplot));
    for ti = 1:tm
        load(strcat(data_dir,num2str(ti,'%05d'),'.mat'),'A','tnow','inc');
        A_full(ti,:) = A;
    end
    save([data_dir,'matrix.mat'],'A_full','t','zplot');
else
    load([data_dir,'matrix.mat'],'A_full','t','zplot');
end

    % Use phase shift theory to find trajectory
    csoli  = @(a,m) m./a.^2 .* ( (a+m).^2 .* (2*log(1+a./m)-1) + m.^2);
    if ~trapping
        [krat, newc] = ST_phase_shifts(Am, 1, csoli(asoli,Am));
            ps     = zminus*(1-krat);
            zs1     = @(t) (csoli(asoli,Am)*t + zminus) ;
            zs2     = @(t) (newc*t            + zminus - ps);
    else
        [krat, newc] = ST_phase_shifts(Am, 1, csoli(asoli,1));
            ps     = zminus*(1-krat);
            zs1     = @(t) (csoli(asoli,1)*t + zminus) ;
            zs2     = @(t) (newc*t            + zminus - ps);
    end


[foo,zm] = min(abs(max(zplot)-(zs1(t)-zminus+z0)));
zm = zm-5;
f1 = figure(1); clf;
    contourf(zplot,t(1:end-5),A_full(1:end-5,:),100,'edgecolor','none');
        cmap = load('CoolWarmFloat257.csv');
        colormap(cmap); 
        xlabel('z'); ylabel('t'); colorbar;
    hold on;
    if ~trapping | ~ trap_lines
        plot(zs1(t(1:zm))-zminus+z0,t(1:zm),'k--','LineWidth',0.25);
        if ~ trapping
             plot(zs2(t(1:zm))-zminus+z0,t(1:zm),'k--','LineWidth',0.25);
        end
    end
    hold off;
    if numel(hstretches)>1
        input('r');
    end
    set(gca,'FontSize',9)
    if trapping
        savePlot(f1,'soli_DSW_contour_trapping',wplot,hplot,'-dpdf','');
        delete(f1); clear('f1');
        save('soli_DSW_trapping.mat');
    else
        savePlot(f1,'soli_DSW_contour_tunneling',wplot,hplot,'-dpdf','');
        delete(f1); clear('f1');
        save('soli_DSW_tunneling.mat');
    end