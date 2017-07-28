% Looks at soliton tunneling theory in the context of a soliton-RW
% interaction
% Compares numerical soliton trajectory to following 
% the soliton characteristic
plot_trajectory_on = 1; 
% Compares numerical RW to theoretical RW
plot_RWerr_on = 0;


% Plotting parameters
wplot = 8.5*1; % width of plots, cm
hplot = 5.25*1;% height of plots, cm
fs    = 9;     % fontsize, pts

% Load numerics
load('tmax_700_zmax_2500_Nz_10000_order_4_init_condns_soli_tunneling_Amax_1.5_asoli_3_hstretch_1_wave_type_r.mat','zplot','t','A_full','f');
tm = length(t);
hstretch = f.hstretch; % factor by which tanh profile is stretched in z
z0 = f.z0;             % absolute location of soliton
zjump = f.zjump;       % absolute location of jump
zminus = -(f.zjump - f.z0); % RELATIVE location of soliton to the jump
uplus = f.Am;          % Jump factor
asoli = f.asoli;       % Soliton amplitude

% Calculate theoretical soliton position vis characteristic 
% dz_s(t)/d(t) = c_s( q_0,z_s(t)/(2*t) ), q_0 tunneling reciprocity factor
[zs] = soliton_position_fun(1,asoli,zminus,uplus);

[foo,zm] = min(abs(max(zplot)-(zs(t)-zminus+z0)));

% Plot results
if plot_trajectory_on
    f1 = figure(1); clf;
        contourf(zplot(1:25:end),t,A_full(:,1:25:end),100,'edgecolor','none');
            cmap = load('CoolWarmFloat257.csv');
            colormap(cmap); 
            xlabel('z'); ylabel('t'); colorbar;
        hold on;
        plot(zs(t(1:zm))-zminus+z0,t(1:zm),'k--','LineWidth',0.25);
        hold off;
        set(gca,'FontSize',fs)
    savePlot(f1,'soli_RW_contour',wplot,hplot,'-dpdf','')
end

% Compare numerical RW to theoretical RW
if plot_RWerr_on
        % Theoretical RW (zero dispersion)
        RW = @(z,t) ones(size(z))       .* (z<=(2*t))               + ...
                    z./(2*t)             .* (z>2*t) .* (z<2*uplus*t) + ...
                    uplus*ones(size(z)) .* (z>=2*uplus*t);
        RW_IC = @(z) ones(size(z))       .* (z<=0)               + ...
                     uplus*ones(size(z)) .* (z>=0);
    RW_err_thy = zeros(size(A_full));
        RW_err_thy(1,:)     = RW_IC(zplot);
        RW_err_thy(2:tm,:)  = A_full(2:tm,:) - [RW(zplot-zjump,t(2:tm))]';
	% Put a ceiling on error from soliton
    cutoff = 0.1;
    RW_err_thy(RW_err_thy>cutoff) = cutoff; 
    
    RW_err_thy = abs(RW_err_thy);
    f4 = figure(4); clf;
        contourf(zplot(1:25:end),t,RW_err_thy(:,1:25:end),100,'edgecolor','none');
            cmap = load('CoolWarmFloat257.csv');
            colormap(cmap); 
            xlabel('z'); ylabel('t'); colorbar;
            set(gca,'Fontsize',fs)
            savePlot(f4,'contour_err_trunc',wplot,hplot,'-dpdf','')
end         