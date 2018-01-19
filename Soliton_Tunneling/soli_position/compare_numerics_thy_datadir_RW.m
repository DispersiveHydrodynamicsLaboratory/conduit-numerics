% Compare theory to exp't
reformat_on = 1; %reformat numerics into matrix; only needs to be done once per trial
wave_type = 'r';
trapping  = 1;
% Plotting parameters
wfac  = 2;
wplot = 8.5*wfac;
hplot = 5.25*wfac;

if trapping
    data_dir = '/Volumes/APPM-DHL/data/conduit_eqtn/_tmax_250_zmax_1000_Nz_4000_order_4_init_condns_soli_tunneling_Amax_2_asoli_2_hstretch_5_wave_type_r_bndry_condns_time_dependent/';
else
    data_dir = '/Volumes/APPM-DHL/data/conduit_eqtn/_tmax_400_zmax_2000_Nz_8000_order_4_init_condns_soli_tunneling_Amax_2_asoli_7_hstretch_5_wave_type_r_trapping_bndry_condns_time_dependent/';
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

    % Calculate theoretical soliton trajectory
    [zs] = soliton_position_fun(1,asoli,zminus,Am);

[foo,zm] = min(abs(max(zplot)-(zs(t)-zminus+z0)));
f1 = figure(1); clf;
    contourf(zplot,t,A_full,100,'edgecolor','none');
        cmap = load('CoolWarmFloat257.csv');
        colormap(cmap); 
        xlabel('z'); ylabel('t'); colorbar;
    hold on;
    plot(zs(t(1:zm))-zminus+z0,t(1:zm),'k--','LineWidth',0.25);
    hold off;
    if numel(hstretches)>1
        input('r');
    end
    set(gca,'FontSize',9)
    if trapping
        savePlot(f1,'soli_RW_contour_trapping',wplot,hplot,'-dpdf','');
        delete(f1);clear('f1');
                save('soli_RW_trapping.mat');
    else
        savePlot(f1,'soli_RW_contour_tunneling',wplot,hplot,'-dpdf','');
        delete(f1);clear('f1');
                save('soli_RW_tunneling.mat');
    end
return;
% Compare RW to theoretical and numerical options
%     % Pick random time step in the middle
%         tchosen = tm-110;%floor(tm/2);
%         Achosen = A_full(tchosen,:);
%         figure(3); clf;
%             plot(zplot,Achosen);
    % Theoretical RW (zero dispersion)
    RW = @(z,t) ones(size(z))       .* (z<=(2*t))               + ...
                z/(2*t)             .* (z>2*t) .* (z<2*uplus*t) + ...
                uplus*ones(size(z)) .* (z>=2*uplus*t);
	RW_IC = @(z) ones(size(z))       .* (z<=0)               + ...
                 uplus*ones(size(z)) .* (z>=0);
% 	% Numerical RW (full conduit EQ)
%     RW_full = load([data_dirRW,'matrix.mat'],'A_full');
%         figure(3); hold on;
%             plot(zplot,RW(zplot-zjump,t(tchosen)),'LineWidth',2);
RW_err_thy = zeros(size(A_full));
% RW_err_num = A_full - RW_full.A_full;
for ti = 2:tm
	Achosen = A_full(ti,:);
    RW_err_thy(ti,:)  = Achosen' - RW(zplot-zjump,t(ti));
%     RW_thy_num(ti,:)  = [RW_full.A_full(ti,:)]'- RW(zplot-zjump,t(ti));
end
cutoff = 0.1;
RW_err_thy(RW_err_thy>cutoff) = cutoff; 
RW_err_thy = abs(RW_err_thy);
f4 = figure(4); clf;
    contourf(zplot(1:25:end),t,RW_err_thy(:,1:25:end),100,'edgecolor','none');
        cmap = load('CoolWarmFloat257.csv');
        colormap(cmap); 
        xlabel('z'); ylabel('t'); colorbar;
        set(gca,'Fontsize',9)
        savePlot(f4,'contour_err_trunc',wplot,hplot,'-dpdf','')
        
% figure(6); clf;
%     plot(zplot,A_full(1,:),zplot,RW_IC(zplot-zjump));
%     axis([zplot(1) zplot(end) 0.9 2]); drawnow; pause(0.1);
%     for ti = 2:tm
%         figure(6); clf;
%         plot(zplot,A_full(ti,:),zplot,RW(zplot-zjump,t(ti))); axis([zplot(1) zplot(end) 0.9 2]);drawnow; pause(0.1);
%     end
        
        
% figure(5); clf;
%     contourf(zplot(1:25:end),t,RW_err_num(:,1:25:end),100,'edgecolor','none');
%         cmap = load('CoolWarmFloat257.csv');
%         colormap(cmap); 
%         xlabel('z'); ylabel('t'); colorbar;
% figure(6); clf;
%     contourf(zplot(1:25:end),t,RW_thy_num(:,1:25:end),100,'edgecolor','none');
%         cmap = load('CoolWarmFloat257.csv');
%         colormap(cmap); 
%         xlabel('z'); ylabel('t'); colorbar;  
            