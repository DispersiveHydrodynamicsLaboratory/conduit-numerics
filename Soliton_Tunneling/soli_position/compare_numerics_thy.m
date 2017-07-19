% Compare theory to exp't
%% IC Parameters
Am     = 2;
asoli = 3;
hstretch = 2;
zjump = zmax/2;
z0    = zjump - 50; % true soliton position
wave_type = 'r'; % r for RW, d for DSW
zminus = -(zjump-z0); % distance between soliton and RW initially
                                 % from soli_tunneling_IC
reformat_on = 1; %reformat numerics into matrix; only needs to be done once per trial

data_dir = '/Volumes/APPM-DHL/data/conduit_eqtn/_tmax_125_zmax_700_Nz_7000_order_4_init_condns_soli_tunneling_Amax_2_asoli_3_hstretch_2_wave_type_r_bndry_condns_time_dependent/';
load([data_dir,'parameters.mat'],'t','dz','Nz');

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

[zs] = soliton_position_fun(1,asoli,zminus,Am);

figure(1); clf;
    contourf(zplot(1:25:end),t,A_full(:,1:25:end),100,'edgecolor','none');
        cmap = load('CoolWarmFloat257.csv');
        colormap(cmap); 
        xlabel('z'); ylabel('t'); colorbar;
    hold on;
    plot(zs(t)-zminus+z0,t,'k','LineWidth',2);
    hold off;

    
