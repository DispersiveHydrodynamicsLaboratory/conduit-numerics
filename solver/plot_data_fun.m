function plot_data_fun(loaddir,varargin)
load([loaddir,'parameters.mat'],'t','Nz','dz','zmax'); %,'f');
z = dz*[1:Nz];
if nargin>1
    nplts = varargin{1};
    if nargin>2
        tmax = varargin{2};
    else
        tmax = Inf;
    end
else
    nplts = 6;
    tmax  = Inf;
end

% Find maximum t index
tind = length(t)-1;
for ii=1:length(t)-1
  [fid,foo] = fopen(strcat(loaddir,num2str(ii,'%05d'),'.mat'),'r');
  if fid == -1 % File does not exist
    tind = ii-1;
    disp(['Maximum time = ',num2str(t(tind))]);
    break;
  end
  fclose(fid);
end

t = t(1:min(tind,tmax));

% Plotting stuff
plot_figs = [1   % Time sequence of n plots
             0   % Boundary deviations
             1   % Initial vs final time
             0]; % Conserved quantity


fontsize = 12;
[foo,zmini] = min(abs(z - 0));
[foo,zmaxi] = min(abs(z - zmax));
zp = zmini:zmaxi;

if plot_figs(1)
    %toutind = 1:6;
    toutind = round(linspace(round((length(t)-1)/(nplts-1)),length(t)-1,nplts-1));
    % Find max and min of solution to plot
    if exist('f','var')
        Amax = max(f(z(zp)));
        Amin = min(f(z(zp)));
    else
        load(strcat(loaddir,num2str(0,'%05d'),'.mat'),'A_init');
        Amax = max(A_init); %max(f(z(zp)));
        Amin = min(A_init); %min(f(z(zp)));
    end
    for ii=1:nplts-1
        load(strcat(loaddir,num2str(toutind(ii),'%05d'),'.mat'),'A','tnow','inc');
        if max(A) > Amax
            Amax = max(A);
        end
        if min(A) < Amin
            Amin = min(A);
        end
    end
    
    figure(4)
    clf()
    % Plot initial condition
    load(strcat(loaddir,num2str(0,'%05d'),'.mat'),'A_init');
    h(1)=subplot(nplts,1,1);
    plot(z(zp),A_init(zp),'b-');
    set(gca,'fontsize',fontsize,'fontname','times');
    ylabel('$A$','interpreter','latex');
    axis([z(zmini),z(zmaxi),Amin-0.2,Amax+0.2]);
    title({['$t = 0$']},'interpreter','latex');
    for ii=1:nplts-1
        load(strcat(loaddir,num2str(toutind(ii),'%05d'),'.mat'),'A','tnow','inc');
        h(ii+1) = subplot(nplts,1,ii+1);
        plot(z(zp),A(zp),'b-');
        set(gca,'fontsize',fontsize,'fontname','times');
        ylabel('$A$','interpreter','latex');
        if ii == nplts-1
            xlabel('$z$','interpreter','latex');
        end
        axis([z(zmini),z(zmaxi),Amin-0.2,Amax+0.2]);
        title({['$t = ',num2str(t(toutind(ii)+1)),'$']},'interpreter','latex');
    end
end
linkaxes(h,'xy')

if plot_figs(2)
    % Extract solution values adjacent to boundaries
    A1 = zeros(1,tind);
    AN = zeros(1,tind);
    tp = zeros(1,tind);
    for ii=1:tind
        load(strcat(loaddir,num2str(ii,'%05d'),'.mat'),'A','tnow','inc');
        A1(ii) = A(1);
        AN(ii) = A(Nz);
        tp(ii) = tnow;
    end
    
    figure(4)
    clf()
    subplot(1,2,1);
    plot(tp,g0(tp),'b-',...
         tp,A1,'r--');
    set(gca,'fontsize',fontsize,'fontname','times');
    legend({['$A(0,t)$'],['$A(z_1,t)$']},'interpreter','latex',...
           'location','northwest');
    xlabel('$t$','interpreter','latex');
    subplot(1,2,2);
    plot(tp,g1(tp),'b-',...
         tp,AN,'r--');
    set(gca,'fontsize',fontsize,'fontname','times');
    legend({['$A(L,t)$'],['$A(z_N,t)$']},'interpreter','latex',...
           'location','northwest');
    xlabel('$t$','interpreter','latex');
end

% Initial and final times
if plot_figs(2)
    load(strcat(loaddir,num2str(tind,'%05d'),'.mat'),'A','tnow','inc');
    
    figure(4)
    clf()
    subplot(1,2,1);
    plot(z,f(z),'b-',...
         z,A,'r--');
    set(gca,'fontsize',fontsize,'fontname','times');
    legend({['$A(z,0)$'],['$A(z,t_{\rm f})$']},'interpreter','latex',...
           'location','northwest');
    xlabel('$z$','interpreter','latex');
    subplot(1,2,2);
    plot(z,f(z)'-A,'b-');
    set(gca,'fontsize',fontsize,'fontname','times');
    xlabel('$z$','interpreter','latex');
    ylabel('$A(z,0)-A(z,t_{\rm f})$','interpreter','latex');
end


