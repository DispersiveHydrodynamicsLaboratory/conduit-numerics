data_dir = '/Volumes/Data Storage/data/conduit_eqtn/_tmax_100_zmax_1000_Nz_4000_Amax_2_z0_500_w_100_v_10_order_4_init_condns__bndry_condns_periodic/';

load([data_dir,'00000.mat'],'A_init');
load([data_dir,'parameters.mat'],'z','t');
load([data_dir,'matrix.mat'],'Amat');

figure(1); clf;
% Plot IC
filename = 'initial_box';
fontsize = 36;
fontname = 'helvetica';
msize    = fontsize/2;
lwidth   = 3;
w1       = 30*1.618;%(8.25-0.6*2)/4*2.54; % width of plot in cm
h1       = 20; %1.25*2.54;%w1/1.618; % height of plot in cm
units    = 'centimeters';
outerpos = 0;

zmin = 400;
zmax = 700;
Amin = 0;
Amax = 5;

znew = z - zmin;

zminnew = zmin - zmin;
zmaxnew = zmax - zmin;
 

fh2 = figure(1); clf;
    p1 =     plot(znew,A_init,'b-',...
                  'LineWidth',lwidth);
         axis([zminnew zmaxnew Amin Amax]);
            set(fh2,'Units',units,'PaperUnits',units,'Color','White')
            set(gca,'fontsize',fontsize,'fontname',fontname,'LineWidth',lwidth);
            pos = get(fh2,'position');
            set(fh2,'Position',[pos(1),pos(2),w1,h1],...
                   'outerposition',[pos(1)-outerpos,pos(2)-outerpos,w1+2*outerpos,h1+2*outerpos],...
                   'PaperPosition',[0,0,w1,h1],...
                   'PaperSize',[w1,h1])
            xlabel('Vertical Length (ND)','HorizontalAlignment','center',...
                            'Interpreter','latex','Fontsize',fontsize);
            yname = ylabel('Area (ND)',...
                            'Interpreter','latex','Fontsize',fontsize);
%             legend('Numerics','Experiment','Location','Best');
savePlot(fh2,filename,w1,h1);

% Plot a specific point in space
tchosen  = 100;
filename = 'box_numerics_evolved';
fontsize = 36;
fontname = 'helvetica';
msize    = fontsize/2;
lwidth   = 3;
w1       = 30*1.618;%(8.25-0.6*2)/4*2.54; % width of plot in cm
h1       = 20; %1.25*2.54;%w1/1.618; % height of plot in cm
units    = 'centimeters';
outerpos = 0;

zmin = 600;
zmax = 1000;
Amin = 0;
Amax = 8;

[~,tind] = min(abs(t-tchosen));

znew = z - 400;

zminnew = zmin - 400;
zmaxnew = zmax - 400;
 

fh2 = figure(2); clf;
    p1 =     plot(znew,Amat(tind,:),'b-',...
                  'LineWidth',lwidth);
         axis([zminnew zmaxnew Amin Amax]);
            set(fh2,'Units',units,'PaperUnits',units,'Color','White')
            set(gca,'fontsize',fontsize,'fontname',fontname,'LineWidth',lwidth);
            pos = get(fh2,'position');
            set(fh2,'Position',[pos(1),pos(2),w1,h1],...
                   'outerposition',[pos(1)-outerpos,pos(2)-outerpos,w1+2*outerpos,h1+2*outerpos],...
                   'PaperPosition',[0,0,w1,h1],...
                   'PaperSize',[w1,h1])
            xlabel('Vertical Length (ND)','HorizontalAlignment','center',...
                            'Interpreter','latex','Fontsize',fontsize);
            yname = ylabel('Area (ND)',...
                            'Interpreter','latex','Fontsize',fontsize);
%             legend('Numerics','Experiment','Location','Best');
savePlot(fh2,filename,w1,h1);

