%% Plots Soliton Tunneling theory, numerics, and experiment
close all;
    load('stdata_plot.mat');
    Ttop = 4:Ttop(end);
    Tbot = Tbot;
% ----- plot parameters -----
fs   = 9; %fontsize
w    = 8.5; %cm
h    = 5.25; %cm
n    = 10;
qaxis = [ 3 10.5 2 7.3];
xtick = [4:2:10];
ytick = [2:2:6];
   fontsize = fs;
   msize = 3; %markersize
   % line specs:  DSWs    RWs    Thy     NumDSW NumRW 
   linespecs = [ {'bs'}, {'rv'}, {'k-'}, {'ks'}, {'kv'} ];
   numcolor  = 'none'; % Numerics markerface color
    %% Figure of thy, numerics, and exp't (corrected diams)
    f3 = figure(3); clf;hold on;
    % Plot theory line and numerics (DSW then RW)
        h3all = plot(ARW,ADSW,linespecs{3});%,...      ST1data.pos.ampl(dsw1).*ST1_corr(dsw1),ST1data.pre.ampl(dsw1)              ,'ks',...ST1data.pre.ampl(rw1)                 ,ST1data.pos.ampl(rw1).*ST1_corr(rw1),'ko',...SBdata.pos.ampl(dsw2).*ST2_corr      ,SBdata.pre.ampl(dsw2)              ,'ks',...SBdata.pre.ampl(rw2)                 ,SBdata.pos.ampl(rw2).*ST2_corr     ,'ko',...
%                   );
            set(h3all(1),'LineWidth',2);
    %   errorbarxy_for_publication(x, y, lerrx, uerrx, lerry, uerry) plots the data with errorbars on both 
    %   x and y axes with error bars [x-lerrx, x+uerrx] and [y-lerry, y+uerry]. If there is
    %   no error on one axis, set corresponding lower and upper bounds to [].
    %   For some reason, need hold on in between each call to this fcn
    % Plot Experimental data and errorbars
                                     hold on;
	% DSWs
    h32d = errorbarxy_for_publication(SBdata.top.ampl(Ttop)*SB_top_corr,SBdata.mid.ampl(Ttop),...
                                sqrt( (SBdata.top.amplStd(Ttop)).^2 + SB_corrstd),sqrt( (SBdata.top.amplStd(Ttop)).^2 + SB_corrstd),...
                                SBdata.mid.amplStd(Ttop),(SBdata.mid.amplStd(Ttop)),...
                                      {linespecs{1}, 'k', 'k'});
           h32d.hMain.MarkerFaceColor = linespecs{1}(1);
           h32d.hMain.MarkerSize = msize;
                                      hold on;
	%RWs
    h32r = errorbarxy_for_publication(SBdata.bot.ampl(Tbot),SBdata.mid.ampl(Tbot),...
                                      SBdata.bot.amplStd(Tbot),SBdata.bot.amplStd(Tbot),...
                                      SBdata.mid.amplStd(Tbot),SBdata.mid.amplStd(Tbot),...
                                      {linespecs{2}, 'k', 'k'});
           h32r.hMain.MarkerFaceColor = linespecs{2}(1);
           h32r.hMain.MarkerSize = msize;
                                      hold on;
                                      
h3num = plot(numdata.backsoliamp(1,:),numdata.leadsoliamp(1,:),linespecs{4},...
             numdata.leadsoliamp(2,:),numdata.backsoliamp(2,:),linespecs{5},...[min(ARW) max(ARW)], polyval(polyfit(numdatax,numdatay,1),[min(ARW) max(ARW)]),linespecs{4},...[min(ARW) max(ARW)], polyval(polyfit(exptx,expty,1),[min(ARW) max(ARW)])        ,'k-.',...
             'Markersize',msize);
            set(h3num(1),'MarkerFaceColor',numcolor,'Linewidth',0.75);
            set(h3num(2),'MarkerFaceColor',numcolor,'Linewidth',0.75);

	% General plotting things
            axis(qaxis);
            set(gca,'XTick',xtick,'YTick',ytick,'Fontsize',fontsize,...
                    'TickLabelInterpreter','latex');
            
            xlabel('$a ~ (\mathrm{smaller}~\overline{u})$','interpreter','latex')
            ylabel('$a ~ (\mathrm{larger}~\overline{u})$','interpreter','latex')
            [leg,cobjs] = legend([h3all(1) h32d.hMain h32r.hMain h3num(1) h3num(2)],{'Theory','Expt-DSWs','Expt-RWs',...
                   'Numerics-DSWs','Numerics-RWs'},'interpreter','latex','FontSize',7,'Location','SouthEast');
            cobjs(6).XData = [0.15 0.27];
            leg.Position = [0.57 0.18 0.3 0.45]; 
            leg.Box = 'off';
                set(gca,'FontSize',fontsize)
                savePlot(f3,'thy_num_expt',w,h)
	hold off;