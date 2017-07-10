%soli_tunneling_process_solis
% soli_tunneling_find_solis
aDSW = [5.25 6.25] - 1.75;

aRW  = [7.25 8.53];
m_minus = 1;
m_plus = 1.75;

debug_on = 0;
ADSW = m_plus;
tcutoffdsw = 24;
tcutoffrw  = 10;
solicutoff  = 10;
inds = 1:length(aDSW);
runstyles = 1:2;
run_types = [{'DSW'} {'RW'}];
leadsoliamp = zeros(2,length(inds));
backsoliamp = zeros(2,length(inds));
leadsolispd = zeros(2,length(inds));
backsolispd = zeros(2,length(inds));
for ind = inds
    for rind = runstyles
        run_type = char(run_types(rind));
        data_dir = ['H:\MATLABdata\data\conduit_eqtn\_tmax_125_zmax_800_',...
                    'Nz_80000_order_4_init_condns_soli_tunneling_',...
                    run_type,'_ADSW_',num2str(ADSW),...
                    '_aDSW_',num2str(aDSW(ind)),...
                    '_aRW_' ,num2str(aRW(ind)),...
                    '_bndry_condns_time_dependent/'];
%         load([data_dir, 'parameters.mat']);    
        load([data_dir,'peaks.mat'],'pks','lcs','t','z');
        if findstr(run_type,'DSW');
            pks = pks(1:end-tcutoffdsw);
            lcs = lcs(1:end-tcutoffdsw);
            t   =   t(1:end-tcutoffdsw);
        else
            pks = pks(1:end-tcutoffrw);
            lcs = lcs(1:end-tcutoffrw);
            t   =   t(1:end-tcutoffrw);
        end
        leadsoliamp(rind,ind) = mean(pks(1:solicutoff));
        backsoliamp(rind,ind) = mean(pks(end-solicutoff:end));
        plead = polyfit(t(1:solicutoff),lcs(1:solicutoff),1);
        pback = polyfit(t(end-solicutoff:end),lcs(end-solicutoff:end),1);
        leadsolispd(rind,ind) = plead(1);
        backsolispd(rind,ind) = pback(1);
                % Figure for debugging
        if debug_on
            figure(1); clf;
                subplot(2,1,1);
                    plot(t,lcs,'k.',...
                         t(1:solicutoff),lcs(1:solicutoff),'b.',...
                         t(end-solicutoff:end),lcs(end-solicutoff:end),'r.');
                subplot(2,1,2);
                    plot(t,pks,'k.',...
                         t(1:solicutoff),pks(1:solicutoff),'b.',...
                         t(end-solicutoff:end),pks(end-solicutoff:end),'r.');
                drawnow; input('r');
        end
    end
end
% Plots
fontsize = 12;
% Plot amplitude differences
figure(3); clf;
    plot(aRW+1,aDSW+ADSW,'kx',...
         backsoliamp(1,:),leadsoliamp(1,:),'bx',...  %aDSW+2,aRW+1,'kx',...
         leadsoliamp(2,:),backsoliamp(2,:),'rx',...
         [backsoliamp(1,1)  aRW(1)+1 leadsoliamp(2,1)],...
         [leadsoliamp(1,1) aDSW(1)+ADSW backsoliamp(2,1)],'k-');
     axis equal
     legend('Expected','DSW (out vs in)','RW (in vs out)','Reference Line','Location','SouthEast');
     xlabel('Amplitude on m=1');
     ylabel(['Amplitude on m=',num2str(ADSW)]);
     title('Soliton Amplitudes (pre vs post)');
     set(gca,'FontSize',fontsize)
     print('soliamps','-dpng');
% Plot speeds
csoli = @(a,m) m./a.^2 .* ( (a+m).^2 .* (2*log(1+a./m)-1) + m.^2);
figure(4); clf;
    plot(csoli(aRW,1),csoli(aDSW,ADSW),'kx',...
         backsolispd(1,:),leadsolispd(1,:),'bx',... %DSW
         leadsolispd(2,:),backsolispd(2,:),'rx',...   %RW
         [backsolispd(1,1)  csoli(aRW(1),1) leadsolispd(2,1)],...
         [leadsolispd(1,1)  csoli(aDSW(1),ADSW) backsolispd(2,1)],'k-');    
     axis equal
     legend('Expected','DSW (out vs in)','RW (in vs out)','Reference Line','Location','SouthEast');
     xlabel('Speed on m=1');
     ylabel(['Speed on m=',num2str(ADSW)]);
     title('Soliton Speeds (Pre vs post)');
	 set(gca,'FontSize',fontsize)
     print('solispeeds','-dpng');
% Calculate transmission coeff's
q = @(c,m) c .* (c + 2*m) ./ m;
lead_qs = q(leadsolispd,[ADSW*ones(1,inds(end));   ones(1,inds(end))]);
back_qs = q(backsolispd,[  ones(1,inds(end)); ADSW*ones(1,inds(end))]);
% Plot transmission coeff's
figure(5); clf;
    plot(q(csoli(aRW,1),1), q(csoli(aDSW,ADSW),ADSW),'kx',...
         back_qs(1,:),lead_qs(1,:),'bx',...          %DSW
         lead_qs(2,:),back_qs(2,:),'rx',...          % RW
         [back_qs(1,1)  q(csoli(aRW(1),1),1) lead_qs(2,1)],...
         [lead_qs(1,1)  q(csoli(aDSW(1),ADSW),ADSW) back_qs(2,1)],'k-');    
     axis equal
    legend('Expected','DSW (out vs in)','RW (in vs out)','Reference Line','Location','SouthEast');
    xlabel('q on m=1');
    ylabel(['q on m=2',num2str(ADSW)]);
    title('Transmission Coeffs (pre vs post)');
    set(gca,'FontSize',fontsize)
    print('soliqs','-dpng');
save('soli_tunneling_numerics_large.mat','run_types','m_minus','m_plus',...
                                   'aDSW','aRW','q','csoli',...
                                   'leadsoliamp','leadsolispd','lead_qs',...
                                   'backsoliamp','backsolispd','back_qs');