% soli_tunneling_find_solis
% load('soli_pairs.mat','aDSW','aRW','m_minus','m_plus');
aDSW = [5.25 6.25] - 1.75;
aRW  = [7.25 8.53];
m_minus = 1;
m_plus = 1.75;

ADSW = m_plus;
debug_on = 1;
for ind = 1:length(aDSW)
    for run_type = [{'DSW'} {'RW'}]
        data_dir = ['H:\MATLABdata\data\conduit_eqtn\_tmax_125_zmax_800_',...
                    'Nz_80000_order_4_init_condns_soli_tunneling_',...
                    char(run_type),'_ADSW_',num2str(ADSW),...
                    '_aDSW_',num2str(aDSW(ind)),...
                    '_aRW_' ,num2str(aRW(ind)),...
                    '_bndry_condns_time_dependent/'];
        load([data_dir, 'parameters.mat']);     
        z = dz*[1:Nz]';
        pks = zeros(size(t));
        lcs = zeros(size(t));
        for tind = 1:length(t)
            load([data_dir,sprintf('%05d.mat',tind)],'A','tnow');
            [pk, loc] = findpeaks(A,z,'SortStr','descend','NPeaks',1);
            pks(tind) = pk;
            lcs(tind) = loc;
        end
        save([data_dir,'peaks.mat'],'pks','lcs','t','z');
        % Figure for debugging
        if debug_on
            figure(1); clf;
                subplot(2,1,1);
                    plot(t,lcs,'.');
                subplot(2,1,2);
                    plot(t,pks,'.');
                drawnow; input('r');
        end
    end     
end