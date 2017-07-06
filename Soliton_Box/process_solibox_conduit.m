save_on       = 1;
plot_debug_on = 0;
plot_solis_on = 0;
process_on    = 0;
choose_tmax   = 1;
tmax_chosen = 800;

%% Preallocate arrays
trialnum = 3^3;
ps     = zeros(1,trialnum);
ampls  = zeros(1,trialnum);
ws     = zeros(1,trialnum);
asolis = zeros(1,trialnum);
boxspd = zeros(1,trialnum);
boxpks = struct;
% ctr = 1;
 amplsa = [0.25 0.5 1];
 wsa = [50 100 200];
asolisa = [3 5 7];
[AMPLS,WS,ASOLIS] = meshgrid(amplsa,wsa,asolisa);

for runi = 1:27
    ampl  = AMPLS(runi);
    w     = WS(runi);
    asoli = ASOLIS(runi);
            if ismember(runi,[1,19,20,22,25])
                tmax_chosen = 200;
            elseif ismember(runi,[10 11 13 14 16 23 26])
                tmax_chosen = 250;
            elseif ismember(runi,[8 27])
                tmax_chosen = 450;
            else
                tmax_chosen = 900;
            end
%             if ampl==1 && w==200 && asoli==3
%                 tmax = 1200;
%             elseif w==200
%                 tmax = 800;
            if ismember(runi,[3 6 7 8 18 27])
                tmax = 800;
            elseif runi==9
                tmax = 900;
            else
                tmax = 400;
            end
%             vstretch = 5;
%             f1 = make_smooth_box(vstretch,ampl,w,500);
%             % Parameters for conduit soliton
%                     Nexp = 12;
%                     L    = 200;
%                 if ampl==1 && ismember(w,[100 200]) && asoli==3
%                     zsoli = 500-w;
%                     [ zphi,phi ] = conduit_soliton_newton_cg( asoli, Nexp, L, 0 );
%                 else
%                     [ zphi,phi ] = conduit_soliton_newton_cg( asoli, Nexp, L, 0 );
%                     zsoli = 500-2*w;
%                 end
%         
%         f2 = @(zq) interp1(zphi+zsoli,phi,zq,'spline',0);
%             f = @(z) f1(z) + f2(z);
                ic_type = ['ampl_',num2str(ampl),...
                           '_w_',num2str(w),...
                           '_asoli_',num2str(asoli)];    
            data_dir = ['/Users/dhl/Documents/MATLAB/conduit/data/_tmax_',num2str(tmax),...
                        '_zmax_1000_Nz_10000_order_4_init_condns_',...
                        ic_type,'_bndry_condns_periodic/'];
            if runi==9
                data_dir = '/Users/dhl/Documents/MATLAB/conduit/data/_tmax_900_zmax_1400_Nz_14000_order_4_init_condns_ampl_1_w_200_asoli_3_bndry_condns_periodic/';
            end
                       load([data_dir,'parameters.mat'],'t','zmax','dz');
            disp(ic_type);
            if process_on
                disp('Calculating maximum time increment in saved data files...');
                for ii=1:length(t)+1
                    [fid,foo] = fopen(strcat(data_dir,num2str(ii,'%05d.mat')),'r');
                    if fid == -1 % File does not exist
                        tm = ii-1;
                        disp(['Maximum time = ',num2str(t(tm))]);
                        break;
                    end
                    tm = ii-1;
                    fclose(fid);
                end
                    disp(['Maximum time = ',num2str(t(tm))]);
                    if choose_tmax
                        tm1 = find(t>=tmax_chosen,1,'first');
                        if isempty(tm1)
                            tm1 = tmax;
                        end
                        tm = min(tm,tm1);
                    end
                % Get rid of larger t values
                t = t(1:tm);

                % Find soliton for all time
                z       = dz:dz:zmax;
                solipks = zeros(size(t));
                solilcs = zeros(size(t));
                boxpk   = struct;
                ctr0 = 1;
                for ti = 1:tm
                    load(strcat(data_dir,num2str(ti,'%05d.mat')),'A','tnow');
                    [pk,lc] = findpeaks(A,z,'SortStr','descend','MinPeakProminence',0.05); 
                    solipks(ti) = pk(1);
                    solilcs(ti) = lc(1);
                    if ti>tm-10
                        boxpk.(['t',num2str(ctr0)]) = pk(2:end);
                        ctr0 = ctr0+1;
                    end
                    
                    if plot_debug_on
                        % Figure for debugging
                        figure(2); clf;
                        plot(z,A,'-',lc,pk,'*');
                        title(num2str(tnow));
                        drawnow; %pause(0.05);
                    end
                end
            [ solilcs ] = numunwrap( solilcs, zmax ); % unwraps locations (periodic BCs)
            badinds = find(diff(solilcs)<0) + 1; % checks for bad indices
            solipks(badinds) = [];
            solilcs(badinds) = [];
            solitms = t;
            solitms(badinds) = [];
            %% Compare pre and post box
            thresh  = mean(solipks(1:10))*0.975;
            boxtime = find(solipks<thresh);
            prebox  = 1:min(boxtime);
            if tmax-max(boxtime)<10 % Soliton starts going through the box again
                d = diff(boxtime);
                inds = find(d>10);
                boxtime = boxtime(1:min(inds));
            end
            posbox  = max(boxtime):tm-numel(badinds);
            
            %% Find speed inside box
            trueboxtm = max(prebox):min(posbox);
            pbox = polyfit(solitms(trueboxtm),solilcs(trueboxtm),1);

            %% Calculate char. line for soliton
            prexi = polyfit(solitms(prebox),solilcs(prebox),1);
            posxi = polyfit(solitms(posbox),solilcs(posbox),1);
                dphi = posxi(2) - prexi(2);

            if plot_solis_on
            %% Plot results
                figure(1); clf;
                    h(1) = subplot(2,1,1);
                                plot(solitms,solilcs,'k.',...
                                     t,polyval(prexi,t),'b-',...
                                     t,polyval(posxi,t),'r-');
                                legend('soliton','pre','post');
                                title(['Trial ',num2str(runi),' Phase Shift: ',num2str(dphi)]);
                    h(2) = subplot(2,1,2);
                                plot(solitms,solipks,'k.',...
                                     solitms(prebox),solipks(prebox),'b.',...
                                     solitms(posbox),solipks(posbox),'r.',...
                                     solitms(trueboxtm),solipks(trueboxtm),'g.',...
                                     'MarkerSize',8);
                    linkaxes(h,'x');
                    drawnow; pause(0.5);
            end

                save([data_dir,'ps.mat'],'prebox','posbox','solilcs','solipks','solitms','dphi','boxpk','pbox');
            else
                load([data_dir,'ps.mat'],'prebox','posbox','solilcs','solipks','solitms','dphi','boxpk','pbox');
            end
            ps(runi)     = dphi;
            ampls(runi)  = ampl;
            ws(runi)     = w;
            asolis(runi) = asoli;
            boxpks.(['run',num2str(runi)]) = boxpk;
            boxspd(runi) = pbox(1);
%             ctr = ctr+1;
end

    
% figure(1); clf;
%     boxparam = ampls.^(1/2).*ws;
%     plot(boxparam,ps,'*');
    if save_on
        save('box_num_results.mat','ps','ampls','ws','asolis','boxpks','boxspd');%,'-append');
    end