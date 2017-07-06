load('box_num_results.mat','ps','ampls','ws','asolis');
for ii=10%1:27
        vstretch = 5; ampl = ampls(ii); w = ws(ii); x0 = 500;
        asoli = asolis(ii);
            if asoli == 7 
                tmax_chosen = 500;
            elseif (ampl==1 && w==200 && asoli==5)
                tmax_chosen = 600;
            else
                tmax_chosen = 800;
            end
            if ampl==1 && w==200 && asoli==3
                tmax = 1200;
            elseif w==200
                tmax = 800;
            else
                tmax = 400;
            end
                ic_type = ['ampl_',num2str(ampl),...
                           '_w_',num2str(w),...
                           '_asoli_',num2str(asoli)];    
            data_dir = ['/Users/dhl/Documents/MATLAB/KdV/data/tmax_',num2str(tmax),...
                        '_zmax_1000_Nz_2048_order_2_init_condns_',...
                        ic_type,'_bndry_condns_periodic/'];
            load([data_dir,'ps.mat'],'prebox','posbox','solilcs','solipks','solitms','dphi','boxpk');
            
            disp('');
            
            
end
