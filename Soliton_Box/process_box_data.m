clear data_dirs
ii = 1;
for ampl = 0.5:0.5:3
    for boxmax = [50 100:100:500]
        vstretch = 30;
        x0 = 500;
        ic_type = ['smooth_box_vstretch_',num2str(vstretch),...
                   '_ampl_',num2str(ampl),'_boxmax_',num2str(boxmax),...
                   '_x0_',num2str(x0)];
        data_dirs{ii} =  ['H:\MATLABdata\data\conduit_eqtn\_tmax_50_zmax_1000_Nz_100000_order_4_init_condns_',...
                     ic_type,'_bndry_condns_periodic/'];
             ii = ii+1;
    end
end
script_extract_envelope_fun(data_dirs)