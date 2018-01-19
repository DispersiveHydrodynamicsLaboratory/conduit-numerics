save_on = 1;
jumps = 8;%[1.5:0.5:5];

for jump = jumps
    data_dir = ['D:\Numerics\conduit_numerics_current\data\conduit_eqtn_tmax_200_zmax_800_Nz_80000_order_4_init_condns_DSW_like_jump_4_bndry_condns_time_dependent\'];
        
	load([data_dir,'parameters.mat'],'t','Nz','dz','zmax');
    z = dz*[1:Nz];
	% Find maximum t index
    tmax = length(t);
    for ii=1:length(t)-1
      [fid,foo] = fopen(strcat(data_dir,num2str(ii,'%05d'),'.mat'),'r');
      if fid == -1 % File does not exist
        tmax = ii-1;
        disp(['Maximum time = ',num2str(t(tmax))]);
        break;
      end
      fclose(fid);
    end
    if save_on
        for ti = 1:tmax
            load(strcat(data_dir,num2str(ti,'%05d'),'.mat'),'A','tnow','inc');
            figure(1); clf;
            plot(z,A);
            drawnow; pause(0.1);
        end
    end
end