% Computes theoretical soliton trajectory through a rarefaction wave
% Requires soliton_position_fun.m to run

uminus = 1;
asoli  = 5;
zminus = 50; % RW starts at 0, soliton starts at -zminus
uplus  = 2;

% Calculate theoretical soliton position via characteristic 
% dz_s(t)/d(t) = c_s( q_0,z_s(t)/(2*t) ), q_0 tunneling reciprocity factor
[zs,tmax] = soliton_position_fun(uminus,asoli,zminus,uplus);

% Plot leading and trailing edges of RW as well as soliton trajectory
figure(1); clf;
    t = 0:tmax;
    plot(t,2*t,'b--',...
         t,2*uplus*t,'b-.',...
         t,zs(t),'r-',...
         'LineWidth',2);
     xlabel('t'); ylabel('z');
     legend('RW-trailing edge','RW-leading edge','Soliton Trajectory',...
            'Location','SouthEast');
     


