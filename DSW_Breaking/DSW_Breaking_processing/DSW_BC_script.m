z0 = 50;
Aplus = 2; Aminus = 1;
t0 = z0/(2*Aplus);
% deltat = dzinit/100;
t = 0:1/10:2*t0;
z = 0:1/10:z0;
    g0      = @(t) BC_DSW( t, t0 , z, z0, Aplus , Aminus);         % BCs at z=0
    dg0     = @(t) BC_DSW_deriv( t, t0 , z, z0, Aplus , Aminus);   % derivative of BC
    g1      = @(t) ones(size(t)) ;              % BCs at z=zmax
    dg1     = @(t) zeros(size(t));              % derivative of BC
    f       = @(z) ones(size(z));               % ICs at t=0
A = g0(t);
size(A)
figure(4);    
contourf(z,t,g0(t)); colorbar;
xlabel('z'); ylabel('t');
