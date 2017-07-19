function[f] = soli_tunneling_IC_varied( A, asoli, zmax, hstretch, zjump, z0, wave_type)
debug_on = 0;
% Generate jump conditions
    vscale   = ((A-1)/2); 
    upshift  = 1+vscale;
    if strcmp(wave_type,'d')
        wave      = @(z) upshift - vscale*tanh(1/hstretch*(z-zjump));
        [ zDSW,phiDSW ] = conduit_soliton_newton_cg( asoli/A, zmax, debug_on );
        soli = @(z) interp1(zDSW*sqrt(A)+z0,phiDSW*A,z,'spline',0);
    elseif strcmp(wave_type,'r')
        wave       = @(z) upshift + vscale*tanh(1/hstretch*(z-zjump));
        [  zRW,phiRW  ] = conduit_soliton_newton_cg( asoli,  zmax, debug_on );
        soli  = @(z) interp1(zRW +z0 ,phiRW ,z,'spline',0);
    else
        disp('Please use d or r for DSW or RW, respectively');
        return;
    end
     
% Putting it all together
    f  = @(z)  wave(z) + soli(z);

if debug_on
% Figure for debugging
zplot = 0:0.1:zmax;
figure(3); clf;
        plot(zplot, f(zplot));
end
