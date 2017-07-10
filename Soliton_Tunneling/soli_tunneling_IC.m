function[fRW, fDSW] = soli_tunneling_IC( A, aDSW, aRW, zmax)
debug_on = 0;
% Generate jump conditions
    vscale   = ((A-1)/2); 
    upshift  = 1+vscale;
    hstretch = 5;
    zjumpDSW      = 700*1/2;
    zjumpRW       = 700*1/4;
    DSW      = @(z) upshift - vscale*tanh(1/hstretch*(z-zjumpDSW));
    RW       = @(z) upshift + vscale*tanh(1/hstretch*(z-zjumpRW));
    
% Generate Solitons
z0DSW = 2/3*zjumpDSW;
z0RW  = 3/4*zjumpRW;
    [ zDSW,phiDSW ] = conduit_soliton_newton_cg( aDSW/A, zmax, debug_on );
        soliDSW = @(z) interp1(zDSW*sqrt(A)+z0DSW,phiDSW*A,z,'spline',0);
    [  zRW,phiRW  ] = conduit_soliton_newton_cg( aRW,  zmax, debug_on );
        soliRW  = @(z) interp1(zRW +z0RW ,phiRW ,z,'spline',0);
     
% Putting it all together
    fDSW = @(z) DSW(z) + soliDSW(z);
    fRW  = @(z)  RW(z) + soliRW(z);

if debug_on
% Figure for debugging
zplot = 0:0.1:zmax;
figure(3); clf;
    subplot(2,1,1);
        plot(zplot,fDSW(zplot));
	subplot(2,1,2);
        plot(zplot, fRW(zplot));
end
