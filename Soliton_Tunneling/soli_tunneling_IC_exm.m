function[fDSW] = soli_tunneling_IC_exm( A, aDSW, zmax)
debug_on = 0;
% Generate jump conditions
    vscale   = ((A-1)/2); 
    upshift  = 1+vscale;
    hstretch = 5;
    zjumpDSW      = 700*1/2;
    zjumpRW       = 700*1/4;
    DSW      = @(z) upshift - vscale*tanh(1/hstretch*(z-zjumpDSW));
    
% Generate Solitons
z0DSW = 2/3*zjumpDSW;
    [ zDSW,phiDSW ] = conduit_soliton_newton_cg( aDSW/A, zmax, debug_on );
        soliDSW = @(z) interp1(zDSW*sqrt(A)+z0DSW,phiDSW*A,z,'spline',0);
     
% Putting it all together
    fDSW = @(z) DSW(z) + soliDSW(z);

if debug_on
% Figure for debugging
zplot = 0:0.1:zmax;
figure(3); clf;
    subplot(2,1,1);
        plot(zplot,fDSW(zplot));
end
