function[f] = soli_tunneling_IC_RW( A, asoli, zmax, hstretch, zjump, z0, wave_type)
debug_on = 0;
% Generate jump conditions
    vscale   = ((A-1)/2); 
    upshift  = 1+vscale;
        wave       = @(z) upshift + vscale*tanh(1/hstretch*(z-zjump));
     soli = @(z) zeros(size(z));
% Putting it all together
    f  = @(z)  wave(z) + soli(z);

if debug_on
% Figure for debugging
zplot = 0:0.1:zmax;
figure(3); clf;
        plot(zplot, f(zplot));
end
