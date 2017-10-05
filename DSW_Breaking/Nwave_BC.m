

AB = 2; % Jump Height (nd)
ab = 0.5; % Fall Depth
zb = 15; % break height (cm)
w  = 12; % wave width (cm)

beta = 1+w/zb*1/(AB-ab);
m = (ab - AB)/w;

t1 = (AB-1)*zb/(2*AB);
t2 = (w-(1-ab)*zb)/(2*ab);
t3 = w/(2) ;
tb = zb/2;
wcheck = (ab/AB)*(AB-zb+2*AB*zb);

if wcheck>w
    warning(['Violation of parameter restrictions. Consider increasing the wave width to at least ',num2str(wcheck)]);
end

areafun = @(t)    ( 1                                .*( t <= 0 ) +...
                    1./(1-2/zb.*t)                   .*( t >  0    &  t <  t1 ) +...
                    (m*zb+AB)./(1+m*zb-2*m*t)        .*( t >= t1   &  t <  t2 ) +...
                    (1-w/zb)./(1-2/zb.*t)            .*( t >= t2   &  t <  t3 ) +...
                    1                                .*( t >= t3              ) );
                
dt = 0.1;

disp(['t1  t2  t3']);
disp(num2str([t1 t2 t3]));
disp(['Wcheck: ',num2str(wcheck),' w: ',num2str(w)]);

figure(2);
plot(0:dt:(tb+1),areafun(0:dt:tb+1),'.');