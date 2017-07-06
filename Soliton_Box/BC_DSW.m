function [ A, Adot ]  = BC_DSW( t, t0 , z, z0, Aplus , Aminus, toffset)
% t = 0:0.01:50; z = 0; z0 = 20; Aplus = 3.5; Aminus = 1; t0 = z0/(2*Aplus);
% Summary of this function goes here
%   Detailed explanation goes here
% Change Coordinates
zeta = -(z-z0);
tau  = -(t-toffset*t0);
[Zeta,Tau] = meshgrid(zeta,tau);

% Set up vectors for tanh curve fitting
dt1 = 0.3;                       % relative distance across connection
dt2 = dt1*Aplus;
clear b
b(1,:) = (Zeta <= (2*Aminus-dt1/2).*Tau );
b(2,:) = (Zeta >  (2*Aminus-dt1/2).*Tau & Zeta <= (2*Aminus+dt1/2)*Tau);
b(3,:) = ((2*Aminus+dt1/2).*Tau<Zeta & Zeta < (2*Aplus-dt2/2)*Tau);
b(4,:) = (Zeta >= (2*Aplus-dt2/2) .*Tau & Zeta < (2*Aplus+dt2/2) *Tau);
b(5,:) = (Zeta>(2*Aplus+dt2/2).*Tau);
x0 = Tau(b(1.,:));
y0 = Aminus.*ones(size(x0));
x1 = Tau(b(3,:));
y1 = (Zeta./(2*Tau)) .*b(3,:)';
    y1 = y1(y1~=0);
    y1(isnan(y1)) = [];
x2 = Tau(b(5,:));
y2 = Aplus           .*b(5,:);
    y2 = y2(y2~=0);
    
%% Initial conditions for the fit
x01 = [ Aminus/0.1, (1+dt1/2)*Aminus/100,      x0(end),            y0(end)];%[0.1,pi*0.01,x1(1),y1(1)];
if Aplus < 3
    if (Aplus==2.5 && ismember(z0,[10 20]))
        x02 = [ -Aplus/0.1 ,  (1+dt2/2)*Aplus/100, mean([x1(end),x2(1)]), mean([y1(end),y2(1)])]; %[0.8,pi*0.01,x1(end),y1(end)];
    else
        x02 = [ -Aplus/0.02 ,  (1+dt2/2)*Aplus/100, mean([x1(end),x2(1)]), mean([y1(end),y2(1)])]; %[0.8,pi*0.01,x1(end),y1(end)];
    end
elseif (Aplus==4.5 && z0 == 15) || (Aplus==4 && z0==10)
    x02 = [ -1.1872    2.9732    mean([x1(end),x2(1)]), mean([y1(end),y2(1)])]; %[0.8,pi*0.01,x1(end),y1(end)];
elseif (Aplus==5 && z0 == 5)
    x02 = [ -0.5700   15.9492    mean([x1(end),x2(1)]), mean([y1(end),y2(1)])]; %[0.8,pi*0.01,x1(end),y1(end)];
elseif (Aplus>=3.5 && z0==20) || (Aplus==4 && z0 ==30) || (Aplus==3 && z0==30) || (ismember(Aplus,[3 4]) && z0==50)
    x02 = [ -0.5700   3.9492    mean([x1(end),x2(1)]), mean([y1(end),y2(1)]) ]; %[0.8,pi*0.01,x1(end),y1(end)];
else
    x02 = [ -Aplus/0.2 ,  (1+dt2/2)*Aplus/100, mean([x1(end),x2(1)]), mean([y1(end),y2(1)])]; %[0.8,pi*0.01,x1(end),y1(end)];
end

disp('Fit1');
[a1,b1,c1,d1] = connector_MM(x0,y0,x1,y1,x01,2);
    f1    = @(x)  a1*   log(cosh(b1*(x-c1))) + d1;
    f1dot = @(x)  -a1*b1*tanh(b1*(x-c1));
disp('Fit2');
[a2,b2,c2,d2] = connector_MM(x1,y1,x2,y2',x02,3);
    f2    = @(x)  a2*   log(cosh(b2*(x-c2))) + d2;
    f2dot = @(x)  -a2*b2*tanh(b2*(x-c2));
A = Aminus                                     .*b(1,:)' +...
    f1(Tau)                                    .*b(2,:)' +... % Correction for discontinuous first derivative
    (Zeta./(2*Tau))                            .*b(3,:)' +....
	f2(Tau)                                    .*b(4,:)' +... % Correction for discontinuous first derivative
    Aplus                                      .*b(5,:)';

Adot = 0                 .*b(1,:)' +...
       f1dot(Tau)        .*b(2,:)' +... % Correction for discontinuous first derivative
       (Zeta./(2*Tau.^2)).*b(3,:)' +...
       f2dot(Tau)        .*b(4,:)' +... % Correction for discontinuous first derivative
       0                 .*b(5,:)';

% Plot for debugging
figure(6); clf;
    subplot(2,1,1);
        plot(t,A,'*');
    subplot(2,1,2);
        plot(t,Adot,'*');


% end

