function [a,b,c,d] = connector_MM(x1,y1,x2,y2,x0,fignum)
%% Documentation
%--------------------------------------------------------------------------
% Purpose: This script fits a tanh a*tanh(b(x-c))+d to the pump rate 
%          profile to smooth the sharp corner(s)
%
% Author:  Dalton Anderson (edited by Michelle Maiden)
% Date:    June 2016
%--------------------------------------------------------------------------
% clear; close all; clc
% x1: z of initial conduit profile
% y1: corresponding conduit values
% x2: z values of what to connect to 
% y2: corresponding conduit values
% x0: best guess of (a,b,c,d) for a*log(cosh(b*(x-c)))+d;
% %% set default profiles
% % exponential
% x1 = -5:0.01:0.2;
% y1 = exp(x1);
% 
% % horizontal line
% x2 = 5:0.01:10;
% y2 = 2*ones(size(x2));

% derivitives
dydx1 = diff(y1)./diff(x1);
dydx2 = diff(y2)./diff(x2);

m1 = dydx1(end);

%% fit tanh between profiles
% x0 = [1,1,1,1]; % initial guess (a,b,c,d)
parameters = [x1(end),x2(1),y1(end),y2(1),dydx1(end),dydx2(1)];%dydx2(1)]; % conditions
fun = @(x) connect_fit(x,parameters); %
options = optimoptions('fsolve','Display','none','Algorithm','trust-region-reflective','StepTolerance',10^-8);

[A,fval,exitflag] = fsolve(fun,x0,options); % solve for optimal a,b,c,d values
n = 0;
while exitflag==0 && n<5
    x0 = A;
    n = n+1; disp(n);
    [A,fval,exitflag] = fsolve(fun,x0,options); % solve for optimal a,b,c,d values
end
if exitflag<=0
    disp(['Fit did not converge. Exitflag: ',num2str(exitflag),'.']);
end
a = A(1);b = A(2); c = A(3); d = A(4);
% % % Figure for debugging
% f = @(x) a*log(cosh(b*(x-c)))+d;
% X = linspace(x1(end),x2(1),100);
% figure(fignum);
%     plot(x1,y1,'k',x2,y2,'k');
%     hold on
%     plot(X,f(X),'*'); drawnow;
%     hold off
% xx = X;
% yy = f(X);
% % input('R');
end