% function [ g0, dg0 ] = triangle_BC( Ab, zb, toff, tmax )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
dg0_on = 0;
    zb = 10;
    Ab = 1.5;
%     zb = 10;

    t1 = zb*(Ab - 1)/(2*Ab);
    m  = -(Ab - 1)/zb;
    t2 = zb/2;
    tm = -1/(2*m);
    
    if t2>tm
        disp('Breaking occurs before profile finished');
    end
    

     g0 = @(t) ones(size(t))                      .* (t<=0 | t>= t2) + ...
               1./( 1 - 2*(t)/zb )                .* (t>0  & t<t1) +...
               Ab*ones(size(t))                   .* (t==t1) +...
               (zb*m+Ab)./(1+zb*m-2*t*m) .* (t>t1 & t<t2); 
               
    dg0 = @(t) zeros(size(t))                            .* (t<=0 | t>= t2 | t==t1) + ...
               (2/zb)./( 1 - 2*t/zb ).^2                 .* (t>0  & t<t1) +...
               -(2/zb)./( (Ab - 2)/(Ab - 1) - 2*t/zb ).^2 .* (t>t1 & t<t2);

    tvec = -zb/2:0.01:zb;
    figure(1); clf; 
    if dg0_on
        subplot(2,1,1);
    end
        plot(tvec,g0(tvec),'k-'); title(zb); 
        axis([-zb/2 max(tvec) min(g0(tvec))-0.25 max(g0(tvec))+0.25])
        hold on;
            plot(t1*ones(1,2),[min(g0(tvec))-0.25 max(g0(tvec))+0.25],'r--',...
                 t2*ones(1,2),[min(g0(tvec))-0.25 max(g0(tvec))+0.25],'b--');
        hold off;
	if dg0_on
        subplot(2,1,2);
            plot(tvec,dg0(tvec),'k-');
            axis([-zb/2 max(tvec) min(dg0(tvec))-0.25 max(dg0(tvec))+0.25])
            hold on;
                plot(t1*ones(1,2),[min(dg0(tvec))-0.25 max(g0(tvec))+0.25],'r--',...
                     t2*ones(1,2),[min(dg0(tvec))-0.25 max(dg0(tvec))+0.25],'b--');
            hold off;
    end
	drawnow;
%     input('r');

% end

