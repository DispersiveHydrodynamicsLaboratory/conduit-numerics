function [ g0, dg0 ] = triangle_BC( Ab, zb, toff, tmax )
%triangle_BCs generates boundary conditions 
%   that result in a triangle wave in the conduit
%   of maximum height Ab and width zb
%   toff is the amount of time (ND) before starting 
%   the ramping up in the profile
%   tmax is time profile is allowed to rundg0_on = 0;
%     zb = 10;

    t1 = zb*(Ab - 1)/(2*Ab);
    m  = -(Ab - 1)/zb;
    t2 = zb/2;
    tm = -1/(2*m);
    
    if t2>tm
        disp('Breaking occurs before profile finished');
    end
    

     g0 = @(t) ones(size((t-toff)))                      .* ((t-toff)<=0 | (t-toff)>= t2) + ...
               1./( 1 - 2*((t-toff))/zb )                .* ((t-toff)>0  & (t-toff)<t1) +...
               Ab*ones(size((t-toff)))                   .* ((t-toff)==t1) +...
               (zb*m+Ab)./(1+zb*m-2*(t-toff)*m) .* ((t-toff)>t1 & (t-toff)<t2); 
               
    dg0 = @(t) zeros(size((t-toff)))                            .* ((t-toff)<=0 | (t-toff)>= t2 | (t-toff)==t1) + ...
               (2/zb)./( 1 - 2*(t-toff)/zb ).^2                 .* ((t-toff)>0  & (t-toff)<t1) +...
               (-2*m)*(zb*m+Ab)./(1+zb*m-2*(t-toff)*m).^2 .* ((t-toff)>t1 & (t-toff)<t2);

    tvec = -zb/2:0.01:zb;
    
    if 0
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
    end
%     input('r');

% end

