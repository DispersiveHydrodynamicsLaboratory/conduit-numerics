function[zs] = soliton_position_fun(uminus,aminus,zminus,uplus);
%% Code to solve for a soliton's position 
%  as it travels up a rarefaction wave
debug_on = 0;
plot_on  = 0;

% % PARAMETERS
% uminus = 1;    % Initial conduit mean
% aminus = 3;    % Initial soliton amplitude
% zminus = -50;  % Initial soliton position
% uplus  = 2; % Final conduit mean

% Calculate q0 = q(c_s(a_-,u_-),u_-) based on initial soliton amplitude and initial conduit mean
csoli  = @(a,m) m./a.^2 .* ( (a+m).^2 .* (2*log(1+a./m)-1) + m.^2);
q      = @(cs,m) (cs.^2 + 2*m.*cs)./(4*m);
q0     = q(csoli(aminus,uminus),uminus);

% Ensure tunneling can happen
cminus = csoli(aminus,uminus);
cDSW   = -uminus + sqrt(uminus^2 + 8*uplus*uminus);
if cminus <= cDSW
    disp('Tunneling criterion is not satisfied.');
    disp('This is not going to work...');
    input('Return to do it anyways.');
end

% ODE right-hand side to be solved for the RW
dzsdt = @(t,zs) -1+sqrt(1+4*q0*zs/t);

% First part: calculate soliton position before interacting with the RW
    zsM = @(t) (-1+sqrt(1+4*q0))*t+zminus;

    % Calculate where soliton catches up with RW
    t1 = zminus./(3-sqrt(1+4*q0));
    z1 = 2*t1;
    
% Second part: calculate soliton position throughout the RW
    % Calculate maximum possible time spent in the RW (based on initial speed)
    t2est = zminus./(2*uplus + 1 - sqrt(1+4*q0));
    if t2est<t1
        % Possible for this to not work
        disp('Bad end time. Using t2 = t1 + 100');
        t2est = t1+100;
    end

    % Use ODE 45 to solve the "hard part"- RW
    [tout,zsout] = ode45(@(t,zs) dzsdt(t,zs), [t1,t2est], z1);
    zsRW = @(t) interp1(tout,zsout,t,'spline','extrap');

    if debug_on
        % Plot for debugging
        figure(1); clf;
            plot(tout,zsout);
            if sum(imag(zsout))>0
                hold on;
                    plot(tout,imag(zsout),'r--');
                hold off
            end
    end

    % Find where this solution actually intersects leading edge of RW
    t2 = fzero(@(t) zsRW(t)-2*uplus*t, t2est);
        % If had to extrapolate madly, run ODE solver again
        while (t2-t2est) > 1
            t2est = t2+10;
            [tout,zsout] = ode45(@(t,zs) dzsdt(t,zs), [t1,t2est], z1);
            zsRW = @(t) interp1(tout,zsout,t,'spline','extrap');
            t2 = fzero(@(t) zsRW(t)-2*uplus*t, t2est);
        end
        zsRW = @(t) interp1(tout,zsout,t,'spline',0);
    
% Third part: calculate soliton position after interaction with RW
    zsP = @(t) (-1 + sqrt(1+4*uplus*q0))*(t-t2) + zsRW(t2);
   

% Combine all parts together to get full solution
zs = @(t) zsM(t) .* (t<=t1)       +...
          zsRW(t).* (t>t1 & t<t2) +...
          zsP(t) .* (t>=t2);
      
if plot_on
    % Plot results
    figure(2); clf;
        tplot = 0:0.01:(t2+50);
        plot(tplot,zs(tplot));
        xlabel('t'); ylabel('z_s(t)');
end
