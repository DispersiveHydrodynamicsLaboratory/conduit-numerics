function[zs,tmax] = soliton_position_fun(uminus,aminus,zminus,uplus);
%% Code to solve for a soliton's position 
%% INPUTS
% uminus 	% Initial conduit mean
% aminus    % Initial soliton amplitude
% zminus    % RW starts at 0, soliton starts at -zminus
% uplus     % Final conduit mean



%  as it travels up a rarefaction wave
debug_on = 0;
plot_on  = 0;
zminus = -zminus;  % Initial soliton position


% Calculate q0 = q(c_s(a_-,u_-),u_-) based on initial soliton amplitude and initial conduit mean
csoli  = @(a,m) m./a.^2 .* ( (a+m).^2 .* (2*log(1+a./m)-1) + m.^2);
q      = @(cs,m) cs*(cs+2*m)/(4*m);
q0     = q(csoli(aminus,uminus),uminus);

% Ensure tunneling can happen
cminus = csoli(aminus,uminus);
cDSW   = -uminus + sqrt(uminus^2 + 8*uplus*uminus);
if cminus <= cDSW
    disp('Tunneling criterion is not satisfied.');
    disp('Soliton will be trapped');
    trapped = 1;
else
    trapped = 0;
end

% ODE right-hand side to be solved for the RW
cs = @(q,m) -m+sqrt(m*(m+4*q0));

% First part: calculate soliton position before interacting with the RW
    zsM = @(t) cs(q0,1)*t + zminus;

    % Calculate where soliton catches up with RW
    t1 = zminus/(2-cs(q0,1));
    z1 = 2*t1;
    disp(['Soliton reaches trailing edge of RW at t=',num2str(round(t1))]);
    
% Second part: calculate soliton position throughout the RW
if ~trapped
    % Calculate maximum possible time spent in the RW (based on initial speed)
    t2est = zminus/(2*uplus-cs(q0,1));
    if t2est<t1
        % Possible for this to not work
        disp('Initial guess: tunneling ends at t2 = t1 + 100');
        t2est = t1+100;
    else
        disp(['Initial guess: tunneling ends at t2 = ',num2str(round(t2est))]);
    end
else
    t2est = t1+100;
    disp(['Will track soliton through RW for 100 (ND time), for a total time of ',num2str(round(t2est))]);
end

    % Use ODE 45 to solve the "hard part"- RW
    [tout,zsout] = ode45(@(t,zs) cs(q0,zs/(2*t)), [t1,t2est], z1);
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
    if ~trapped
        % Find where this solution actually intersects leading edge of RW
        t2 = fzero(@(t) zsRW(t)-2*uplus*t, t2est);
            % If had to extrapolate madly, run ODE solver again
            while (t2-t2est) > 1
                t2est = t2+10;
                [tout,zsout] = ode45(@(t,zs) cs(q0,zs/(2*t)), [t1,t2est], z1);
                zsRW = @(t) interp1(tout,zsout,t,'spline','extrap');
                t2 = fzero(@(t) zsRW(t)-2*uplus*t, t2est);
            end
            zsRW = @(t) interp1(tout,zsout,t,'spline',0);
            z2 = 2*uplus*t2;
            disp(['Soliton is through the RW at t=',num2str(round(t2))]);
        % Third part: calculate soliton position after interaction with RW
            zsP = @(t) cs(q0,uplus)*(t-t2) + z2;
    else
        t2 = t2est; 
        zsRW = @(t) interp1(tout,zsout,t,'spline','extrap');
    end
   

% Combine all parts together to get full solution
if ~trapped
    zs = @(t) zsM(t) .* (t<=t1)       +...
              zsRW(t).* (t>t1 & t<t2) +...
              zsP(t) .* (t>=t2);
else
    zs = @(t) zsM(t) .* (t<=t1)       +...
              zsRW(t).* (t>=t1);
end
      
if plot_on
    % Plot results
    figure(2); clf;
        tplot = 0:0.01:(t2+75);
        plot(tplot,zs(tplot));
        hold on;
            plot(t1,z1,'b*',t2,z2,'rx');
        xlabel('t'); ylabel('z_s(t)');
end

tmax = t2+75;
