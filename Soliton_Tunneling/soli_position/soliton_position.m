%% Code to solve for a soliton's position 
%  as it travels up a rarefaction wave
debug_on = 1;
% PARAMETERS
uminus = 1;    % Initial conduit mean
aminus = 3;    % Initial soliton amplitude
zminus = -50;  % Initial soliton position
uplus  = 2; % Final conduit mean

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
t1 = zminus./(3-sqrt(1+4*q0));
z1 = 2*t1;
% Calculate maximum possible time spent in the RW (based on initial speed)
t2 = zminus./(2*uplus + 1 - sqrt(1+4*q0));
if t2<t1
    disp('Bad end time. Using t2 = t1 + 100');
    t2 = t1+100;
end

% Use ODE 45 to solve the "hard part"- RW
[tout,zsout] = ode45(@(t,zs) dzsdt(t,zs), [t1,t2], z1);

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


% Combine with the "easy parts-" constant conduit

