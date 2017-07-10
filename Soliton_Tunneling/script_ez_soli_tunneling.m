save_on  = 1;  % Set to nonzero if you want to run the solver, set
                % to 0 if you want to plot
periodic = 0; % set to nonzero to run periodic solver (no BCs need)
              % set to 0 to run solver with time-dependent BCs      
check_IC = 0; % Set to 1 to only plot ICs
plot_on  = 0;  % Set to 1 if you want to plot just before and just
                % after (possibly) calling the solver
% load('soli_pairs.mat','aDSW','aRW','m_minus','m_plus');
aDSW = [5.25 6.25] - 1.75;
aRW  = [7.25 8.53];
m_minus = 1;
m_plus = 1.75;
if m_minus==1
    ADSW = m_plus;
else
    error('Please check jump criteria. Non-unity background detected.');
end
%% Numerical Parameters
tmax     = 125;    % Solver will run from t=0 to t=tmax
zmax     = 800;     % Solver will solve on domain z=0 to z=zmax
numout   = round(tmax) ;           % Number of output times
t        = linspace(0,tmax,numout);  % Desired output times
dzinit   = 1/100; % Set to 1/500 for optimum
Nz       = round(zmax/dzinit);
h        = 4;

if periodic
    dz       = zmax/Nz;    % Spatial  discretization
else
    dz       = zmax/(Nz+1);    % Spatial  discretization
end

for ind = 1:length(aDSW)          
    disp(['Soliton pair: ',num2str(aDSW(ind)),' ',num2str(aRW(ind))]);
    driver_conduit_solver_soli_tunneling_DSW;
        if plot_on && ~save_on
            input('r');
        end
%     driver_conduit_solver_soli_tunneling_RW;
%             if plot_on && ~save_on
%             input('r');
%             end
end