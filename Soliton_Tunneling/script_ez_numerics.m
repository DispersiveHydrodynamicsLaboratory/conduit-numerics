save_on  = 1;  % Set to nonzero if you want to run the solver, set
                % to 0 if you want to plot
periodic = 0; % set to nonzero to run periodic solver (no BCs need)
              % set to 0 to run solver with time-dependent BCs      
check_IC = 1; % Set to 1 to only plot ICs
plot_on  = 1;  % Set to 1 if you want to plot just before and just
                % after (possibly) calling the solver

% Choose parameters here
asolis = 2:2:10; 
m_minus = 1; % note this one is the same throughout
m_pluses = 1.75:0.25:2.5;

% Make a meshgrid of variables that change
% If only one variable is changing, skip this
% If have more than 2 variables change, usde ndgrid
[ASOLIS, MPLUSES] = meshgrid(asolis,m_pluses);

% Set up for loop to run over all options
for ind = 1:numel(ASOLIS)
    %% Choose parameters based on loop iteration
    aDSW = ASOLIS(ind);
    m_plus = MPLUSES(ind);
    
    %% Numerical Parameters--can be made to change based on parameters
        % Example of using if-then statements to do so
            if aDSW < 3 
                tmax = 125;    % Solver will run from t=0 to t=tmax
            else
                tmax = 500;
            end
        % Example of using a parameter value to do so
            zmax     = 800*m_plus;     % Solver will solve on domain z=0 to z=zmax

    disp(['Soliton: ',num2str(aDSW),' Jump: ',num2str(m_plus)]);
	driver_conduit_solver_ez_exm;
    % Can also use this setup for other scripts, provided they use 
    % the same format for data directory
%     soli_tunneling_find_solis;
end