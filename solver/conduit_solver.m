function conduit_solver( t, zmax, Nz, h, g0, dg0, g1, dg1, f, directory );
% Solves the scalar conduit equation using a FD method in space
% and Matlab's ode45 solver in time
% Solves:
%
%  A_t = A*P
%  A*P - (A^2*P_z)_z = -(A^2)_z
%
% Inputs:
%
% t        :  1D array of output times
% zmax, dz :  spatial grid parameters 
% h        :  order of finite difference method 
% g0, dg0  :  boundary condition in time at z=0 and its derivative
% g1, dg1  :  boundary condition in time at z=zmax and its derivative
% f        :  initial condition in space at t=0
% directory - output directory where data at each output timestep
%             is written to a file conduit#####.mat
%
% Outputs:  NONE except for data written to files
%
  global inc dir tout start;

  % Set global variables for ode solver output function
  tout = t;
  dir = directory;
  inc = 0;
  dt = 1e-1;

  % Setup grid
  dz = zmax/(Nz+1);
  z = dz*[1:Nz]';


  %
  % Output what we are about to do
  %
  disp(['Solving conduit eqtn.']);
  disp([' Time interval:  [', num2str(t(1)),...
        ',',num2str(t(end)),']']);
  
  %
  % Construct initial condition on spatial domain
  %
  A_init = f(z);
  save(strcat(dir,num2str(inc,'%05d')),'A_init','inc');
  inc = inc + 1;
  A = A_init; tnow = 0;
  save(strcat(dir,num2str(inc,'%05d')),'A','tnow','inc');
  disp(['Time = ',num2str(tout(inc)),', inc = ',int2str(inc),...
          '/',int2str(length(tout))]);
  
  %Set appropriate ODE solver options.
  opts = odeset('OutputFcn',@output,'Refine',1,'RelTol',1e-6, ...
                'AbsTol',1e-6,'Stats','on','maxstep',dt);
  start = tic;

  % Get soln using Matlab RK timestepper (myode23 or myode45 both work)
  if h ==2 % Use second order FD method
      %
      % Precompute differentiation matrices 
      %
      D1  = 1/(2*dz) * spdiags([ -ones(Nz,1),...
                                  zeros(Nz,1),...
                                  ones(Nz,1)],...
                               -1:1, Nz,Nz);
      D2  = 1/(dz^2) * spdiags([  ones(Nz,1),...
                               -2*ones(Nz,1),...
                                  ones(Nz,1) ],...
                               -1:1,Nz,Nz);
    ode45(@(t,y) compute_deriv2( t, y, dz, D1, D2, g0, dg0, g1, dg1 ),...
                                 tout, A_init, opts);
  elseif h ==4 % Use fourth order FD method
      %
      % Precompute differentiation matrices 
      %
      D1  = 1/(dz) * spdiags([ [ 1/12*ones(Nz-3,1);  1/2; 0; 0],...
                              [ -2/3 *ones(Nz-2,1); -3/2; 0],...
                        [ -5/6;       zeros(Nz-2,1); 5/6],...
                      [0;  3/2;  2/3 *ones(Nz-2,1)],...
                        [0;  0; -1/2; -1/12*ones(Nz-3,1)] ],...
                                -2:2,Nz,Nz);
      D1(Nz,Nz-3) = 1/(dz) * -1/12;
      D1(1,4)     = 1/(dz) *  1/12;
      D2  = 1/(dz^2) * spdiags([ [-1/12*ones(Nz-3,1);  1/3; 0; 0],...
                                 [ 4/3 *ones(Nz-2,1);  1/2; 0],...
                           [-5/3; -5/2 *ones(Nz-2,1); -5/3],...
                         [0; 1/2;  4/3 *ones(Nz-2,1)],...
                      [0; 0; 1/3; -1/12*ones(Nz-3,1)] ],...
                                -2:2,Nz,Nz);
      D2(Nz,Nz-3) = 1/(dz^2) * -1/12;
      D2(1,4)     = 1/(dz^2) * -1/12;
    ode45(@(t,y) compute_deriv4( t, y, dz, D1, D2, g0, dg0, g1, dg1 ),...
                                 tout, A_init, opts);
  end
  % Finish and clean up
  finish = toc(start);
  disp('Calculation Finished');
  days_left = datenum([0 0 0 0 0 toc(start)]);
  time_left=datevec(days_left-floor(days_left));
  disp(['Computation time = ',...
        int2str(floor(days_left)),'d ',...
        int2str(time_left(4)),'h ',...
        int2str(time_left(5)),'m ',...
        num2str(time_left(6)),'s']);
  
%
% Called by odesolver at every output time
%
function status = output(t,y,flag)
  global inc dir tout start;
  status = 0;
  if strcmp(flag,'')
    for ii=1:length(t)
      % Extract solution
      A = y(:,ii);
      tnow = t(ii);
      % Increment and output estimate of time left
      inc = inc + 1;
      % Output results to file for subsequent analysis
      save(strcat(dir,num2str(inc,'%05d')),'A','tnow','inc');
      percentage_of_work = (inc)/(length(tout));
      seconds_left = ((1/percentage_of_work)-1)*(toc(start));
      days_left = datenum([0 0 0 0 0 seconds_left]);
      time_left=datevec(days_left-floor(days_left));
      disp(['Time = ',num2str(tout(inc)),...
          ', inc = ',int2str(inc),'/',int2str(length(tout)),...
          ', computation time left = ',...
          int2str(floor(days_left)),'d ',...
          int2str(time_left(4)),'h ',...
          int2str(time_left(5)),'m ',...
          num2str(time_left(6)),'s']);
    end
  end





  