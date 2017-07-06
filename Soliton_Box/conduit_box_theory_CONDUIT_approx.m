    load('box_num_results.mat','ps','ampls','ws','asolis','N_num','ampls_num','boxspd');
     amplsa = [0.25 0.5 1];
     wsa = [50 100 200];
    asolisa = [3 5 7];
    [AMPLS,WS,ASOLIS] = meshgrid(amplsa,wsa,asolisa);
    dphi_num = ps;
    numruns = numel(dphi_num);
    N_box = zeros(1,numruns);
    dphi_box = zeros(1,numruns);
    % Function for number of solitons determination (conduit)
    %     \km^2 = \frac{1}{2}\left( \lambda - \frac{2}{\phibar} + \sqrt{\frac{\lambda}{\phibar}(4+\phibar\lambda)} \right)
    km = @(phi,lambda,z) 1/2*(lambda - 2./phi(z) + sqrt(lambda/phi(z)*(4+phi(z)/lambda)));
    
    for ii = 1:numruns
        %% Initial condition properties
        ampl  = AMPLS(jj)+1;
        boxmax     = WS(jj);
        asoli = ASOLIS(jj);
        % Generate IC
            vstretch = 5;
            x0 = 500;
            v = make_smooth_box(vstretch,ampl,boxmax,x0);
            v = @(z) v(z) + ones(size(z));



    end