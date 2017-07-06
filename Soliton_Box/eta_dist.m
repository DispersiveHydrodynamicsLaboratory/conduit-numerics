function[f,etaf] = eta_dist(v,xmin,xmax,amin,ampl)        
        %% Distribution in spectral parameter, eta
        %  v: initial condition *usually box-like)
        %  xmin/xmax: minimum/maximum of domain
        %  amin: considered minimum of box
        %  ampl: considered maximum of box
        
        eta_a = @(a) sqrt(a/12); % spectral parameter
        I2 = @(v,x,eta) sqrt(6)/(pi) * eta ./ sqrt(v(x) - 6*eta.^2); % Integrand, based on \tilde{k}
        
        % Find bounds based on spectral parameter
        % Where v_0(x) = a = 6*eta^2
        x1zero = @(eta) fzero(@(x) v(x)-eta_a(2*ampl),[xmin,xmax/2]);
        x2zero = @(eta) fzero(@(x) v(x)-eta_a(2*ampl),[xmax/2,xmax]);
        eta = eta_a(amin:10^(-3):2*ampl);
        x2   = zeros(size(eta));
        x1   = zeros(size(eta));
        % Generate function of eta for bounds
        for ii = 1:length(eta)
            x1(ii) = x1zero(eta(ii));
            x2(ii) = x2zero(eta(ii));
        end
        feta = zeros(size(eta));
        for ii = 1:length(eta)
            feta(ii) = integral(@(x) I2(v,x,eta(ii)),x1(ii),x2(ii)); %*1/(pi*sqrt(6))
        end
f = @(etaq) interp1(eta,feta,etaq,'spline',0);
etaf = eta;