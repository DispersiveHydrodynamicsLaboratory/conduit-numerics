% Try to explain conduit SoliBox with KdV SoliBox Theory
save_on = 1;
plot_debug_on = 0;

if save_on 
    load('box_num_results.mat','ps','ampls','ws','asolis','N_num','ampls_num','boxspd');
     amplsa = [0.25 0.5 1];
     wsa = [50 100 200];
    asolisa = [3 5 7];
    [AMPLS,WS,ASOLIS] = meshgrid(amplsa,wsa,asolisa);
    dphi_num = ps;
    numruns = numel(dphi_num);
    N_box = zeros(1,numruns);
    dphi_box = zeros(1,numruns);
    dphi_gas = zeros(1,numruns);
    boxspd_gas = zeros(1,numruns);
    for jj=1:numruns
        disp(jj);
        % Amplitude distribution function for a box
        %% Box problem asymptotics
        %% Initial condition properties
        ampl  = AMPLS(jj)+1;
        boxmax     = WS(jj);
        asoli = ASOLIS(jj);
        % Generate IC
            vstretch = 5;
            x0 = 500;
            v = make_smooth_box(vstretch,ampl,boxmax,x0);
            v = @(x) v(x);
            xmin = 0;       % Domain start
            xmax = 1000;     % Domain end
        % IC integration conditions
        bmin = fzero(@(x) v(x)-(10^(-3)),[xmin,x0]); % z where box begins
        bmax = fzero(@(x) v(x)-(10^(-3)),[x0,xmax]); % z where box ends
        amin = v(bmin);     % considered amplitude minimum
        amax = 2*ampl-10^(-6);        % maximum soliton, approx. leading edge of DSW of same jump

        %% Determination of number of solitons
        km = @(v,x,lambda) sqrt(2/3*(v(x) - lambda)); %k_-, trailing edge wavenum
        N = round(1/(sqrt(6)*pi)*integral(@(x) km(v,x,0),bmin,bmax) + 1/2); % integral %

        %% Distribution in spectral parameter, eta
        I2 = @(v,x,eta) sqrt(6)/(pi) * eta ./ sqrt(v(x) - 6*eta.^2) ;            % Integrand, based on \tilde{k}
%         I3 = @(v,x,eta) sqrt(v(x)-eta.^2);     % Integrand whose derivative gives I2
        % Find bounds based on spectral parameter
        % Where v_0(x) = a = 2*eta^2
        eta_a = @(a) sqrt(a/12);
        x1zero = @(eta) fzero(@(x) v(x)-6*eta.^2,[xmin,xmax/2]);
        x2zero = @(eta) fzero(@(x) v(x)-6*eta.^2,[xmax/2,xmax]);
        eta = eta_a(amin:10^(-3):2*ampl-10^(-3));
        x2   = zeros(size(eta));
        x1   = zeros(size(eta));
        feta = zeros(size(eta));

        % Generate function of eta for bounds, integrate over said bounds
        for ii = 1:length(eta)
            x1(ii) = x1zero(eta(ii));
            x2(ii) = x2zero(eta(ii));
            feta(ii) = integral(@(x) I2(v,x,eta(ii)),x1(ii),x2(ii));
        end
        
        if plot_debug_on
            % Check where integration bounds are
            figure(1); clf;subplot(3,1,1);
            plot(xmin:10^(-2):xmax,v(xmin:10^(-2):xmax),...
                 x1,v(x1),'b*',x2,v(x2),'r*'); drawnow;
        end
        
        % Plot comparing thy soliton distribution to num soli dist
        num_etas = eta_a(ampls_num.(['run',num2str(jj)])/2);
        f = @(etaq) interp1(eta,feta,etaq,'spline',0);
        if plot_debug_on
            figure(1); 
                subplot(3,1,2);
                fnorm = integral(@(eta) f(eta),min(eta),max(eta));
                plot(eta,f(eta)/fnorm);
                title('Comparison of normalized f(\eta) to normalized histogram of soliton amplitudes');
                hold on; histogram(num_etas,'Normalization','pdf'); hold off;
        end

        % Use the eta distribution to find the phase shift of a soliton 
        % going through the box
        ps_thy = @(eta1,eta2) 1./eta1.*log( abs((eta1+eta2)./(eta1-eta2))); % 2-soli phase shift
        I1 = @(eta0, eta) ps_thy(eta0, eta).*f(eta); % integrand
        eta0 = eta_a(asoli);       % Soliton spectral parameter
        eta1 = 0;               % lower bound of integration
        eta2 = eta_a(2*ampl); % upper bound of integration
            
            dphi = integral(@(eta) I1(eta0,eta), eta1, eta2);
        
        
        
        N_box(jj) = N;
        dphi_box(jj) = dphi;
        
        disp(['Predictions for box with width ',num2str(boxmax),...
              ' and ampl ',num2str(ampl)]);
        disp(['Box, N: ',num2str(N),' Box, Phase Shift: ',num2str(dphi)]);
        
        % Plot looking at I1 integrand characteristics
        etaplot  = eta1:10^-3:eta2;
        if plot_debug_on
            figure(1); 
            subplot(3,1,3);
                plot(etaplot,I1(eta0,etaplot));
                title('Phase Shift Integrand Characteristics');
                drawnow; pause(0.1);
        end
    end
    save('ps_num_thy.mat','dphi_num','dphi_box','ampls','ws','asolis');
else
    load('ps_num_thy.mat','dphi_num','dphi_box','ampls','ws','asolis');
end 
% Now, assume box is a 1-component gas to find the effective speed,
        % phase shift due to the 'gas.' 
        ps_thy = @(eta1,eta2) 1./eta1.*log( abs((eta1+eta2)./(eta1-eta2))); % 2-soli phase shift
        s = @(eta1,eta0,f0) 4*( eta1.^2 - ps_thy(eta1,eta0).*eta0.^2.*f0(eta0) )./...
                              (    1    - ps_thy(eta1,eta0)        .*f0(eta0) );
%         etas = eta_a(assoli);   % Spectral parameter of soliton
%         etaf = eta_a(2*ampl); % Spectral parameter of gas component
        f0   = @(eta0) eta0/3;
        % critical density
        dphi_gas = ws.*( 1-4*eta_a(asolis).^2./s(eta_a(asolis),eta_a(2*ampls),f0) );
        save('ps_num_thy.mat','dphi_gas','-append');
%         % DSW theory density
%         dphi_gaspi = ws.*( 1-4*eta_a(asolis).^2./s(eta_a(asolis),eta_a(2*ampls),@(eta0)eta0/pi) );
pbox = polyfit(dphi_num,dphi_box,1);
pgas = polyfit(dphi_num,dphi_gas,1);
disp(['BoxProb fit: ',num2str(pbox(1)),'x + ',num2str(pbox(2))]);
disp(['Soligas fit: ',num2str(pgas(1)),'x + ',num2str(pgas(2))]);

figure(2); clf;
        h1 = plot(dphi_num,dphi_box,'bo',...
                  dphi_num,dphi_gas,'ks',...
                  dphi_num,dphi_num,'r-',...
                  'MarkerSize',8);
            set(h1(1),'MarkerFaceColor','b');
            set(h1(2),'MarkerFaceColor','k');
            xlabel('Phase Shift (from numerics)');
            ylabel('Phase Shift (from theory)');
            legend('Box Prediction','Gas Prediction','y=x');
            set(gca,'FontSize',12);
            title('Phase shift of test soliton through a box');
%     figure(3);
%         plot(N_num,N_box,'bo',...
%              N_box,N_box,'r-','MarkerFaceColor','b');
%          title('Number of Solitons in the Box');
%          xlabel('Numerics'); ylabel('Theory');
figure(4); clf;
    h2 = plot(boxspd,s(eta_a(asolis+1),eta_a(2*ampls),f0),'ks',...
              boxspd,boxspd,'r-');
          set(h2(1),'MarkerFaceColor',ones(1,3)*0.5,'MarkerSize',8);
         xlabel('Test soliton speed (from numerics)');
         ylabel('Test soliton speed (from theory)');
         legend('Gas Prediction','y=x');