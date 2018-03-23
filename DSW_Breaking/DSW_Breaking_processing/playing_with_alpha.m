%% Write conversion factors in terms of unfitted and fitted values
load('/Volumes/Data Storage/2018_02_23_DSW_Breaking/quantities.mat',...
          'g','Delta','mui','bottomPixToCm','alpha','epsilon','dsw','L0fac','T0fac','myU');
load('/Volumes/Data Storage/2018_02_23_DSW_Breaking/Processing_Get_Breaking/fig_quants.mat')

alpham = alpha;
    muiDm  = mui/Delta;
    epsilonm = epsilon; % from original fluids
    Q0 = 0.25; % mL/min
    R0m = (alpham*Q0^0.25) / 2; % cm (alpha is in cg min units)
    Lm = R0m/sqrt(8);           % cm
    Um = g*R0m.^2/8./muiDm;   % cm/min
    Tm = sqrt(8)/g./R0m.*muiDm;    % min

alphat = 0.3:0.0001:0.34;
    muiDt = alphat.^4 * pi * g / 2^7; % mui/Delta ratio 
    R0t = (alphat*Q0^0.25) / 2; % cm (alpha is in cg min units)
    Lt = R0t/sqrt(8);           % cm
    Ut = g*R0t.^2/8./muiDt;   % cm/min
    Tt = sqrt(8)/g./R0t.*muiDt;    % min

%% Calculate expected tb's based on minimizing residuals
tb = [dsw.tb];
tb_res = mean([out([2:10,12:15]).tb] - [dsw([2:10,12:15]).tb]);
tb_fit = tb + tb_res; % this might have the wrong sign

epsilont = epsilonm * (Tt./Tm .* mean(tb)./mean(tb_fit)).^2;

zb = (15:25)'; % incorrect dimensional zb
zbt = repmat(zb,1,length(Lt)) .* repmat((Lt./sqrt(epsilont) * sqrt(epsilonm)/Lm),length(zb),1); % correct dimensional zb, based on fitted parameters

[ZB,ALPHAT] = meshgrid(zb,alphat);
figure(1); clf;
    contourf(ZB,ALPHAT,zbt',25,'edgecolor','none');
    colorbar;

%% GOAL: minimize zb error using alpha
%% (equivalent to using muiD)
Zcoeff = sqrt(Q0*8*g/(pi)) * Tm*tb_fit./(Lm*tb) ;
muiDfac = [out.zb]/(Zcoeff.*[in.zb]);

mDfac = (-1/2);
muiD_fit = muiDfac^mDfac;
disp(muiD_fit);
alpha_fit = ( (2^7*muiD_fit)/(pi*g) )^(1/4);
disp(alpha_fit);
%% Test results
zbfit = (Zcoeff.*[in.zb])*(muiD_fit)^(1/mDfac);
figure(2);
plot(zbfit,[out.zb],'o',[16,27],[16 27],'k--');


    
    