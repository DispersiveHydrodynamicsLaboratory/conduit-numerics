% Script for extracting leading edge speed from 
% Full image contour data
% Trials 21-3, 25, 31-33, 42, 45 good to watch ampl over time
% EDITED by MM to determine max oscillations for jump size

data_dir = 'F:/DATA/Photos/2015_06_30_and_07_01_and_07_02/dsw_full/';
process_data = 1;
save_on      = 1;
restr_osc    = 0; % Set nonzero to restrict number of oscillations to exp't
noise_on     = 1; % Set to nonzero to process DSWs with noisy ICs

% set_quantities;
trials     = [7:13,17:50]; 
ntrials    = length(trials);      
fontsize = 16;

expt_jumps   = zeros(size(trials));
expt_times   = zeros(size(trials));
expt_maxoscs = zeros(size(trials));
expt_speeds  = zeros(size(trials));
expt_amps    = zeros(size(trials));


for ii = 1:ntrials
    source_dir = [data_dir,'Trial',sprintf('%02d',trials(ii)),'/'];
%     load([source_dir,'contour_data.mat'],'D','source_dir','T',...
%          'num_files','vert','horiz','interpfactor','zconversion');
    load([source_dir,'full_camera_processed_data.mat'],...
            'speed','speed_err','amplitude','amplitude_err',...
            'jump','jump_err','background_err','finalt','maxosc');
    expt_jumps(ii)   = jump;
    expt_times(ii)   = finalt;
    expt_maxoscs(ii) = maxosc;
    expt_speeds(ii)  = speed;
    expt_amps(ii)    = amplitude;
end


% Script for extracting leading edge speed from 
% Full conduit equations numerics of DSWs
if restr_osc
    trials = [1.4 1.4 1.6 1.6 1.8 2 2 2.2 2.2 2.4 2.4 2.4 2.6 2.6 2.6,...
              2.8 2.8 2.8 2.8 3 3.2 3.4 3.4 3.4 3.4 3.8 3.8 3.8 4 4 4,...
              4.2 4.4 4.6];
    maxoscs = [2 2 3 4 4 5 4 6 8 6 5 7 10 6 10 7 9 8 7 8 7 8 9 10 11 9,...
               7 9 10 9 10 9 11 9];
else
    trials = [1.5:0.5:3,4,5,6,8];%1.4:0.2:4.6;
end
ntrials    = length(trials);   
    num_jumps   = zeros(size(trials));
    num_times   = zeros(size(trials));
    num_maxoscs = zeros(size(trials));
    num_speeds  = zeros(size(trials));
    num_amps    = zeros(size(trials));
fontsize = 16;
% DSW amplitude speed relation
Aminus = [1:0.01:20];
thy_speeds = @(z) sqrt(1+8*z)-1;
Asoli = zeros(size(Aminus));
Asoli(1) = 1;
Asoli(2) = fzero(@(z) thy_speeds(Aminus(2))-(z^2*(2*log(z)-1)+1)/(z-1)^2,...
                     1.01);
for ii=3:length(Asoli)
   Asoli(ii) = fzero(@(z) thy_speeds(Aminus(ii))-(z^2*(2*log(z)-1)+1)/(z-1)^2,...
                     Asoli(ii-1));
end
% Build interpolant for ease of evaluation
Adswpp = interp1(Aminus,Asoli,'spline','pp');
if noise_on
    noiseparam = '_with_noise';
else
    noiseparam = '';
end

for ind = 1:ntrials
    if noise_on
        data_dir = ['D:\Numerics\conduit_numerics_current/data/',...
                    'conduit_eqtn_tmax_250_zmax_1400_Nz_140000_order_4_',...
                    'init_condns_DSW_like_jump_',num2str(trials(ind)),...
                    '_with_better_noise_bndry_condns_time_dependent/'];
    else
        data_dir = ['D:\Numerics\conduit_numerics_current',...
            '/data/',...
            'conduit_eqtn_tmax_200_zmax_1200_Nz_120000_order_4_init_condns_',...
            'DSW_like_jump_',num2str(trials(ind)),...
            '_bndry_condns_time_dependent/'];
    end
    disp(['Trial: ',num2str(trials(ind))]);
    source_dir = data_dir;
%     source_dir = [data_dir,'Trial',sprintf('%02d',trials(ind)),'/'];
    load([source_dir,'parameters.mat']);
    if process_data
    z = dz:dz:zmax-dz; %ND
    
        tmax = floor(.6*length(t)); % reduces BC issues

    T = t(1:tmax);
    [X,Y] = meshgrid(z,T);
    D     = zeros(size(X));
    for ti = 1:length(T)
        load(strcat(source_dir,num2str(ti,'%05d')),'A','tnow');
        D(ti,:) = A;
    end
    D = fliplr(D);
    
    % Identify leading edge automatically.  Use findpeaks with a minimum
    % threshold peak height estimated from flow rate data
    % Assumes background is 1
    expected_jump = sqrt(trials(ind)); % Assumes Poiseulle flow
    expected_amplitude = ppval(Adswpp,expected_jump);
    % Use a reduced amplitude to be conservative
    expected_amplitude = 0.5*(expected_amplitude-1) + 1;
    % Rescale so that leading amplitude is 1
    expected_background = 1/expected_amplitude;
    expected_jump = expected_jump/expected_amplitude;
    % Scale leading amplitude to 1
    Dmax = max(D(:));
    d = (D/Dmax); % TRANSPOSE REMOVED TO ENSURE PARALLEL WITH EXP'T 
    A = d;        % SQUARE REMOVED BECAUSE ALREADY AREA
    % Now use findpeaks to identify leading edge
    dsw_position = [];
    dsw_time = [];
    dsw_amplitude = [];
    ctr = 1;
    for jj=1:length(T)
%        figure(1);
%            clf();
%            plot(z,A(jj,:),'b-',...
%                 z(locs),pks,'r*'); drawnow;
%             pause(0.25);
       [pks,locs] = findpeaks(A(jj,:),...
           'minpeakheight',0.5,...
           'minpeakprominence',0.3); 
       if jj==length(T)
       disp(['Time: ',num2str(T(jj))]);
       disp(['Peaks Found: ',num2str(length(locs))]);
       end
       if ~isempty(pks)
           dsw_position(ctr) = z(locs(1));
           dsw_amplitude(ctr) = pks(1);
           dsw_time(ctr) = T(jj);
           ctr = ctr+1;
%            figure(1);
%            clf();
%            plot(z,A(jj,:),'b-',...
%                 z(locs),pks,'r*'); drawnow;
%             pause(0.25);
%            input('continue?');
           % Break out of loop if reached experimental number of
           % oscillations
           if restr_osc
                check = (length(pks) > maxoscs(ind));
           else
                check = (z(locs(1)) <= 0.1*z(end));
           end
           if check
               break;
           end
       end
    end
    finalt = T(jj);
    maxosc = length(locs);
    disp(['Final time: ',num2str(finalt),' Max Osc Reached: ',num2str(maxosc)]);
    num_maxoscs(ind) = maxosc;
    % Determine speed and error estimate by fitting to a line
    p = polyfit(dsw_time,dsw_position,1);
    speed = -p(1);
    % Get error in speed
    ybar = mean(dsw_position);
    xbar = mean(dsw_time);
    SSxx = sum((dsw_time-xbar).^2);
    SSxy = sum((dsw_time-xbar).*(dsw_position-ybar));
    m = SSxy/SSxx;
    b = ybar-m*xbar;
    SSe = sum((dsw_position-m*dsw_time-b).^2);
    syx = sqrt(SSe/(length(dsw_time)-2));
    speed_err = 2*syx/sqrt(SSxx);
    
    % Now estimate background conduit area
    cutoff = 50; % distance in cm away from fitted leading edge position
    inds = find(Y > dsw_time(1) & Y < dsw_time(end) & X < polyval(p,Y)-cutoff);
    Abackground = A(inds);
    Aplus = mean(Abackground,'omitnan');
    Aplusstd = std(Abackground,'omitnan');
    
    % Rescale everything by background conduit area
    A = A/Aplus;
    Aplusstd = Aplusstd/Aplus;
    dsw_amplitude = dsw_amplitude/Aplus;
    Abackground = Abackground/Aplus;
    expected_jump = expected_jump/Aplus;
    expected_amplitude = expected_amplitude/Aplus;
    Aplus = 1;
    
    % Leading edge amplitude, only use latter half of data
    Nt = length(dsw_time);
    tstartind = round(Nt/2);
    Adsw = mean(dsw_amplitude(tstartind:end),'omitnan');
    Adswstd = std(dsw_amplitude(tstartind:end),'omitnan');
    
    % Estimate trailing edge amplitude; use a narrow band of space time to
    % estimate this
    inds = find(Y >= dsw_time(1) & Y <= dsw_time(1)+50 & X > 0.9*X(end));
    Atrail = A(inds);
    Aminus = mean(Atrail,'omitnan');
    Aminusstd = std(Atrail,'omitnan');
    
    % Save results
    amplitude = Adsw;
    amplitude_err = 2*Adswstd;
    jump = Aminus;
    jump_err = 2*Aminusstd;
    background_err = 2*Aplusstd;
    
    figure(1);
    clf();
    subplot(2,2,1);
    plot(dsw_time,dsw_position,'k.',...
         dsw_time,polyval(p,dsw_time),'r--');
    xlabel('$t$','interpreter','latex');
    ylabel('$z_+$','interpreter','latex');
    title(['$s_+ = ',num2str(-p(1)),...
           ' \pm ',num2str(speed_err),'$ cm/s'],...
           'interpreter','latex');
    axis tight;
    subplot(2,2,2);
    plot(dsw_time,dsw_amplitude,'k.');
    hold on;
    errorbar([dsw_time(tstartind(1)),dsw_time(end)],...
              Adsw*[1,1],2*Adswstd*[1,1],'r-');
    %plot([dsw_time(1),dsw_time(end)],expected_amplitude*[1,1],'b--');
    hold off;
    xlabel('$t$ (s)','interpreter','latex');
    ylabel('$A_+$','interpreter','latex');
    title(['$A_+ = ',num2str(Adsw),' \pm ',...
            num2str(2*Adswstd),'$'],'interpreter','latex');
    axis tight;
    subplot(2,2,3);
    plot([1:length(Abackground)],Abackground,'k.');
    hold on;
    errorbar([1,length(Abackground)],Aplus*[1,1],2*Aplusstd*[1,1],'r-');
    hold off;
    ylabel('$A$','interpreter','latex');
    title(['Background conduit area = $',num2str(Aplus),' \pm ',...
    num2str(2*Aplusstd),'$'],'interpreter','latex');    
    axis tight;
    subplot(2,2,4);
    plot([1:length(Atrail)],Atrail,'k.');
    hold on;
    errorbar([1,length(Atrail)],Aminus*[1,1],2*Aminusstd*[1,1],'r-');
    %plot([1,length(Atrail)],expected_jump*[1,1],'b--');
    hold off;
    ylabel('$A$','interpreter','latex');
    title(['Trailing conduit area = $',num2str(Aminus),' \pm ',...
        num2str(2*Aminusstd),'$'],'interpreter','latex');
    axis tight;
    % Plot in landscape mode
    fig = gcf;
    paperUnits = get(fig, 'PaperUnits');
    set(fig,'PaperUnits','inches');
    set(fig,'paperorientation','landscape');
    paperSize = get(fig,'PaperSize');
    paperPosition = [.25 .25 paperSize - .25];
    set(fig,'PaperPosition', paperPosition);
    set(fig,'PaperUnits',paperUnits);
%     print('-dpdf',[source_dir,'data_fits.pdf']);
    drawnow;
%     input('Return')
%     % Plot contour now appropriately normalized
%     num_t_pts = 200;
%     num_z_pts = 200;
%     t = linspace(0,T(end),num_t_pts);
%     zp = linspace(z(1),z(end),num_z_pts);
%     [Xi,Yi] = meshgrid(zp,t);
%     a = interp2(X,Y,A,Xi,Yi,'spline');
%     zp = -zp;
%     zp = zp - zp(end);
        if restr_osc
            save(['./data/','full_camera_processed_data_trial_',num2str(ind),'.mat'],...
                'speed','speed_err','amplitude','amplitude_err',...
                'jump','jump_err','background_err','finalt','maxosc');
        else
            save([source_dir,'full_camera_processed_data.mat'],...
                'speed','speed_err','amplitude','amplitude_err',...
                'jump','jump_err','background_err','finalt','maxosc');
        end
    else
        if restr_osc
            load(['./data/','full_camera_processed_data_trial_',num2str(ind),'.mat'],...
                'speed','speed_err','amplitude','amplitude_err',...
                'jump','jump_err','background_err','finalt','maxosc');
        else
            load([source_dir,'full_camera_processed_data.mat'],...
                'speed','speed_err','amplitude','amplitude_err',...
                'jump','jump_err','background_err','finalt','maxosc');
        end
    end
        num_jumps(ind)   = jump;
        num_times(ind)   = finalt;
        num_maxoscs(ind) = maxosc;
        num_speeds(ind)  = speed;
        num_amps(ind)    = amplitude;
end
if save_on
if restr_osc
    save('dsw_soliton_numerics_restr.mat',...
         'num_jumps','num_times','num_maxoscs','num_speeds','num_amps');
else
    save('dsw_soliton_numerics_no_restr.mat',...
         'num_jumps','num_times','num_maxoscs','num_speeds','num_amps');
end
end
% Compare to Theory
% Whitham theory
thy_jumps = [1.01:0.01:2.3];
thy_amps = zeros(size(thy_jumps));
thy_speeds = sqrt(1+8*thy_jumps.^2)-1;
f = @(a) (2*a^2*log(a)-a^2+1)/(a-1)^2;
thy_amps(1) = sqrt(fzero(@(a) f(a)-thy_speeds(1),2));
for ii=2:length(thy_amps)
   thy_amps(ii) = sqrt(fzero(@(a) f(a)-thy_speeds(ii),thy_amps(ii-1)^2));
end
% Previous values for Numerics
% Conduit equation numerics
old_jumps = sqrt([1.5:0.5:3,4,5,6,8]);
old_amps  = sqrt([2.0903,3.4281,4.9479,6.6327,10.4537,14.7972,19.6816,30.8412]);
old_speeds = [2.5856,3.0878,3.5165,3.8918,4.5279,5.0533,5.5073,6.2578];

%display in diameter
num_jumps = sqrt(num_jumps);
num_amps  = sqrt(num_amps);

figure(4);clf;
% subplot(2,2,1);
%     plot(num_jumps,num_maxoscs,'*');%,...
%          %expt_jumps,expt_maxoscs,'o')
%     legend('Numerics','Experiment');
%     xlabel('Nondimensional Jump (Area)'); ylabel('Number of Oscillations');
subplot(2,2,2);
    plot(num_jumps,num_speeds,'*',...
         old_jumps,old_speeds,'md',...%expt_jumps,expt_speeds,'o',...
         thy_jumps,thy_speeds,'-')
    legend('Numerics','Old Numerics','Theory');
    xlabel('Nondimensional Jump (Diameter)'); ylabel('Leading Edge Speed');
subplot(2,2,3)
    plot(num_jumps,num_amps,'*',...
         old_jumps,old_amps,'md',...%expt_jumps,expt_amps,'o',...
         thy_jumps,thy_amps,'-')
    legend('Numerics','Old Numerics','Theory');
    xlabel('Nondimensional Jump (Diameter)'); ylabel('Leading Edge Diameter');
subplot(2,2,4)
    plot(num_amps,num_speeds,'*',...
         old_amps,old_speeds,'md',...%expt_amps,expt_speeds,'o',...
         thy_amps,thy_speeds,'-')
    legend('Numerics','Old Numerics','Theory');
    xlabel('Leading Edge Diameter'); ylabel('Leading Edge Speed');

