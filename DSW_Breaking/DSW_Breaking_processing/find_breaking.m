%% Breakheight/Times Determination
function[Abreak] = find_breaking(A,z,t,Aback,Afront,z0,t0,toffset,plot_slopes_on, plot_debug_on, is_numerics)
        t        = t-t(1);     % Time starts at 0
        numpics  = length(t);  % numpics
        zmax     = max(z);
        zmin     = min(z);
        Amax     = max(A(:));
        t0 = t0*toffset; % Expected breaking time (RIGHT)
    % Determine slope between front and back
    Aslope = [];
    tslope = [];
    thresh = 0.1; % Changed for experiments; value for numerics is 10^(-1);
    
    for ri = 1:numpics
        if ri==numpics
            disp('');
        end
        if plot_slopes_on
        figure(6); clf;
                plot(z,A(ri,:)); 
                axis([0 zmax 1 Amax]);
        end
%         Skip if back conduit isn't large enough yet
        if (Aback - max(A(ri,:))) > thresh
            if plot_slopes_on
                disp(['At time: ', num2str(t(ri)),' Slopefinder broke due to threshold violation']);
            end
            continue;
        end
        % Find lower bound for z+/- (wave front)
        [peak , locp] = findpeaks(fliplr(A(ri,10:end-10)),z(10:end-10),'MinPeakProminence',0.2,'SortStr','none','NPeaks',1);
            locp = zmin+zmax-locp;
        if isempty(peak) % peak is too small
            [peak, ip] = max(A(ri,:));
            locp = z(ip);
        end
        if plot_slopes_on
            figure(6); hold on;
                plot(locp,peak,'r*'); hold off;
                pause(0.1);
                drawnow;            
        end
        %% Following bound is changed for experiments, as near upper camera boundary,
        %% not near conduit boundary
        if (zmax-locp) < 3*Aback % wave is near upper boundary
            break;
        end
%         if abs(mean(A(ri,end-10:end))-Afront)>0.5
%             if plot_slopes_on
%                 disp(['At time: ', num2str(t(ri)),'DSW is leaving imaging field']);
%             end
%             continue;
%         end

        % Calculate slope of difference
        % To determine z+/-, cut off wave at first peak
            bndminus = 0.8*peak; %Aback;%max(A(ri,:));%0.9*Aback;
            bndplus  = 1+0.2*peak; %(1+0.2*Aback);
            zbnds    = z(z>locp);
            Abnds    = A(ri,z>locp);
            Afunbnds = @(zq) interp1(zbnds,Abnds,zq,'spline'); % Extrapolation value altered from numerics
            % Initial guesses 
            zminus0 = min([1.25*locp zmax]);
            if isempty(locp)
                zminus0 = mean([zmin zmax]);
            end
            [zminus, foo, exitflag1] = fzero(@(z) Afunbnds(z)-bndminus,[min(zbnds) max(zbnds)]);
            if exitflag1~=1
                disp('fzero problem');
            end
            [zplus, foo, exitflag2] = fzero(@(z) Afunbnds(z)-bndplus,[min(zbnds) max(zbnds)]);
            if exitflag2~=1
                disp('fzero problem');
            end
        % If z+/- is too far from value, skip
        if abs(Afunbnds(zminus)-bndminus)>0.1 || abs(Afunbnds(zplus)-bndplus)>0.1  || exitflag2~=1
            disp(['At time: ', num2str(t(ri)),' Slopefinder broke due to inaccurate zeros']);
            continue;
        end;
        if plot_slopes_on && mod(ri,50)==0
            disp(['Index: ',num2str(ri), ' Time: ',num2str(t(ri))]);
        end
        % Slope calculations
            Afun = @(zq) interp1(z,A(ri,:),zq,'spline',1);
            Aslope = [ Aslope, (Afun(zminus) - Afun(zplus))/(zminus - zplus) ];
            tslope = [ tslope, t(ri) ];
%         % Figure for debugging
          if plot_slopes_on
            figure(6); hold on;
                plot([zminus zplus], Afun([zminus zplus]),'b.','MarkerSize',8);
                plot(locp,peak,'r*');
                title(['Time: ', num2str(t(ri)),' Slope: ',num2str(Aslope(end))]);
                drawnow; pause(0.5);
          end
        if abs(Aslope(end)-mean(Aslope))>10^3 || isnan(Aslope(end))
            disp('Uh-oh');
        end
    end
    told = t;
    t = t(t<tslope(end));
	%% Interpolate breakheight-determining function
        newt = linspace(min(tslope), max(tslope), 200);
        if ~is_numerics
         Aslope = wiener2(Aslope,[3,3]); % Smooth data first
        end
         Aslopefun = csapi(tslope,Aslope);    % Function for cubic spline approx
         As = @(tq) fnval(Aslopefun,tq);    % Function evaluation
        dAs = @(tq) fnval(fnder(Aslopefun),tq); % First Derivative
       ddAs = @(tq) fnval(fnder(fnder(Aslopefun)),tq);   % Second Derivative
       
       %% Determine breaking time based on data
       if ~is_numerics
           newt1 = (newt(end)+newt(1))/2;
       else
           newt1  = newt(1);
       end
       Abreak.time = fminbnd(dAs,newt1,newt(end));
        
    
%         time2 = newt(dAs(newt)==min(dAs(newt)));
%         if abs(time1-time2)>1
%             Abreak.time = time2;
%         else
%             Abreak.time = time1;
%         end
        
        
       % Plot breaktime-determining function and its derivatives (from interpolation)
       if plot_debug_on
           figure(7); clf;
           title(['Aback: ',num2str(Aback),' zb: ',num2str(z0)]);
           subplot(3,1,1);
            plot(newt, As(newt));
            hold on; 
                plot([t0 t0], [min(As(newt)) max(As(newt))], '-'); 
                plot([Abreak.time Abreak.time], [min(As(newt)) max(As(newt))], '-');
                    legend('A_{slope}(t)','t_b Expected','t_b Observed');
            hold off;
           subplot(3,1,2);
            plot(newt,dAs(newt));
            hold on; 
                plot([t0 t0], [min(dAs(newt)) max(dAs(newt))], '-'); 
                plot([Abreak.time Abreak.time], [min(dAs(newt)) max(dAs(newt))], '-');
            hold off;
           subplot(3,1,3);
            plot(newt,ddAs(newt));
            hold on; 
                plot([t0 t0], [min(ddAs(newt)) max(ddAs(newt))], '-'); 
                plot([Abreak.time Abreak.time], [min(ddAs(newt)) max(ddAs(newt))], '-');
            hold off;
            drawnow;
       end
       
        %% Interpolate area function to breaking time
        Abreaktime    = interp2(z,told,A,z,Abreak.time);
        %% Convert to function; look at derivatives
        if ~is_numerics
         Abreaktime = smooth(z,Abreaktime,15); % Smooth data first
        end
         Abreakfun = csapi(z,Abreaktime);    % Function for cubic spline approx
         Ab = @(zq) fnval(Abreakfun,              zq);   % Function evaluation
%         dAb = @(zq) fnval(fnder(Abreakfun),       zq);   % First Derivative
%         % For experiments, Aslopefun is very messy, causing a lot of problems
%         % Will try to correct for this by smoothing derivative function as well
%         % as Aslope
%         dAbvec = smooth(dAb(z),5);
%         dAb    = @(zq) fnval(csapi(z,dAbvec),zq);
%        ddAb = @(zq) fnval(fnder(fnder(Abreakfun)),zq);   % Second Derivative (too messy to be useful)

        %% Determine breaking height based on data at breaking time
        % look ahead of peak
            [peak , locp] = findpeaks(fliplr(Abreaktime),z,'MinPeakProminence',0.1,'SortStr','none','NPeaks',1);
%             if ~is_numerics % experiment is backwards
                locp = zmin+zmax-locp;
%             end
        breakarea = Aback;%(peak+1)/2;
        Abreak.height = fzero(@(z) Ab(z)-breakarea,[locp max(z)]);
        
        
        
        % Plot breaktime-area function and its derivatives (from interpolation)
        if plot_debug_on
           figure(8); clf;
           title(['Aback: ',num2str(Aback),' zb: ',num2str(z0)]);
%            subplot(3,1,1);
            plot(z, Ab(z));
            hold on; 
                plot([z0 z0], [min(Ab(z)) max(Ab(z))], '-'); 
                plot([Abreak.height Abreak.height], [min(Ab(z)) max(Ab(z))], '-');
                    legend('A(t_b,z)','z_b Expected','z_b Observed');
%             hold off;
%            subplot(3,1,2);
%             plot(z,dAb(z));
%             hold on; 
%                 plot([z0 z0], [min(dAb(z)) max(dAb(z))], '-'); 
%                 plot([Abreak.height Abreak.height], [min(dAb(z)) max(dAb(z))], '-');
%             hold off;
%            subplot(3,1,3);
%             plot(z,ddAb(z));
%             hold on; 
%                 plot([z0 z0], [min(ddAb(z)) max(ddAb(z))], '-'); 
%                 plot([Abreak.height Abreak.height], [min(ddAb(z)) max(ddAb(z))], '-');
%             hold off;
            drawnow;
        end
        
        if plot_debug_on
        disp('Plotting...');
        figure(2); clf; % Line plot of moment of breaking
                plot(z,Ab(z),'b');
                    xlabel('z');ylabel('Area (Normalized to Initial Background)');
                    title(['Jump ',num2str(Aback) , ' time ', num2str(Abreak.time)]);
                hold on
                    plot([Abreak.height,Abreak.height], [0,Aback+0.5], 'b-');
                    plot([z0, z0],[0,Aback+0.5], 'r-');
                    legend('A(t_b,z)','z_b Observed','z_b Expected');
                hold off
                drawnow; 
                input('Return to continue.');
        end
