load('boxdata.mat','AMAXES','AMPLS','NS','WS');

for ii = 1:16 %Trial number, current maximum is 16
    Amax = AMAXES(ii);
    w(ii)    = WS(ii);
    N(ii)    = NS(ii);
    %ampl = AMPLS.(['run',num2str(ii,'%05d')]);
    n(ii) = amp_dist(Amax,w(ii));
    %disp(['A box with amplitude ',num2str(Amax(ii)),...
          %' and width ',num2str(w),...
          %' has ',num2str(N(ii)),' solitons. ',sprintf('\n')]);...
          %'Amplitudes of solitons follow.']);
      %pause(2);
	%disp(ampl');
    %disp(['Theory says this box should have ', num2str(n(ii)), ' solitons.']);
    
end

figure(1);
plot(w(1:4), N(1:4),'o',w(1:4), n(1:4), 'x','linewidth',2);
title('N Values for Different Widths at Amax = 1');
set(gcf,'color','white');
set(gca,'fontsize',13);
xlabel('Box Width');
ylabel('N Value');
legend('Experiments','Theory');
hold off

figure(2);
plot(w(5:8), N(5:8),'o',w(5:8), n(5:8), 'x','linewidth',2);
title('N Values for Different Widths at Amax = 2');
set(gcf,'color','white');
set(gca,'fontsize',13);
xlabel('Box Width');
ylabel('N Value');
legend('Experiments','Theory');
hold off

figure(3);
plot(w(9:12), N(9:12),w(9:12),'o', n(9:12), 'x','linewidth',2);
title('N Values for Different Widths at Amax = 3');
set(gcf,'color','white');
set(gca,'fontsize',13);
xlabel('Box Width');
ylabel('N Value');
legend('Experiments','Theory');
hold off

figure(4);
plot(w(13:16), N(13:16),'o',w(13:16), n(13:16), 'x','linewidth',2);
title('N Values for Different Widths at Amax = 4');
set(gcf,'color','white');
set(gca,'fontsize',13);
xlabel('Box Width');
ylabel('N Value');
legend('Experiments','Theory');
hold off
