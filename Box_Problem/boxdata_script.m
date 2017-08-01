load('boxdata.mat','AMAXES','AMPLS','NS','WS');

for ii = 1:16 %Trial number, current maximum is 16
    Amax = AMAXES(ii);
    w    = WS(ii);
    N    = NS(ii);
    ampl = AMPLS.(['run',num2str(ii,'%05d')]);
    disp(['A box with amplitude ',num2str(Amax),...
          ' and width ',num2str(w),...
          ' has ',num2str(N),' solitons. ',sprintf('\n'),...
          'Amplitudes of solitons follow.']);
      pause(2);
	disp(ampl');
end