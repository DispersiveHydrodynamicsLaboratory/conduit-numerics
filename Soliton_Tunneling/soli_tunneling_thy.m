
%% Soliton amplitude relation
%% m is mean of conduit, a is m-total soliton amplitude
csoli = @(a,m) m./a.^2 .* ( (a+m).^2 .* (2*log(1+a./m)-1) + m.^2);

trans_coeff = @(c,m) c .* (c + 2*m) / m;

m_minus = 1;
m_plus  = 1.75;

a = 0.05:0.05:7;
aDSW = [6.2 7.0];
qi = [];
for ii = 1:length(aDSW)
    qi(ii) = find(abs(a-aDSW(ii))<10^(-14));
end

trans_coeff_plus  = trans_coeff(csoli(a,m_plus),m_plus);
trans_coeff_minus = zeros(size(trans_coeff_plus));
trans_soli_size   = zeros(size(trans_coeff_plus));

myfun = @(tcp,am,mm) tcp-trans_coeff(csoli(am,mm),mm);
for ii = 1:numel(trans_coeff_plus)
    am = fsolve(@(am) myfun(trans_coeff_plus(ii),am,m_minus),a(ii));
    trans_soli_size(ii)   = am;
    trans_coeff_minus(ii) = trans_coeff(csoli(am,m_minus),m_minus);
end
aRW = trans_soli_size(qi);

%% Check to make sure coefficients are the same
figure(1); clf;
plot(trans_coeff_plus,trans_coeff_minus,'-',...
     trans_coeff_plus(qi),trans_coeff_minus(qi),'*');
xlabel('Coefficient for (c_+,a_+)');
ylabel('Coefficient for (c_-,a_-)');
%% Plot actual soliton amplitudes
figure(2); clf;
plot(a,trans_soli_size,'-',...
     aDSW,aRW,'*');
xlabel('Amplitude on A_{DSW} background');
ylabel('Amplitude on A_{0} background');

%% Display amplitude pairs to be used
figure(3); clf;
uitable('Data',[aDSW' trans_soli_size(qi)'],...
        'ColumnName',{'DSW Soliton','Rare Soliton'});
    
save('soli_pairs.mat','aDSW','aRW','m_minus','m_plus');


