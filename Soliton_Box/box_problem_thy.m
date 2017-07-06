%% Box problem asymptotics

%% Load previous numerics files
% q = load('box_results_num.mat');
% for ii = 1:length(q.Amax_num)
%     Amax = q.Amax_num(ii);
%     w    = q.widths_num(ii);
%     pks  = q.peak_heights(ii,:);
%% Initial condition properties
% Generate IC
    vstretch = 15;
    ampl = 0.5;
    boxmax = 100;
    z0 = 500;
    f = make_smooth_box(vstretch,ampl,boxmax,z0);
    f = @(z) f(z) + ones(size(z));
    zmin = 0;       % Domain start
    zmax = 1000;%q.zmax;  % Domain end
%     z0   = q.z0;    % Approximate middle of one-hump IC
% f    = @(z) ((Amax-1)/2+1) .* ones(size(z))              + ...
%             (Amax-1)/2 .*+tanh((z-(z0-w/2))/4) .* (z<=z0) + ...
%             (Amax-1)/2 .*-tanh((z-(z0+w/2))/4) .* (z>z0);
% IC integration conditions
bmin = fzero(@(z) f(z)-(1+10^(-3)),[zmin,z0]); % z where box begins
bmax = fzero(@(z) f(z)-(1+10^(-3)),[z0,zmax]); % z where box ends
amin = f(bmin);     % box minimum
amax = ampl;        % box maximum

% k_-, the trailing edge distinguished limit
km = @(phi,z,lambda) sqrt(1/2*(lambda-2./phi(z)+...
        sqrt(lambda./phi(z).*(4+phi(z).*lambda))));
% Derivative of km wrt lambda
dkdl = @(phi,z,lambda) (1 + (2 + lambda.*phi(z))./...
    ( sqrt(lambda.*(lambda + 4./phi(z))).*phi(z)))./...
    (2*sqrt(2).*sqrt( lambda + sqrt(lambda.*...
    (lambda + 4./phi(z))) - 2./phi(z)));
   
% lambda(A), the relationship between total soliton 
%   amplitude and wavelength
lambda = @(A) 4.*(A-1).^4./((-A.^2 + 2.*A.^2.*log(A) + 1).*...
                             (A.^2 + 2.*A.^2.*log(A)-4.*A+3));
% Derivative of lambda wrt a = A-1
dlda = @(a) (4*a.^4*(-4 + 4*(1 + a) + 4*(1 + a)*log(1 + a)))./...
    ((1 - (1 + a).^2 +     2*(1 + a).^2*log(1 + a))*(3 - 4*(1 + a) +...
    (1 + a).^2 +     2*(1 + a).^2*log(1 + a)).^2) +...
    ( 16*a.^4*(1 + a) *log(   1 + a))./((1 - (1 + a).^2 +...
    2*(1 + a).^2*log(1 + a)).^2*(3 -     4*(1 + a) +...
    (1 + a).^2 + 2*(1 + a).^2*log(1 + a))) - ( 16*a.^3)./...
    ((1 - (1 + a).^2 + 2*(1 + a).^2*log(1 + a))*(3 -     4 *(1 + a) +...
    (1 + a).^2 + 2*(1 + a).^2 *log(1 + a)));
% N, the total number of solitons
N = 1/(2*pi) *integral(@(z) km(f,z,1/2),bmin,bmax);

% f(a), the amplitude distribution of solitons
z1fun = @(a) fzero(@(z) f(z)-1./(2*lambda(a+1)),[zmin,zmax/2]);
z2fun = @(a) fzero(@(z) f(z)-1./(2*lambda(a+1)),[zmax/2,zmax]);
a = amin:10^(-3):amax;
fa = zeros(size(a));
for ii = 1:length(a)
    z1 = z1fun(a(ii));
    z2 = z2fun(a(ii));
    fa(ii) = 1/(2*pi)*integral(@(z) dkdl(f,z,lambda(a(ii)+1)).*dlda(a(ii)),z1,z2);
end

figure(1); clf;
z = zmin:10^(-3):zmax;
    subplot(2,1,1);
        plot(z,f(z));
        xlabel('z');
        ylabel('A_0(z)+1');
        title(['Expected Number of Solitons: ',round(num2str(N))]);
    subplot(2,1,2);
        plot(a,fa);
        title('Amplitude Distribution of Solitons');
        xlabel('a');
        ylabel('f(a)');
    disp(N);
%     input('r')
% end