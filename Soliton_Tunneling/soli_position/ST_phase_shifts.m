% Function for determining soliton phase shift from Whitham theory
function [krat, newc] = ST_phase_shifts(m_minus, m_plus, csoli)

% Theory of phase shifts (m is mean)
V      = @(m)    2*m;                % Nonlinear coefficient
q      = @(cs,m) cs*(cs+2*m)./(4*m);  % Riemann invariant
cs     = @(q,m) -m+sqrt(m.*(m+4*q)); % speed of soliton
q0     = q(csoli, m_minus);

newc   = fzero(@(cs) q(cs,m_plus)-q0, csoli+1);

% Derivative of soliton speed wrt mean m
dcsdm  = @(q,m) -1 + (2*m + 4*q0)./(2*sqrt(m.^2 + 4*m*q0));

% Integrand for determining phase shift
I1 = @(q,m) dcsdm(q,m)./(cs(q,m)-V(m));

% Integrate and take exponential
logkrat  = integral(@(m) I1(q0,m),m_minus,m_plus);
krat = exp(logkrat);






