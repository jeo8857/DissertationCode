function pS_in = sedimentation(n,r,L,phi,par,dp,Q)
% sedimentation
% r
r = r(:);
L = L(:);
% v = abs(Q)./A;
v = abs(Q)./(n.*(pi*r.^2));

% Q(1)
% figure(1)
% hold on
% plot(1000*r,1000*v)
% pause

% v = abs(Q)./(pi*n.*r.^2)
tau = par.rho*dp^2/(18*par.mu); % relaxation time
Sp = v.*tau; % stopping distance of the particle 
vs = par.g*tau; % settling/terminal velocity
t_res = L./v; % residence time in an individual airway for each generation
% figure(2)
% plot(t_res)
% stop
% t_rem = t_res-tau;

e = (3*vs.*L)./(8*v.*r);
% e.^(2/3)
% 1-e.^(2/3)

% pause
% cunningham
A1 = 1.257;
A2 = 0.4;
A3 = 0.55;
% lam = Sp.*(1033.2274./(1033.2274+P(2:end)-P(1)))*(par.temp/296.2); % mean free path (Raabe)
% lam = Sp.*(1033.2274./(1033.2274+P))*(par.temp/296.2); % mean free path (Raabe)

% lam = (0.0653*10^-6)*(1033.2274./(1033.2274+P(2:end)-P(1)))*(par.temp/296.2); % mean free path (Raabe)
% c = 1+2.63*10^-6*(6.23+2.01*exp(-8.32*10^4*dp))/dp;
% MFP = 10^-6;
% c = 1+(2*MFP./dp).*(A1+A2*exp(-A3*dp/MFP));
% c = 1+(2*L/dp).*(A1+A2*exp(-A3*dp./L))
% c = 1+(2*(vs.*(t_res))/dp).*(A1+A2*exp(-A3*dp./(vs.*t_res)));
% c = 1+(2*lam/dp).*(A1+A2*exp(-A3*dp./lam));

% P_Pa = (1033.2274+P)*.0980665; % convert pressure FROM cmh2o TO Pascals
% lam = (par.mu/0.4999)*sqrt(pi/(8*par.rho)).*sqrt(1./P_Pa); % mean free path of air in meters
% Kn = lam/dp; % Knudsen number
% Kn = Kn(:);
% c = 1+Kn.*(2.514+0.8*exp(-0.55./Kn)); % coefficients from John Cimbala video - link in meeting notes

lam = 0.0712*10^-6; % @ 37 deg Cel, 100% humidity and 76 cmhg atmospheric pressure (icrp 1994)
c = 1+(lam/dp)*(2.514+0.8*exp(-0.55*(dp/lam)));

S = 1-exp(-(2/pi)*cos(phi).*c.*vs.*t_res./r);

% if any(Q==0)
%     t_res
%     size(S)
%     S = ones(length(v));
% end

% plot(S)
% stop
% S = (2/pi)*(2*e.*sqrt(1-e.^(2/3))-e.^(1/3).*sqrt(1-e.^(2/3))+asin(e.^(1/3)))

% Ps(1) = S(1);
% Ps(2) = (1-S(1))*S(2);
% Ps(3) = (1-S(1))*(1-S(2))*S(3);
% Ps(4) = (1-S(1))*(1-S(2))*(1-S(3))*S(4);

% pS_in = Ps.*pV;
pS_in = S;
% plot(pS_in)
% hold on
% pause
% clc
% pS_in

% pause

% size(S)
% size(pS_in)
% pause
end

