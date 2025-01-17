function [pD_in,pD_out] = G_diffusion(n,r,L,par,dp,Q)
%Brownian diffusion
r  = r(:);
L = L(:);
% sqrt(2)*pi*dp^2*par.rho*287.05.*2*r
% Kn = par.k*par.temp./(sqrt(2)*pi*dp^2*1032.955718028069)
% pause
% v = abs(Q)./(pi*r.^2);
% v = abs(Q)./A;
v = abs(Q)./(n.*(pi*r.^2));
% v = abs(Q)./(pi*r.^2);
tau = par.rho*dp^2/(18*par.mu); % relaxation time
% Sp = v.*tau; % stopping distance of the particle 
t_res = L./v; % residence time

% cunningham
A1 = 1.257;
A2 = 0.4;
A3 = 0.55;
% lam = 0.01*L.*(1033.2274./(1033.2274+P(2:end)-P(1)))*(par.temp/296.2); % mean free path (Raabe)
% lam = (0.0653*10^-6)*(1033.2274./(1033.2274+P(2:end)-P(1)))*(par.temp/296.2); % mean free path (Raabe)
% lam = 0.01*r.*(1033.2274./(1033.2274+P(2:end)-P(1)))*(par.temp/296.2); % mean free path (Raabe)
% lam = 0.01*r.*(1033.2274./(1033.2274+P))*(par.temp/296.2); % mean free path (Raabe)

% lam = Sp.*(1033.2274./(1033.2274+P))*(par.temp/296.2); % mean free path (Raabe)
% c = 1+2.63*10^-6*(6.23+2.01*exp(-8.32*10^4*dp))/dp;
% MFP = 10^-6;
% c = 1+(2*MFP./dp).*(A1+A2*exp(-A3*dp/MFP));
% c = 1+(2*L/dp).*(A1+A2*exp(-A3*dp./L))
% c = 1+(2*Sp/dp).*(A1+A2*exp(-A3*dp./Sp))
% c = 1+(2*lam/dp).*(A1+A2*exp(-A3*dp./lam));

% P_Pa = (1033.2274+P)*.0980665; % convert pressure FROM cmh2o TO Pascals
% lam = (par.mu/0.4999)*sqrt(pi/(8*par.rho)).*sqrt(1./P_Pa); % mean free path of air in meters
% Kn = lam/dp; % Knudsen number
% Kn = Kn(:);
% c = 1+Kn.*(2.514+0.8*exp(-0.55./Kn)); % coefficients from John Cimbala video - link in meeting notes

lam = 0.0712*10^-6; % @ 37 deg Cel, 100% humidity and 76 cmhg atmospheric pressure (icrp 1994)
c = 1+(lam/dp)*(2.514+0.8*exp(-0.55*(dp/lam)));

Dc = c.*par.k*par.temp/(3*pi*par.mu*dp); % diffusion coefficient
% L(20)^2
% Dc(20)
% Dtime = L(20)^2/Dc(20)
% 2*sqrt(Dc(20)*t_res(20))
% t_res(20)
% Dtime/t_res(20)
% pause
% Dc = par.mu*par.k*par.temp*(1-exp(-t/()))
Sd = 2*sqrt(Dc.*t_res); % diffusion length
% plot(Sd)
% stop
h = Sd.^2./(2*(2*r).^2); % ratio of diffusion length to airway diameter
% h = 2*L.*Dc./(v.*(2*r).^2);

D = ones(length(r),1);

% turbulent
D(1:3) = 2.828*h(1:3).^(1/2).*(1-0.314*h(1:3).^(1/2));

% laminar
D = 1-0.819*exp(-7.315*h)-0.0976*exp(-44.61*h)-0.0325*exp(-114*h)-0.0509*exp(-79.31*h.^(2/3));% yeh 1980

% if any(Q==0)
%     D = 1;
% end

% Pd(1) = D(1);
% Pd(2) = (1-D(1))*D(2);
% Pd(3) = (1-D(1))*(1-D(2))*D(3);
% Pd(4) = (1-D(1))*(1-D(2))*(1-D(3))*D(4);

% figure(2)
% hold on
% plot(dp,h(4),'ro')
% pD_in = Pd.*pV;
pD_in = D;

end

