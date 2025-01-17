function [pI_in,pI_out] = G_impaction(r,theta,n,par,dp,Q)

% plot(2*r)
% stop
% v = abs(Q)./A;
v = abs(Q)./(n.*(pi*r.^2)); % velocity of air through each individual airway in each generation, m/s

tau = par.rho*dp^2/(18*par.mu); % relaxation time
Sp = v.*tau; % stopping distance of particle in any individual airway for each generation, m
% if t>1.5
%     figure(1)
%     hold on
%     plot(Q*10^3)
% 
%     figure(2)
%     hold on
%     plot(V*10^3)
% 
%     figure(4)
%     hold on
%     plot(v)
%     
% end

% cunningham
A1 = 1.257;
A2 = 0.4;
A3 = 0.55;

% lam = Sp.*(1033.2274./(1033.2274+P-Ppl))*(par.temp/296.2); % mean free path (Raabe)
% lam = Sp.*(1033.2274./(1033.2274+P))*(par.temp/296.2); % mean free path (Raabe)
% lam = (par.mu/0.4999).*sqrt(pi./(8*par.rho*(1033.2274+P)))

% P_Pa = (1033.2274+(s(:,2)))*.0980665; % convert pressure FROM cmh2o TO Pascals
% lam = 10^6*(par.mu/0.4999)*sqrt(pi/(8*par.rho)).*sqrt(1./P_Pa);
% P_Pa = (1033.2274+P)*.0980665; % convert pressure FROM cmh2o TO Pascals
% lam = (par.mu/0.4999)*sqrt(pi/(8*par.rho)).*sqrt(1./P_Pa); % mean free path of air in meters
% Kn = lam/dp; % Knudsen number
% Kn = Kn(:);
% c = 1+Kn.*(2.514+0.8*exp(-0.55./Kn)); % coefficients from John Cimbala video - link in meeting notes

lam = 0.0712*10^-6; % @ 37 deg Cel, 100% humidity and 76 cmhg atmospheric pressure (icrp 1994)
c = 1+(lam/dp)*(2.514+0.8*exp(-0.55*(dp/lam)));

% lam = (0.0653*10^-6)*(1033.2274./(1033.2274+P(2:end)-P(1)))*(par.temp/296.2); % mean free path (Raabe)
% c = 1+2.63*10^-6*(6.23+2.01*exp(-8.32*10^4*dp))/dp;
% MFP = 10^-6;
% c = 1+(2*MFP./dp).*(A1+A2*exp(-A3*dp/MFP));
% c = 1+(2*L/dp).*(A1+A2*exp(-A3*dp./L))
% c = 1+(2*Sp/dp).*(A1+A2*exp(-A3*dp./Sp))
% c = 1+(2*lam/dp).*(A1+A2*exp(-A3*dp./lam)); % coefficients from Raabe
stk = c.*Sp./(2*r); % stokes number for an individual airway in each generation

% stk = c.*Sp./(L);
I = ones(length(n),1);
% I(find(stk.*theta < 1)) = 1-(2/pi)*acos(theta.*stk)+(1/pi)*sin(2*acos(theta.*stk));

% figure(1)
% hold on
% plot(stk)
% pause
% for i = 1:length(n)
%     if stk(i)*theta(i) < 1
%         I(i) = 1-(2/pi)*acos(theta(i)*stk(i))+(1/pi)*sin(2*acos(theta(i)*stk(i)));
%     else
%         I(i) = 1;
%     end
% end

I = 1-exp(-4*theta.*stk);
% if t>1.5
%     figure(3)
%     hold on
%     plot(I)
% 
%     figure(5)
%     plot(theta.*stk)
%     stop
% end
% Pi(1) = I(1);
% Pi(2) = (1-I(1))*I(2);
% Pi(3) = (1-I(1))*(1-I(2))*I(3);
% Pi(4) = (1-I(1))*(1-I(2))*(1-I(3))*I(4);

% pI_in = Pi.*pV;
pI_in = I; % the probability of impaction in an individual airway in each generation
% pI_in = Pi;
% plot(pI_in)
% hold on
% pause
end

