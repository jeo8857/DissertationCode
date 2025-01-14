% RUN DRIVER BEFORE RUNNING THIS SCRIPT
% close all;clear;clc
A_Driver

% compare relaxation time to residence time

% If the residence time is greater than the settling time, set the particle
% velocity to the settling velocity

% If the residence time is less than the settling time, then the particle
% does not have enough time to reach its settling velocity, so its velocity
% remains equal to the flow velocity

dp = 10*10^-6;

morph = F_LungMorphometry(data_tog,vol,res,dp);
par = F_Parameters();
L = morph.L;
r = morph.r;
theta = morph.Ba;
phi = morph.Ga;
n = morph.n;
V = morph.V;
pV = V./sum(V);

c1_ind = 1;
c2_ind = 2:5;
c3_ind = 6:16;
c4_ind = 17:25;


Qc = vol.TV/res.TI; % characteristic volumetric flow in m^3/s
v = Qc./(pi*r.^2);

L./v

stop
% cunningham
A1 = 1.257;
A2 = 0.4;
A3 = 0.55;
c = 1+2.63*10^-6*(6.23+2.01*exp(-8.32*10^4*(dp*10^6)))/(dp*10^6);

% impaction
tau = c.*par.rho*dp^2/(18*par.mu); % relaxation time
Sp = v.*tau; % stopping distance of the particle during inhale
stk = Sp./(2*r);
Pi = zeros(25,1);
for i = 1:25
    if stk(i)*theta(i) < 1
        Pi(i) = 1-(2/pi)*acos(theta(i)*stk(i))+(1/pi)*sin(2*acos(theta(i)*stk(i)));
    else
        Pi(i) = 1;
    end
end

% Pi = pV.*Pi;
Pi = G_impaction(r,theta,V,par,vol,res,dp);

% sedimentation
vs = par.g*tau; % settling/terminal velocity
t_res = L./v;
Ps = 1-exp(-(2/pi)*cos(phi).*vs.*t_res./r);
% Ps = Ps.*pV;
Ps = G_sedimentation(r,L,phi,V,par,vol,res,dp);

% diffusion
D = par.k*par.temp*c/(3*pi*par.mu*dp); % diffusion coefficient
Sd = 2*sqrt(D.*t_res); % diffusion length
h = Sd.^2./(2*r); % ratio of diffusion length to airway diameter
Pd = 1-0.819*exp(-7.315*h)-0.0976*exp(-44.61*h)-0.0325*exp(-114*h);
% Pd = Pd.*pV;
Pd = G_diffusion(r,L,V,par,vol,res,dp);

P = Pi+Ps+Pd-Pi.*Ps-Pi.*Pd-Ps.*Pd+Pi.*Ps.*Pd;
fr = cond(P,1-P);

% figure(1)
% hold on
% plot([c1_ind c1_ind(end)+1],L([c1_ind c1_ind(end)+1])./v([c1_ind c1_ind(end)+1]),'-b.', ...
%     [c2_ind c2_ind(end)+1],L([c2_ind c2_ind(end)+1])./v([c2_ind c2_ind(end)+1]),'-r.', ...
%     [c3_ind c3_ind(end)+1],L([c3_ind c3_ind(end)+1])./v([c3_ind c3_ind(end)+1]),'-k.', ...
%     c4_ind,L(c4_ind)./v(c4_ind),'-m.')
% yline(tau,'r')
% set(gca,'YScale','log')
% title('time')

% figure(2)
% hold on
% plot(1:25,v,'-b.')
% yline(vs,'r')
% set(gca,'YScale','log')
% title('velocity')

% figure(3)
% hold on
% plot(1:25,Ps,'-r.')
% set(gca,'YScale','log')
% title('probability')
% 
% figure(4)
% hold on
% plot(1:25,Pi,'-r.')
% set(gca,'YScale','log')
% 
% figure(5)
% hold on
% plot(1:25,Pd,'-r.')
% set(gca,'YScale','log')

% figure(6)
% hold on
% % plot(1:25,Pi,'b.')
% % plot(1:25,Ps,'k.')
% % plot(1:25,Pd,'r.')
% plot(1:25,P,'-m.')
% plot(1:25,fr,'k.')
% set(gca,'YScale','log')

function C = cond(v,v_compl)
% calculate conditional probability for each generation
    C(1) = v(1);
    for i = 2:length(v)
        pnot = prod(v_compl(1:i-1));
%         C(i) = n(i)*pnot*v(i);
        C(i) = pnot*v(i);
    end
end

function c = comb(v)
    for k = 1:length(v)
        comb = nchoosek(v,k);
        pcomb = prod(comb,2);
        scomb(k) = (-1)^(k+1)*sum(pcomb);
    end
    c = sum(scomb);
end