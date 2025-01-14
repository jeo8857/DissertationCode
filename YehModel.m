% MORPHOMETRY
sheet_name1 = 'Yeh1980';

Y = readtable('Morphometry','Sheet',sheet_name1,'VariableNamingRule','preserve');

% save each column as array
g = Y{:,1}; % generation number
n = Y{:,2}; % number of airways in each generation
L = Y{:,3}; % length of airways in each generation (cm)
R = Y{:,4}/2; % radius of airways in each generation (cm)
theta = Y{:,5}.*pi/180; % branching angle converted FROM degrees TO radians
phi = Y{:,6}.*pi/180; % gravity angle converted FROM degrees TO radians
V = Y{:,8}./10^6; % volume

dp = [0.01 0.06 0.2 0.6 1 2 3 4 6 10];
%Fi = 1; %inspiratory flow rate (l/min)
%Fe = 1;
%Ni = -0.62+0.475.*log(Da.^2*Fi);
%Ne = -0.62+0.475.*log(Da.^2*Fe);

Vt = 750; % tidal volume (cm^3)
Vn = 50; % volume of NP region (cm^3)
bf = 12; % breathing frequency (bpm)

% Assuming Ni=Ne=N
F = Vt*bf;
N = -0.62+0.475.*log(dp.^2*F);

% plot(log(dp.^2*F))
% stop
plot(dp,N)
stop
TBp = .3;
Pp = .2;
% NP = Ni+Ne*(1-Ni)*(1-(TBp+Pp)*(1-Vn/Vt));
NP = N+N.*(1-N).*(1-(TBp+Pp).*(1-Vn/Vt));

plot(dp,NP)
hold on
plot(dp,N)
stop
v = 1;
g = 980.665; % gravitational acceleration (cm/s^2)
rho_p = 1; % unit particle density (g/cm^3)
temp = 310.15; % body temperature in Kelvins (equivalent to 37 deg celcius)
k = 1;
mu = 1.9*10^-4; %g/cm/s
lam = 0.0712; % @ 37 deg Cel, 100% humidity and 76 cmhg atmospheric pressure (icrp 1994)
C = 1+(lam/dp)*(2.514+0.8*exp(-0.55*(dp/lam)));
D = C.*k*temp/(3*pi*mu*dp);

% laminar diffusion
x = L*D./(2*R.^2*v);
Pd = 1-0.819*exp(-7.315*x)-0.0976*exp(-44.61*x)-0.0325*exp(-114*x)-0.0509*exp(-79.31*x.^(2/3));

% sedimentation
Ps = 1-exp((-4*g*C*rho_p*rp*2*L*cos(phi))/(9*pi*mu*R*v));

% impaction
St = C*rho_p*(dp/2).^2*v./(9*mu*R);
if theta*St<1
    Pi = 1-(2/pi)*acos(theta*St)+(1/pi)*sin(2*acos(theta*St));
else
    Pi = 1;
end

TB = (1-N)*(1-Vn/Vt)*TBp;
P = (1-N)*(1-Vn/Vt)*Pp;