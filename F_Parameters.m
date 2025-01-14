function par = F_Parameters()

%--------------------------------------------------------------------------
par.Cin = 3*10^8;%.2; % concentration determined by device, mass/volume

% Air properties
mu_scale = 0;
mu = 1.86*10^(-5); % viscosity of air at 25 degrees celcius in P*s (pascal = Newtons/m^2 = kg/(s^2*m))
par.mu = mu*(1+mu_scale); % viscosity of air at 25 degrees celcius in P*s (pascal = Newtons/m^2 = kg/(s^2*m))

par.g = 9.8; % gravitational acceleration m/s^2
par.temp = 300; % normal adult body temperature in kelvins

% LUNG MECHANICS-----------------------------------------------------------

% % compliances (L/cmH2O) 
% par.Ct = .00127; 
% par.Cb = .00238;
% par.Cc = .0131;
% par.Cra = .208;

% compliances (m^3/cmH2O) 
par.Ct = .00127/1e3; 
par.Cb = .00238/1e3;
par.Cc = .0131/1e3;
par.Cra = .208/1e3;

% unstressed volume (L)
% par.Vut = 34.4/1e3;
% par.Vub = 6.63/1e3;
% par.Vuc = 18.7/1e3;
% par.Vura = 740/1e3;
% % par.Vura = 1.263;

% unstressed volume (m^3)
par.Vut = 0.0344/1e3; %34.4/1e3;
par.Vub = 0.00663/1e3; %6.63/1e3;
par.Vuc = 0.0187/1e3; %18.7/1e3;
par.Vura = 1.263/1e3;%0.74/1e3; %740/1e3;
% par.Vura = 0.74/1e3; %740/1e3;
% par.Vura = 2.112/1e3;

% additional parameters
par.CCW = .2445/1e3;% L/cmH2O %*vcf/pcf;
% par.CCW = .2985;%.2445/1e3;% m^3/cmH2O %*vcf/pcf;
par.Pao = 0; % airway opening pressure - zero above atmospheric

% Particle properties
% par.rho = 1; % particle density in kg/L
par.rho = 1000; % particle density in kg/m^3
par.k = 1.38*10^(-23); % Boltzmann constant in Joules = kg*m^2/s^2
end
