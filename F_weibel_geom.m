function par = weibel_geom()


% All values are in mL,s,and mmHg unless otherwise stated

%--------------------------------------------------------------------------
par.Cin = 200;%.2; % concentration determined by device, mass/volume

% without scaling
cscale = 1; % scaling factor for bronchioles and conducting airways (Schum and Yeh 1980) - multiplies radii
rscale = 1; % scaling factor for respiratory airways (Schum and Yeh 1980) - multiplies radii and lengths

% with scaling factors calculated based on Schum and Yeh
% cscale = 1.4986; % scaling factor for bronchioles and conducting airways (Schum and Yeh 1980) - multiplies radii
% rscale = 1.3095; % scaling factor for respiratory airways (Schum and Yeh 1980) - multiplies radii and lengths

% playing with values
% cscale = 40; % scaling factor for bronchioles and conducting airways (Schum and Yeh 1980) - multiplies radii
% rscale = 40; % scaling factor for respiratory airways (Schum and Yeh 1980) - multiplies radii and lengths

par.RR = 12; % breaths/min
par.IE = 0.6; % dimensionless (shifts the transition between inspiration and expiration in Pmus)
par.T = 60/par.RR; % respiratory period in seconds
par.TE = par.T/(1+par.IE);
par.TI = par.TE*par.IE;
par.tau = par.TE/5; %par.TE/5; 
par.HR = 72; % beats/min
par.per = 60/par.HR; % time to complete 1 beat in seconds
par.Tsys0 = .5; 
par.ksys = .075;
par.TV = .5; % L

% Air properties
par.mu = 1.86*10^(-5); % viscosity of air at 25 degrees celcius in pascals = Newtons/m^2 = kg/(s*m)
par.g = 9.8; % gravitational acceleration m/s^2
par.temp = 310; % normal adult body temperature in kelvins

% LUNG MECHANICS-----------------------------------------------------------

% compliances (L/cmH2O) 
par.Ct = .00127;%*vcf/pcf; 
par.Cb = .00238;%*vcf/pcf;
par.Cc = .0131;%*vcf/pcf;
par.Cra = .2;%*vcf/pcf;
% par.Cra = 2;%*vcf/pcf;
par.Cve = 0.5;


%%% resistances (cmH2O*s*1/L)

% GEOMETRY

% unstressed volume (L)
par.Vut = 34.4/1e3;
par.Vub = 6.63/1e3;
par.Vuc = 18.7/1e3;
par.Vura = 1.263;%*1e3; 

% additional parameters
par.CCW = .2445;% L/cmH2O %*vcf/pcf;
% par.FRC = 2.4*vcf; 
% par.PplEE = -5; % cmH2O %*pcf;
par.Pmusmin = -5; %-5; % cmH2O %*pcf; (this parameter changes the amplitude of Pmus)
par.Pao = 0; % airway opening pressure - zero above atmospheric

% Particle properties
par.rho = 1; % particle density in kg/L
% rho_p = 1000; % particle density in kg/m^3
% par.dp = .000015; % particle diameter in meters
par.k = 1.38*10^(-23); % Boltzmann constant in Joules = kg*m^2/s^2

% Weibel geometry

% number of tubes in generation
N0 = 1;
N1 = 2;
N2 = 4;
N3 = 8;
N4 = 16;
N5 = N4*2;
N6 = N5*2;
N7 = N6*2;
N8 = N7*2;
N9 = N8*2;
N10 = N9*2;
N11 = N10*2;
N12 = N11*2;
N13 = N12*2;
N14 = N13*2;
N15 = N14*2;
N16 = N15*2;
N17 = N16*2;
N18 = N17*2;
N19 = N18*2;
N20 = N19*2;
N21 = N20*2;
N22 = N21*2;
N23 = N22*2;

% vessel radius in cm
r0 = 1.8/2;

% bronchioles
r1 = 1.22/2*cscale;
r2 = 0.83/2*cscale;
r3 = 0.56/2*cscale;
r4 = 0.45/2*cscale;

% conducting airways
r5 = 0.35/2*cscale;
r6 = 0.28/2*cscale;
r7 = 0.23/2*cscale;
r8 = 0.186/2*cscale;
r9 = 0.154/2*cscale;
r10 = 0.13/2*cscale;
r11 = 0.109/2*cscale;
r12 = 0.095/2*cscale;
r13 = 0.082/2*cscale;
r14 = 0.074/2*cscale;
r15 = 0.066/2*cscale;

% respiratory airways
r16 = 0.06/2*rscale;
r17 = 0.054/2*rscale;
r18 = 0.05/2*rscale;
r19 = 0.047/2*rscale;
r20 = 0.045/2*rscale;
r21 = 0.043/2*rscale;
r22 = 0.041/2*rscale;
r23 = 0.041/2*rscale;

% make radius of last generation 5% smaller
% r23 = r23-0.05*r23;

% vessel length in cm
L0 = 12;

% bronchioles and conducting airways
L1 = 4.76;
L2 = 1.9;
L3 = 0.76;
L4 = 1.27;
L5 = 1.07;
L6 = 0.9;
L7 = 0.76;
L8 = 0.64;
L9 = 0.54;
L10 = 0.46;
L11 = 0.39;
L12 = 0.33;
L13 = 0.27;
L14 = 0.23;
L15 = 0.2;

% respiratory airways
L16 = 0.165*rscale;
L17 = 0.141*rscale;
L18 = 0.117*rscale;
L19 = 0.099*rscale;
L20 = 0.083*rscale;
L21 = 0.07*rscale;
L22 = 0.059*rscale;
L23 = 0.05*rscale;

% L16 = 0.165;
% L17 = 0.141;
% L18 = 0.117;
% L19 = 0.099;
% L20 = 0.083;
% L21 = 0.07;
% L22 = 0.059;
% L23 = 0.05;

% parallel resistances in each generation (assuming resistance in each
% parallel vessel is the same)
R0 = (8/pi)*(par.mu/98.0665)*((L0/r0^4)*1000)/N0;

% bronchioles and conducting airways
R1 = (8/pi)*(par.mu/98.0665)*((L1/r1^4)*1000)/N1; 
R2 = (8/pi)*(par.mu/98.0665)*((L2/r2^4)*1000)/N2; 
R3 = (8/pi)*(par.mu/98.0665)*((L3/r3^4)*1000)/N3;
R4 = (8/pi)*(par.mu/98.0665)*((L4/r4^4)*1000)/N4;
R5 = (8/pi)*(par.mu/98.0665)*((L5/r5^4)*1000)/N5;
R6 = (8/pi)*(par.mu/98.0665)*((L6/r6^4)*1000)/N6;
R7 = (8/pi)*(par.mu/98.0665)*((L7/r7^4)*1000)/N7;
R8 = (8/pi)*(par.mu/98.0665)*((L8/r8^4)*1000)/N8;
R9 = (8/pi)*(par.mu/98.0665)*((L9/r9^4)*1000)/N9;
R10 = (8/pi)*(par.mu/98.0665)*((L10/r10^4)*1000)/N10;
R11 = (8/pi)*(par.mu/98.0665)*((L11/r11^4)*1000)/N11;
R12 = (8/pi)*(par.mu/98.0665)*((L12/r12^4)*1000)/N12;
R13 = (8/pi)*(par.mu/98.0665)*((L13/r13^4)*1000)/N13;
R14 = (8/pi)*(par.mu/98.0665)*((L14/r14^4)*1000)/N14;
R15 = (8/pi)*(par.mu/98.0665)*((L15/r15^4)*1000)/N15;
% respiratory airways
R16 = (8/pi)*(par.mu/98.0665)*((L16/r16^4)*1000)/N16;
R17 = (8/pi)*(par.mu/98.0665)*((L17/r17^4)*1000)/N17;
R18 = (8/pi)*(par.mu/98.0665)*((L18/r18^4)*1000)/N18;
R19 = (8/pi)*(par.mu/98.0665)*((L19/r19^4)*1000)/N19;
R20 = (8/pi)*(par.mu/98.0665)*((L20/r20^4)*1000)/N20;
R21 = (8/pi)*(par.mu/98.0665)*((L21/r21^4)*1000)/N21;
R22 = (8/pi)*(par.mu/98.0665)*((L22/r22^4)*1000)/N22;
R23 = (8/pi)*(par.mu/98.0665)*((L23/r23^4)*1000)/N23;




% R1 = (8/pi)*(par.mu/98.0665)*(L1/r1^4)/N1; 
% R2 = (8/pi)*(par.mu/98.0665)*(L2/r2^4)/N2; 
% R3 = (8/pi)*(par.mu/98.0665)*(L3/r3^4)/N3;
% R4 = (8/pi)*(par.mu/98.0665)*(L4/r4^4)/N4;
% R5 = (8/pi)*(par.mu/98.0665)*(L5/r5^4)/N5;
% R6 = (8/pi)*(par.mu/98.0665)*(L6/r6^4)/N6;
% R7 = (8/pi)*(par.mu/98.0665)*(L7/r7^4)/N7;
% R8 = (8/pi)*(par.mu/98.0665)*(L8/r8^4)/N8;
% R9 = (8/pi)*(par.mu/98.0665)*(L9/r9^4)/N9;
% R10 = (8/pi)*(par.mu/98.0665)*(L10/r10^4)/N10;
% R11 = (8/pi)*(par.mu/98.0665)*(L11/r11^4)/N11;
% R12 = (8/pi)*(par.mu/98.0665)*(L12/r12^4)/N12;
% R13 = (8/pi)*(par.mu/98.0665)*(L13/r13^4)/N13;
% R14 = (8/pi)*(par.mu/98.0665)*(L14/r14^4)/N14;
% R15 = (8/pi)*(par.mu/98.0665)*(L15/r15^4)/N15;
% % respiratory airways
% R16 = (8/pi)*(par.mu/98.0665)*(L16/r16^4)/N16;
% R17 = (8/pi)*(par.mu/98.0665)*(L17/r17^4)/N17;
% R18 = (8/pi)*(par.mu/98.0665)*(L18/r18^4)/N18;
% R19 = (8/pi)*(par.mu/98.0665)*(L19/r19^4)/N19;
% R20 = (8/pi)*(par.mu/98.0665)*(L20/r20^4)/N20;
% R21 = (8/pi)*(par.mu/98.0665)*(L21/r21^4)/N21;
% R22 = (8/pi)*(par.mu/98.0665)*(L22/r22^4)/N22;
% R23 = (8/pi)*(par.mu/98.0665)*(L23/r23^4)/N23;

% Lumped resistances (series resistances) - compartments based on Fukui (1972)

% trachea
par.Rt = R0;
par.rt = r0/100; % meters
par.Lt = L0/100; % meters
% par.Rt = par.rt;

% main bronchus and bronchioles (generations 1-4)
par.Rb = R1+R2+R3+R4;
par.rb = ((r1+r2+r3+r4)/4)/100; % meters
par.Lb = ((L1+L2+L3+L4)/4)/100; % meters
% par.Rb = par.rb;

% conductive and transitory airways (generations 5-15)
par.Rc = R5+R6+R7+R8+R9+R10+R11+R12+R13+R14+R15;
par.rc = ((r5+r6+r7+r8+r9+r10+r11+r12+r13+r14+r15)/11)/100; % meters
par.Lc = ((L5+L6+L7+L8+L9+L10+L11+L12+L13+L14+L15)/11)/100; % meters
% par.Rc = par.rc;

% respiratory airways (generations 16-23)
par.Rra = R16+R17+R18+R19+R20+R21+R22+R23;
par.rra = ((r16+r17+r18+r19+r20+r21+r22+r23)/8)/100; % meters
par.Lra = ((L16+L17+L18+L19+L20+L21+L22+L23)/8)/100; % meters
% par.Rra = par.Rra+0.05*par.Rra;

%testing

% geom.t0 = (8/pi)*(par.mu/98.0665)*((geom.L0/geom.r0^4)*1000) 
% geom.t1 = (8/pi)*(par.mu/98.0665)*((geom.L1/geom.r1^4)*1000) 
% geom.t2 = (8/pi)*(par.mu/98.0665)*((geom.L2/geom.r2^4)*1000) 
% geom.t3 = (8/pi)*(par.mu/98.0665)*((geom.L3/geom.r3^4)*1000)
% geom.t4 = (8/pi)*(par.mu/98.0665)*((geom.L4/geom.r4^4)*1000)
% geom.t5 = (8/pi)*(par.mu/98.0665)*((geom.L5/geom.r5^4)*1000)
% geom.t6 = (8/pi)*(par.mu/98.0665)*((geom.L6/geom.r6^4)*1000)
% geom.t7 = (8/pi)*(par.mu/98.0665)*((geom.L7/geom.r7^4)*1000)
% geom.t8 = (8/pi)*(par.mu/98.0665)*((geom.L8/geom.r8^4)*1000)
% geom.t9 = (8/pi)*(par.mu/98.0665)*((geom.L9/geom.r9^4)*1000)
% geom.t10 = (8/pi)*(par.mu/98.0665)*((geom.L10/geom.r10^4)*1000)
% geom.t11 = (8/pi)*(par.mu/98.0665)*((geom.L11/geom.r11^4)*1000)
% geom.t12 = (8/pi)*(par.mu/98.0665)*((geom.L12/geom.r12^4)*1000)
% geom.t13 = (8/pi)*(par.mu/98.0665)*((geom.L13/geom.r13^4)*1000)
% geom.t14 = (8/pi)*(par.mu/98.0665)*((geom.L14/geom.r14^4)*1000)
% geom.t15 = (8/pi)*(par.mu/98.0665)*((geom.L15/geom.r15^4)*1000)
% geom.t16 = (8/pi)*(par.mu/98.0665)*((geom.L16/geom.r16^4)*1000)
% geom.t17 = (8/pi)*(par.mu/98.0665)*((geom.L17/geom.r17^4)*1000)
% geom.t18 = (8/pi)*(par.mu/98.0665)*((geom.L18/geom.r18^4)*1000)
% geom.t19 = (8/pi)*(par.mu/98.0665)*((geom.L19/geom.r19^4)*1000)
% geom.t20 = (8/pi)*(par.mu/98.0665)*((geom.L20/geom.r20^4)*1000)
% geom.t21 = (8/pi)*(par.mu/98.0665)*((geom.L21/geom.r21^4)*1000)
% geom.t22 = (8/pi)*(par.mu/98.0665)*((geom.L22/geom.r22^4)*1000)
% geom.t23 = (8/pi)*(par.mu/98.0665)*((geom.L23/geom.r23^4)*1000)
end

