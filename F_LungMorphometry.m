function morph = F_LungMorphometry(data_tog,volumes,respiratory_settings)
% NOTES:
% unit of length = meter
par = F_Parameters();
res = respiratory_settings;
vol = volumes;

if data_tog==1
    sheet_name1 = 'Yeh1980';
    display('Using Yeh data')
end

if data_tog==2
    sheet_name1 = 'Weibel1963';
    sheet_name2 = 'Yeh1980';
    display('Using Weibel data')
end

Y = readtable('Morphometry','Sheet',sheet_name1,'VariableNamingRule','preserve');

cscale = 1;%.6; % scaling factor for bronchioles and conducting airways (Schum and Yeh 1980) - multiplies radii
rscale = 1;%.1; % scaling factor for respiratory airways (Schum and Yeh 1980) - multiplies radii and lengths

% DETERMINE PERCENT CHANGE IN OUTCOME WITH CHANGE IN SCALE FACTOR
% cscale = 0.7407; % scaling factor for bronchioles and conducting airways (Schum and Yeh 1980) - multiplies radii
% rscale = 0.7637; % scaling factor for respiratory airways (Schum and Yeh 1980) - multiplies radii and lengths
% 
% cscale = 0.5; % scaling factor for bronchioles and conducting airways (Schum and Yeh 1980) - multiplies radii
% rscale = 0.5; % scaling factor for respiratory airways (Schum and Yeh 1980) - multiplies radii and lengths
% 
% cscale = 0.3; % scaling factor for bronchioles and conducting airways (Schum and Yeh 1980) - multiplies radii
% rscale = 0.3; % scaling factor for respiratory airways (Schum and Yeh 1980) - multiplies radii and lengths

% cscale = 1.35; % scaling factor for bronchioles and conducting airways (Schum and Yeh 1980) - multiplies radii
% rscale = 1.25; % scaling factor for respiratory airways (Schum and Yeh 1980) - multiplies radii and lengths

% cscale = 10; % scaling factor for bronchioles and conducting airways (Schum and Yeh 1980) - multiplies radii
% rscale = 10; % scaling

% Dead space volume at current lung volume: VD = Vt+Vb+Vc
% Dead space volume at TLC: VDTLC = 0.3TLC (find a source for this)
% Respiratory volume at current lung volume: VR= Vra
% Respiratory volume at TLC: VRTLC = 0.7TLC (find a source for this)
% cscale = sqrt(vol.VDTLC/(.3*(vol.FRC+vol.TV))); % old
% cscale = sqrt(vol.VDTLC/(Vt+Vb+Vc)); % updated
% % rscale = (vol.VRTLC/(.7*(vol.FRC+vol.TV)))^(1/3); % old
% rscale = (vol.VRTLC/Vra)^(1/3); % updated
convf = 100;

% save each column as array
if data_tog == 1 % Yeh
    g = Y{:,1}; % generation number
    n = Y{:,2}; % number of airways in each generation
    L = Y{:,3}/convf; % length of airways in each generation converted FROM centimeters TO meters
    r = Y{:,4}/2/convf; % radius of airways in each generation (diameter cut in half and converted FROM centimeters TO meters)
    Ba = Y{:,5}.*pi/180; % branching angle converted FROM degrees TO radians
    Ga = Y{:,6}.*pi/180; % gravity angle converted FROM degrees TO radians
    A = n.*(pi*r.^2); %Y{:,7}; % cross sectional area calculated using area equation since I already converted the radii
    V = Y{:,8}./10^6; % volume
    CV = Y{:,9}./10^6; % cumulative volume
end

% figure
% plot(g,2*r./L,'-o')
% pause
if data_tog==2 % Weibel
    g = Y{:,1}; % generation number
    n = Y{:,2}; % number of airways in each generation
    L = Y{:,3}/convf; % length of airways in each generation
    r = Y{:,4}/2/convf; % radius of airways in each generation
    A = n.*(pi*r.^2);

    Y = readtable('Morphometry','Sheet',sheet_name2,'VariableNamingRule','preserve');

    Ba = Y{:,5}.*pi/180; % branching angle converted FROM degrees TO radians
    Ga = Y{:,6}.*pi/180; % gravity angle converted FROM degrees TO radians
    
end

% group generations into compartments
% compartment 1 consists of the trachea (generation 1 in Yeh)
c1_ind = 1;

% compartment 2 consists of the main bronchus and bronchioles (generations
% 2-5 in Yeh)
c2_ind = 2:5;
r(c2_ind) = r(c2_ind)*cscale; % SHOULD THIS BE SCALED?

% compartment 3 consists of the conductive and transitory airways (generations
% 6-16 in Yeh)
c3_ind = 6:16;
r(c3_ind) = r(c3_ind)*cscale;

% compartment 4 consists of the respiratory airways and alveoli (generations
% 17-25 in Yeh)
c4_ind = 17:length(g);
r(c4_ind) = r(c4_ind)*rscale;
L(c4_ind) = L(c4_ind)*rscale;

% calculate parallel resistance based on Pouseille flow - this is the 
% resistance in each generation
mu = par.mu/98.0665; % convert viscosity FROM pascal seconds TO cmh2o seconds (1 cmh2o = 98.0665 Pa)

R = (8/pi)*mu*(L./r.^4)./n; % cmh2o s m^-3 - resistance in identical parallel tubes
% R = (8/pi)*mu*(L./r.^4)./(2.^(g-1)); % cmh2o s m^-3 - resistance in identical parallel tubes
% R = (8/pi)*mu*(L./r.^4); % cmh2o s m^-3 - resistance in single tube

morph.c1_ind = c1_ind;
morph.c2_ind = c2_ind;
morph.c3_ind = c3_ind;
morph.c4_ind = c4_ind;
morph.R = R;
morph.r = r;
morph.L = L;
morph.Ba = Ba;
morph.Ga = Ga;
morph.n = n;
% morph.Rt = R(c1_ind); 
% avg_r(1) = r(c1_ind);
% avg_L(1) = L(c1_ind);
% avg_Ba(1) = Ba(c1_ind);
% avg_Ga(1) = Ga(c1_ind);
% avg_A(1) = A(c1_ind);

% main bronchus and bronchioles
% morph.Rb = sum(R(c2_ind)); 
% avg_r(2) = mean(r(c2_ind));
% avg_L(2) = mean(L(c2_ind));
% avg_Ba(2) = mean(Ba(c2_ind));
% avg_Ga(2) = mean(Ga(c2_ind));
% avg_A(2) = mean(A(c2_ind));

% conductive and transitory airways
% morph.Rc = sum(R(c3_ind)); 
% avg_r(3) = mean(r(c3_ind));
% avg_L(3) = mean(L(c3_ind));
% avg_Ba(3) = mean(Ba(c3_ind));
% avg_Ga(3) = mean(Ga(c3_ind));
% avg_A(3) = mean(A(c3_ind));

% respiratory airways and alveoli
% morph.Rra = sum(R(c4_ind));
% avg_r(4) = mean(r(c4_ind));
% avg_L(4) = mean(L(c4_ind));
% avg_Ba(4) = mean(Ba(c4_ind));
% avg_Ga(4) = mean(Ga(c4_ind));
% avg_A(4) = mean(A(c4_ind));

% morph.r = avg_r;
% morph.L = avg_L;
% morph.Ba = avg_Ba;
% morph.Ga = avg_Ga;
% morph.A = avg_A;
% morph.n = n;
% morph.V = V;
% morph.c1_ind = c1_ind;
% morph.c2_ind = c2_ind;
% morph.c3_ind = c3_ind;
% morph.c4_ind = c4_ind;
end
