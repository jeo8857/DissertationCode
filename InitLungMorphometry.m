function morph = InitLungMorphometry(data_tog,FRC)
% NOTES:
% unit of length = meter
par = F_Parameters();

if data_tog==1
    sheet_name = 'Yeh1980';
    display('Using Yeh data')
end

if data_tog==2
    sheet_name = 'Weibel1963';
    display('Using Weibel data')
end

Y = readtable('Morphometry','Sheet',sheet_name,'VariableNamingRule','preserve');

% cscale = 1; % scaling factor for bronchioles and conducting airways (Schum and Yeh 1980) - multiplies radii
% rscale = 1; % scaling factor for respiratory airways (Schum and Yeh 1980) - multiplies radii and lengths

% cscale = 1.4986; % scaling factor for bronchioles and conducting airways (Schum and Yeh 1980) - multiplies radii
% rscale = 1.3095; % scaling factor for respiratory airways (Schum and Yeh 1980) - multiplies radii and lengths

% cscale = 2.5; % scaling factor for bronchioles and conducting airways (Schum and Yeh 1980) - multiplies radii
% rscale = 2.5; % scaling factor for respiratory airways (Schum and Yeh 1980) - multiplies radii and lengths

cscale = sqrt(par.VDTLC/(.3*TV));
rscale = (par.VRTLC/(.7*TV))^(1/3);

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
    V = A.*L; %Y{:,8}; % volume
    CV = Y{:,9}; % cumulative volume
end

if data_tog==2 % Weibel
    g = Y{:,1}; % generation number
    n = Y{:,2}; % number of airways in each generation
    L = Y{:,3}/100; % length of airways in each generation
    r = Y{:,4}/2/100; % radius of airways in each generation
    A = n.*(pi*r.^2);
end

% group generations into compartments
% compartment 1 consists of the trachea (generation 1 in Yeh)
c1_ind = 1;

% compartment 2 consists of the main bronchus and bronchioles (generations
% 2-5 in Yeh)
c2_ind = 2:5;
r(c2_ind) = r(c2_ind)*cscale;

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

% R = (8/pi)*mu*((L./r.^4)*0.001)./n; % cmh2o s L^-1
R = (8/pi)*mu*(L./r.^4)./n; % cmh2o s m^-3

% trachea
morph.Rt = R(c1_ind); 
morph.rt = r(c1_ind); 
morph.Lt = L(c1_ind); 
morph.Bat = Ba(c1_ind);
% morph.Rt = (8/pi)*mu*(morph.Lt/morph.rt); % cmh2o s m^-3

% main bronchus and bronchioles
morph.Rb = sum(R(c2_ind));  
morph.rb = (sum(r(c2_ind))/length(c2_ind)); 
morph.Lb = (sum(L(c2_ind))/length(c2_ind)); 
morph.Bab = sum(Ba(c2_ind))/length(c2_ind);
% morph.Rb = (8/pi)*mu*(morph.Lb/morph.rb);

% conductive and transitory airways
morph.Rc = sum(R(c3_ind));  
morph.rc = (sum(r(c3_ind))/length(c3_ind)); 
morph.Lc = (sum(L(c3_ind))/length(c3_ind)); 
morph.Bac = sum(Ba(c3_ind))/length(c3_ind);
% morph.Rc = (8/pi)*mu*(morph.Lc/morph.rc);

% respiratory airways and alveoli
morph.Rra = sum(R(c4_ind));
morph.rra = (sum(r(c4_ind))/length(c4_ind));
morph.Lra = (sum(L(c4_ind))/length(c4_ind)); 
morph.Bara = sum(Ba(c4_ind))/length(c4_ind);
morph.Rra = (8/pi)*mu*(morph.Lra/morph.rra);

end
