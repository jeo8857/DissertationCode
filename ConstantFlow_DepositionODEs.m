function sys = ConstantFlow_DepositionODEs(t,s,par,Y,res,vol,dp,mech,geom_mode,ratio,Qin,Qex,Vin,Vex)

convf = 100;

% BREATHING DYNAMICS
tt = mod(t,res.T);
if tt>=0 && tt<=res.TI 
    dPmus = -2*res.Pmusmin*tt/(res.TI*res.TE)+res.Pmusmin*res.T/(res.TI*res.TE);
    u = 1;
elseif tt>res.TI && tt<=res.T
    dPmus = res.Pmusmin/(1-exp(-res.TE/res.tau))*(-1/res.tau)*exp(-(tt-res.TI)/res.tau);
    u = 0;
end

% save each column as array
gnum = Y{:,1}; % generation number
n = Y{:,2}; % number of airways in each generation
Ldata = Y{:,3}/convf; % length of airways in each generation converted FROM centimeters TO meters
rdata = Y{:,4}/2/convf; % radius of airways in each generation (diameter cut in half and converted FROM centimeters TO meters)
Ba = Y{:,5}.*pi/180; % branching angle converted FROM degrees TO radians
Ga = Y{:,6}.*pi/180; % gravity angle converted FROM degrees TO radians
V = Y{:,8}./10^6; % volume

% Indices for each compartment
c1_ind = 1; % compartment 1 consists of the trachea (generation 1 in Yeh)
c2_ind = 2:5; % compartment 2 consists of the main bronchus and bronchioles (generations
% 2-5 in Yeh)
c3_ind = 6:16; % compartment 3 consists of the conductive and transitory airways (generations
% % 6-16 in Yeh)
% c4_ind = 17:length(g); % compartment 4 consists of the respiratory airways and alveoli (generations
% % 17-25 in Yeh)
c4_ind = 17:length(gnum)-1; % compartment 4 consists of the respiratory airways and alveoli (generations
% 17-25 in Yeh)
c5_ind = length(gnum);


% SCALING

% Scale to our desired lung volume to initialize
sc = (vol.TLC/sum(V))^(1/3); % scaled to TLC
r0 = rdata*sc;
L0 = Ldata*sc;

if ratio==1
    VDratio = sum(V([c1_ind c2_ind c3_ind]))/sum(V);
%     VRratio = sum(V(c4_ind))/sum(V);
elseif ratio==2
    VDratio = .12;
%     VRratio = .88;
end

VDTLC = VDratio*vol.TLC;
VRTLC = (1-VDratio)*vol.TLC;

% STATE VARIABLES

% Nt = s(1); % mass in trachea
% Dt = s(2); % deposited in trachea
% Nb = s(3); % mass in bronchi
% Db = s(4); % deposited in bronchi
% Nc = s(5); % mass in conducting
% Dc = s(6); % deposited in conducting
% Nra = s(7); % mass in respiratory
% Dra = s(8); % deposited in respiratory

Nt = s(1); % mass in trachea
Dt = s(2); % deposited in trachea
Nb = s(3); % mass in bronchi
Db = s(4); % deposited in bronchi
Nc = s(5); % mass in conducting
Dc = s(6); % deposited in conducting
Nra = s(7); % mass in respiratory
Dra = s(8); % deposited in respiratory
Nalv = s(9);
Dalv = s(10);

if tt>=0 && tt<=res.TI
    Qt = Qin(1);
    Qb = Qin(2);
    Qc = Qin(3);
    Qra = Qin(4);
    Qalv = Qin(4);

    Vt = Vin(1);
    Vb = Vin(2);
    Vc = Vin(3);
    Vra = Vin(4);
    Valv = Vin(5);
    VA = Vin(6);
    
elseif tt>res.TI && tt<=res.T
    Qt = -Qex(1);
    Qb = -Qex(2);
    Qc = -Qex(3);
    Qra = -Qex(4);
    Qalv = -Qex(4);

    Vt = Vex(1);
    Vb = Vex(2);
    Vc = Vex(3);
    Vra = Vex(4);
    Valv = Vex(5);
    VA = Vex(6);
end

if geom_mode==1
    cscale = 1;
    rscale = 1;
elseif geom_mode==2
    % scaled to TLC
    cscale = ((Vt+Vb+Vc)/VDTLC)^(1/2); % updated
%     rscale = (Vra/VRTLC)^(1/3); % updated
    rscale = (VA/VRTLC)^(1/3); % updated
end

r = ones(size(r0));
L = ones(size(L0));

% group generations into compartments
r(c1_ind) = r0(c1_ind)*cscale; 
L(c1_ind) = L0(c1_ind);

r(c2_ind) = r0(c2_ind)*cscale;
L(c2_ind) = L0(c2_ind);

r(c3_ind) = r0(c3_ind)*cscale;
L(c3_ind) = L0(c3_ind);

r(c4_ind) = r0(c4_ind)*rscale;
L(c4_ind) = L0(c4_ind)*rscale;

r(c5_ind) = r0(c5_ind)*rscale;
L(c5_ind) = L0(c5_ind)*rscale;

% calculate parallel resistance based on Pouseille flow - this is the 
% resistance in each generation
mu = par.mu/98.0665; % convert viscosity FROM pascal seconds TO cmh2o seconds (1 cmh2o = 98.0665 Pa)

R = (8/pi)*mu*(L./r.^4)./n; % cmh2o s m^-3 - resistance in identical parallel tubes

Rt = sum(R(c1_ind));
Rb = sum(R(c2_ind));
Rc = sum(R(c3_ind));
% Rra = sum(R(c4_ind));
Rra = sum(R([c4_ind c5_ind]));

% Rt = 1.021*1e3;
% Rb = 0.3369*1e3;
% Rc = 0.3063*1e3;
% Rra = 0.0817*1e3;

% PROBABILITIES

if tt>=0 && tt<=res.TI
    Q(c1_ind) = Qin(1);
    Q(c2_ind) = Qin(2);
    Q(c3_ind) = Qin(3);
    Q(c4_ind) = Qin(4);
    Q(c5_ind) = Qin(4);
    Q = Q(:);

elseif tt>res.TI && tt<=res.T
    Q(c1_ind) = Qex(1);
    Q(c2_ind) = Qex(2);
    Q(c3_ind) = Qex(3);
    Q(c4_ind) = Qex(4);
    Q(c5_ind) = Qex(4);
    Q = Q(:);
end

pI = G_impaction(r,Ba,n,par,dp,Q); % impaction in each compartment
pS = G_sedimentation(n,r,L,Ga,par,dp,Q); % sedimentation in each compartment
pD = G_diffusion(n,r,L,par,dp,Q); % diffusion in each compartment

if mech==1
    pS = 0; pD = 0;
elseif mech==2
    pI = 0; pD = 0;
elseif mech==3
    pI = 0; pS = 0;
elseif mech==5
    pD = 0;
end

p = pI+pS+pD-pI.*pS-pI.*pD-pS.*pD+pI.*pS.*pD;


% pt = sum(cond(p(c1_ind)));
% pb = sum(cond(p(c2_ind)));
% pc = sum(cond(p(c3_ind)));
% pra = sum(cond(p(c4_ind)));

pt = sum(cond(p(c1_ind)));
pb = sum(cond(p(c2_ind)));
pc = sum(cond(p(c3_ind)));
pra = sum(cond(p(c4_ind)));
palv = sum(cond(p(c5_ind)));

%==========================================================================

Cin = par.Cin;

% sys(1) = u*(Qt*Cin-Qt*(Nt/Vt)*pt-Qb*(Nt/Vt)*(1-pt))... %inhale
%     -(1-u)*(Qb*(Nb/Vb)*(1-pb)-Qb*(Nt/Vt)*pt-Qt*(Nt/Vt)*(1-pt)); %exhale
% 
% sys(2) = u*Qt*(Nt/Vt)*pt-(1-u)*Qb*(Nt/Vt)*pt;
% 
% sys(3) = u*(Qb*(Nt/Vt)*(1-pt)-Qb*(Nb/Vb)*pb-Qc*(Nb/Vb)*(1-pb))... %inhale
%     -(1-u)*(Qc*(Nc/Vc)*(1-pc)-Qc*(Nb/Vb)*pb-Qb*(Nb/Vb)*(1-pb)); %exhale
% 
% sys(4) = u*Qb*(Nb/Vb)*pb-(1-u)*Qc*(Nb/Vb)*pb;
% 
% sys(5) = u*(Qc*(Nb/Vb)*(1-pb)-Qc*(Nc/Vc)*pc-Qra*(Nc/Vc)*(1-pc))... %inhale
%     -(1-u)*(Qra*(Nra/Vra)*(1-pra)-Qra*(Nc/Vc)*pc-Qc*(Nc/Vc)*(1-pc)); %exhale
% 
% sys(6) = u*Qc*(Nc/Vc)*pc-(1-u)*Qra*(Nc/Vc)*pc;
% 
% sys(7) = u*(Qra*(Nc/Vc)*(1-pc)-Qra*(Nra/Vra)*pra)... %inhale
%     -(1-u)*(-Qra*(Nra/Vra)*pra-Qra*(Nra/Vra)*(1-pra)); % exhale % NOTE: Am I accounting for what deposited during inhale?
% 
% sys(8) = u*Qra*(Nra/Vra)*pra-(1-u)*Qra*(Nra/Vra)*pra;

sys(1) = u*(Qt*Cin-Qt*(Nt/Vt)*pt-Qb*(Nt/Vt)*(1-pt))... %inhale
    -(1-u)*(Qb*(Nb/Vb)*(1-pb)-Qb*(Nt/Vt)*pt-Qt*(Nt/Vt)*(1-pt)); %exhale

sys(2) = u*Qt*(Nt/Vt)*pt-(1-u)*Qb*(Nt/Vt)*pt;

sys(3) = u*(Qb*(Nt/Vt)*(1-pt)-Qb*(Nb/Vb)*pb-Qc*(Nb/Vb)*(1-pb))... %inhale
    -(1-u)*(Qc*(Nc/Vc)*(1-pc)-Qc*(Nb/Vb)*pb-Qb*(Nb/Vb)*(1-pb)); %exhale

sys(4) = u*Qb*(Nb/Vb)*pb-(1-u)*Qc*(Nb/Vb)*pb;

sys(5) = u*(Qc*(Nb/Vb)*(1-pb)-Qc*(Nc/Vc)*pc-Qra*(Nc/Vc)*(1-pc))... %inhale
    -(1-u)*(Qra*(Nra/Vra)*(1-pra)-Qra*(Nc/Vc)*pc-Qc*(Nc/Vc)*(1-pc)); %exhale

sys(6) = u*Qc*(Nc/Vc)*pc-(1-u)*Qra*(Nc/Vc)*pc;

% lung generations 16-24
sys(7) = u*(Qra*(Nc/Vc)*(1-pc)-Qra*(Nra/Vra)*pra-Qra*(Nra/Vra)*(1-pra))... %inhale
    -(1-u)*(Qra*(Nalv/Valv)*(1-palv)-Qra*(Nra/Vra)*pra-Qra*(Nra/Vra)*(1-pra)); % exhale % NOTE: Am I accounting for what deposited during inhale?

sys(8) = u*Qra*(Nra/Vra)*pra-(1-u)*Qra*(Nra/Vra)*pra;

% lung generation 25 (alveoli)
sys(9) = u*(Qra*(Nra/Vra)*(1-pra)-Qra*(Nalv/Valv)*palv)... %inhale
    -(1-u)*(-Qra*(Nalv/Valv)*palv-Qra*(Nalv/Valv)*(1-palv)); % exhale % NOTE: Am I accounting for what deposited during inhale?

sys(10) = u*Qra*(Nalv/Valv)*palv-(1-u)*Qra*(Nalv/Valv)*palv;

sys = sys';

function c = cond(v)
% calculate conditional probability for each generation
    c(1) = v(1);
    for j = 2:length(v)
        pnot = prod(1-v(1:j-1));
        c(j) = pnot*v(j);
    end
end

end