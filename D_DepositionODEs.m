function sys = D_DepositionODEs(t,s,par,Y,res,vol,dp,mech,geom_mode,output,ratio,Qin,Qex,Vin,Vex)

convf = 100;

% BREATHING DYNAMICS
tt = mod(t,res.T);
if tt>=0 && tt<=res.TI 
    dPmus = -2*res.Pmusmin*tt/(res.TI*res.TE)+res.Pmusmin*res.T/(res.TI*res.TE);
    u = 1;
elseif tt>res.TI && tt<=res.T
    dPmus = res.Pmusmin/(1-exp(-res.TE/res.tau))*(-1/res.tau)*exp(-(tt-res.TI)/res.tau);
%     dPmus = -2*res.Pmusmin*tt/(res.TI*res.TE)+res.Pmusmin*res.T/(res.TI*res.TE);
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
% 6-16 in Yeh)
% c4_ind = 17:length(gnum); % compartment 4 consists of the respiratory airways and alveoli (generations
% 17-25 in Yeh)
c4_ind = 17:length(gnum)-1; % compartment 4 consists of the respiratory airways and alveoli (generations
% 17-25 in Yeh)
c5_ind = length(gnum);

% SCALING

% Scale to our desired lung volume to initialize
sc = (vol.TLC/sum(V))^(1/3); % scaled to TLC
% sc = (vol.FRC/sum(V))^(1/3); % scaled to FRC
% sc = ((vol.FRC+0.5*vol.TV)/sum(V))^(1/3);
% sc = ((vol.FRC)/sum(V))^(1/3);
r0 = rdata*sc;
L0 = Ldata*sc;
% figure(10)
% hold on
% plot(r0)
% pause
% V0 = pi*r0(1:end-1).^2.*L0(1:end-1).*n(1:end-1)
% sum(V0)/vol.TLC
% sum(V(1:end-1))/sum(V)

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

Ppl = s(1); % intrapleural pressure
Pt = s(2); % pressure in larynx
Pb = s(3); % pressure in trachea
Pc = s(4); % pressure in bronchai
Pra = s(5); % pressure in alveoli (lung pressure)

% Nt = s(6); % mass in trachea
% Dt = s(7); % deposited in trachea
% Nb = s(8); % mass in bronchi
% Db = s(9); % deposited in bronchi
% Nc = s(10); % mass in conducting
% Dc = s(11); % deposited in conducting
% Nra = s(12); % mass in respiratory
% Dra = s(13); % deposited in respiratory

Nt = s(6); % mass in trachea
Dt = s(7); % deposited in trachea
Nb = s(8); % mass in bronchi
Db = s(9); % deposited in bronchi
Nc = s(10); % mass in conducting
Dc = s(11); % deposited in conducting
Nra = s(12); % mass in respiratory
Dra = s(13); % deposited in respiratory
Nalv = s(14);
Dalv = s(15);

% LUNG VOLUMES
Vt = par.Ct*Pt+par.Vut; % volume of air in trachea
Vb = par.Cb*(Pb-Ppl)+par.Vub; % volume of air in bronchea
Vc = par.Cc*(Pc-Ppl)+par.Vuc; % volume of air in conducting airways
% Vra = par.Cra*(Pra-Ppl)+par.Vura; % volume of air in respiratory airways
VA = par.Cra*(Pra-Ppl)+par.Vura; % volume of air in respiratory airways
% Vtotal = Vt+Vb+Vc+Vra;
Vtotal = Vt+Vb+Vc+VA;
    
% VDFRC = Vt(1)+Vb(1)+Vc(1); % deadspace volume at FRC
% VRFRC = Vra(1); % respiratory volume at FRC

if geom_mode==1
    cscale = 1;
    rscale = 1;
elseif geom_mode==2
    % scaled to TLC
    cscale = ((Vt+Vb+Vc)/VDTLC)^(1/2); % updated
%     rscale = (Vra/VRTLC)^(1/3); % updated
    rscale = (VA/VRTLC)^(1/3); % updated
elseif geom_mode==3
    cscale = 1;
%     rscale = (Vra./VRFRC).^(1/3);
    rscale = (VA./VRFRC).^(1/3);
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

Vra = sum(pi*r(c4_ind).^2.*L(c4_ind).*n(c4_ind));
Valv = VA-Vra;

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

% Q(c1_ind) = Qt;
% Q(c2_ind) = Qb;
% Q(c3_ind) = Qc;
% Q(c4_ind) = Qra;
% Q = Q(:);

if tt>=0 && tt<=res.TI
    Q(c1_ind) = Qin(1);
    Q(c2_ind) = Qin(2);
    Q(c3_ind) = Qin(3);
    Q(c4_ind) = Qin(4);
    Q(c5_ind) = Qin(4);
    Q = Q(:);
    
    cscale = ((Vin(1)+Vin(2)+Vin(3))/VDTLC)^(1/2); % updated
    rscale = (Vin(4)/VRTLC)^(1/3); % updated
%     cscale
%     stop
elseif tt>res.TI && tt<=res.T
    Q(c1_ind) = Qex(1);
    Q(c2_ind) = Qex(2);
    Q(c3_ind) = Qex(3);
    Q(c4_ind) = Qex(4);
    Q(c5_ind) = Qex(4);
    Q = Q(:);

    cscale = ((Vex(1)+Vex(2)+Vex(3))/VDTLC)^(1/2); % updated
    rscale = (Vex(4)/VRTLC)^(1/3); % updated
    
end

rc = ones(size(r0));
Lc = ones(size(L0));

% group generations into compartments
rc(c1_ind) = r0(c1_ind)*cscale; 
Lc(c1_ind) = L0(c1_ind);

rc(c2_ind) = r0(c2_ind)*cscale;
Lc(c2_ind) = L0(c2_ind);

rc(c3_ind) = r0(c3_ind)*cscale;
Lc(c3_ind) = L0(c3_ind);

rc(c4_ind) = r0(c4_ind)*rscale;
Lc(c4_ind) = L0(c4_ind)*rscale;

rc(c5_ind) = r0(c5_ind)*rscale;
Lc(c5_ind) = L0(c5_ind)*rscale;

pI = G_impaction(rc,Ba,n,par,dp,Q); % impaction in each compartment
pS = G_sedimentation(n,rc,Lc,Ga,par,dp,Q); % sedimentation in each compartment
pD = G_diffusion(n,rc,Lc,par,dp,Q); % diffusion in each compartment

% pI = G_impaction(r,Ba,n,par,dp,Q); % impaction in each compartment
% pS = G_sedimentation(n,r,L,Ga,par,dp,Q); % sedimentation in each compartment
% pD = G_diffusion(n,r,L,par,dp,Q); % diffusion in each compartment

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


pt = sum(cond(p(c1_ind)));
pb = sum(cond(p(c2_ind)));
pc = sum(cond(p(c3_ind)));
pra = sum(cond(p(c4_ind)));
% stop
palv = sum(cond(p(c5_ind)));

%==========================================================================

if output==4
%     pt=0;
%     pb=0;
%     pc=0;
%     pra=0;
    if t<=res.T
        Cin = par.Cin;
    elseif t>res.T
        Cin = 0;
    end
else
    Cin = par.Cin;
end

% VOLUMETRIC FLOWS

Qt = (par.Pao-Pt)/Rt;
Qb = (Pt-Pb)/Rb; 
Qc = (Pb-Pc)/Rc; 
Qra = (Pc-Pra)/Rra;

% PRESSURE ODEs

sys(1) = 1/par.CCW*Qb+dPmus; % Ppl - intrapleural pressure 
sys(2) = 1/par.Ct*(Qt-Qb); % Pt - pressure in trachea 
sys(3) = 1/par.Cb*(Qb-Qc)+sys(1); % Pb - pressure in bronchea
sys(4) = 1/par.Cc*(Qc-Qra)+sys(1); % Pc - bronchial conducting airways
sys(5) = 1/par.Cra*(Qra)+sys(1); % Pra - pressure in respiratory airways 

% particles enter on every breath
% Cin = par.Cin;

% Deposition during inhale and exhale

% particles entering at time t are exempt from deposition
sys(6) = u*(Qt*Cin-Qt*(Nt/Vt)*pt-Qb*(Nt/Vt)*(1-pt))... %inhale
    -(1-u)*(Qb*(Nb/Vb)*(1-pb)-Qb*(Nt/Vt)*pt-Qt*(Nt/Vt)*(1-pt)); %exhale

sys(7) = u*Qt*(Nt/Vt)*pt-(1-u)*Qb*(Nt/Vt)*pt;

sys(8) = u*(Qb*(Nt/Vt)*(1-pt)-Qb*(Nb/Vb)*pb-Qc*(Nb/Vb)*(1-pb))... %inhale
    -(1-u)*(Qc*(Nc/Vc)*(1-pc)-Qc*(Nb/Vb)*pb-Qb*(Nb/Vb)*(1-pb)); %exhale

sys(9) = u*Qb*(Nb/Vb)*pb-(1-u)*Qc*(Nb/Vb)*pb;

sys(10) = u*(Qc*(Nb/Vb)*(1-pb)-Qc*(Nc/Vc)*pc-Qra*(Nc/Vc)*(1-pc))... %inhale
    -(1-u)*(Qra*(Nra/Vra)*(1-pra)-Qra*(Nc/Vc)*pc-Qc*(Nc/Vc)*(1-pc)); %exhale

sys(11) = u*Qc*(Nc/Vc)*pc-(1-u)*Qra*(Nc/Vc)*pc;

% lung generations 16-24
sys(12) = u*(Qra*(Nc/Vc)*(1-pc)-Qra*(Nra/Vra)*pra-Qra*(Nra/Vra)*(1-pra))... %inhale
    -(1-u)*(Qra*(Nalv/Valv)*(1-palv)-Qra*(Nra/Vra)*pra-Qra*(Nra/Vra)*(1-pra)); % exhale % NOTE: Am I accounting for what deposited during inhale?

sys(13) = u*Qra*(Nra/Vra)*pra-(1-u)*Qra*(Nra/Vra)*pra;

% lung generation 25 (alveoli)
sys(14) = u*(Qra*(Nra/Vra)*(1-pra)-Qra*(Nalv/Valv)*palv)... %inhale
    -(1-u)*(-Qra*(Nalv/Valv)*palv-Qra*(Nalv/Valv)*(1-palv)); % exhale % NOTE: Am I accounting for what deposited during inhale?

sys(15) = u*Qra*(Nalv/Valv)*palv-(1-u)*Qra*(Nalv/Valv)*palv;



% sys(6) = u*(Qt*Cin-Qt*(Nt/Vt)*pt-Qb*(Nt/Vt)*(1-pt))... %inhale
%     -(1-u)*(Qb*(Nb/Vb)*(1-pb)-Qb*(Nt/Vt)*pt-Qt*(Nt/Vt)*(1-pt)); %exhale
% 
% sys(7) = u*Qt*(Nt/Vt)*pt-(1-u)*Qb*(Nt/Vt)*pt;
% 
% sys(8) = u*(Qb*(Nt/Vt)*(1-pt)-Qb*(Nb/Vb)*pb-Qc*(Nb/Vb)*(1-pb))... %inhale
%     -(1-u)*(Qc*(Nc/Vc)*(1-pc)-Qc*(Nb/Vb)*pb-Qb*(Nb/Vb)*(1-pb)); %exhale
% 
% sys(9) = u*Qb*(Nb/Vb)*pb-(1-u)*Qc*(Nb/Vb)*pb;
% 
% sys(10) = u*(Qc*(Nb/Vb)*(1-pb)-Qc*(Nc/Vc)*pc-Qra*(Nc/Vc)*(1-pc))... %inhale
%     -(1-u)*(Qra*(Nra/Vra)*(1-pra)-Qra*(Nc/Vc)*pc-Qc*(Nc/Vc)*(1-pc)); %exhale
% 
% sys(11) = u*Qc*(Nc/Vc)*pc-(1-u)*Qra*(Nc/Vc)*pc;
% 
% sys(12) = u*(Qra*(Nc/Vc)*(1-pc)-Qra*(Nra/Vra)*pra)... %inhale
%     -(1-u)*(-Qra*(Nra/Vra)*pra-Qra*(Nra/Vra)*(1-pra)); % exhale % NOTE: Am I accounting for what deposited during inhale?
% 
% sys(13) = u*Qra*(Nra/Vra)*pra-(1-u)*Qra*(Nra/Vra)*pra;

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