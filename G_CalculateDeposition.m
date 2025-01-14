function [f,g] = G_CalculateDeposition(t,s,dp,par,res,Y,vol,mech,geom_mode,breaths,output,ratio,Qin,Qex,Vin,Vex)

tt = mod(t,res.T);
for i = 1:length(tt)
    if tt(i)>=0 && tt(i)<=res.TI
        u(i) = 1;
    elseif tt(i)>res.TI && tt(i)<=res.T
        u(i) = 0;
    end
end


convf = 100;

% save each column as array
gnum = Y{:,1}; % generation number
n = Y{:,2}; % number of airways in each generation
Ldata = Y{:,3}/convf; % length of airways in each generation converted FROM centimeters TO meters
rdata = Y{:,4}/2/convf; % radius of airways in each generation (diameter cut in half and converted FROM centimeters TO meters)
V = Y{:,8}./10^6; % volume

% Indices for each compartment
c1_ind = 1; 
c2_ind = 2:5; 
c3_ind = 6:16; 
c4_ind = 17:length(gnum)-1; 
c5_ind = length(gnum); 

% SCALING

% Scale to our desired lung volume to initialize
sc = (vol.TLC/sum(V))^(1/3); % scaled to TLC
r0 = rdata*sc;
L0 = Ldata*sc;

if ratio==1
    VDratio = sum(V([c1_ind c2_ind c3_ind]))/sum(V);
elseif ratio==2
    VDratio = 0.12;
end

VDTLC = VDratio*vol.TLC;
VRTLC = (1-VDratio)*vol.TLC;

Ppl = s(:,1); % intrapleural pressure
Pt = s(:,2); % pressure in larynx
Pb = s(:,3); % pressure in trachea
Pc = s(:,4); % pressure in bronchai
Pra = s(:,5); % pressure in alveoli (lung pressure)

% LUNG VOLUMES

Vt = par.Ct*Pt+par.Vut; % volume of air in trachea
Vb = par.Cb*(Pb-Ppl)+par.Vub; % volume of air in bronchea
Vc = par.Cc*(Pc-Ppl)+par.Vuc; % volume of air in conducting airways
% Vra = par.Cra*(Pra-Ppl)+par.Vura; % volume of air in respiratory airways
% Vtotal = Vt+Vb+Vc+Vra;
VA = par.Cra*(Pra-Ppl)+par.Vura; % volume of air in respiratory airways
Vtotal = Vt+Vb+Vc+VA;

% SCALE GEOMETRY
if geom_mode==1
    cscale = ones(1,length(t));
    rscale = ones(1,length(t));
elseif geom_mode==2

    % scaled to TLC
    
    cscale = ((Vt+Vb+Vc)./VDTLC).^(1/2); % updated
%     rscale = (Vra./VRTLC).^(1/3); % updated
    rscale = (VA./VRTLC).^(1/3); % updated
    
elseif geom_mode==3
    cscale = ones(1,length(t));
%     rscale = (Vra./VRFRC).^(1/3);
    rscale = (VA./VRFRC).^(1/3);
end

r = ones(length(r0),1);
L = ones(length(L0),1);

for i = 1:length(t)
%     convf = 100;
    
    % group generations into compartments
    % compartment 1 consists of the trachea (generation 1 in Yeh)
    r(c1_ind) = r0(c1_ind)*cscale(i); 
    L(c1_ind) = L0(c1_ind);

%     r = r(:);
%     L = L(:);

    % calculate parallel resistance based on Pouseille flow - this is the 
    % resistance in each generation
    mu = par.mu/98.0665; % convert viscosity FROM pascal seconds TO cmh2o seconds (1 cmh2o = 98.0665 Pa)

    R = (8/pi)*mu*(L./r.^4)./n; % cmh2o s m^-3 - resistance in identical parallel tubes
    
    Rt(i) = sum(R(c1_ind));
    

%     Rt = 1.021*1e3;
% Rb = 0.3369*1e3;
% Rc = 0.3063*1e3;
% Rra = 0.0817*1e3;

    % VOLUMETRIC FLOWS
    Qt(i) = (par.Pao-Pt(i))/Rt(i);
end

Qt = Qt(:);
% plot(t,Qt)
% stop


% Particles
Nt = s(:,6);
Dt = s(:,7);
Nb = s(:,8);
Db = s(:,9);
Nc = s(:,10);
Dc = s(:,11);
Nra = s(:,12);
Dra = s(:,13);
Nalv = s(:,14);
Dalv = s(:,15);
Ntotal = Nt+Nb+Nc+Nra+Nalv;
Dtotal = Dt+Db+Dc+Dra+Dalv;

Tind = find(abs(res.T-t)==min(abs(res.T-t))); % index at end of first breath
TIind = find(abs(res.TI-t)==min(abs(res.TI-t))); % index at end of first inhalation
In1 = trapz(t(1:TIind),u(1:TIind)'.*Qt(1:TIind).*par.Cin); % total amount breathed in

% f(1) = Dt(Tind)/In1;
% f(2) = Db(Tind)/In1;
% f(3) = Dc(Tind)/In1;
% f(4) = Dra(Tind)/In1;
% f(5) = Dtotal(Tind)/In1;
% stop
Tind2 = find(abs((breaths-1)*res.T-t)==min(abs((breaths-1)*res.T-t))); % index at the beginning of the last breathing cycle
In2 = In1+Ntotal(Tind2); % The amount available to deposit during last breathing cycle (particle in + particles suspended at beginning of cycle)

%     Out = abs(trapz(t(Tind2:end),(1-u(Tind2:end))'.*Qt(Tind2:end).*(Nt(Tind2:end)./Vt(Tind2:end)).*(1-pt(Tind2:end))));
%     figure(70)
%     hold on
%     plot(t,Dtotal)
%     plot(t(end),Dtotal(end),'ko')
%     plot(t(Tind2),Dtotal(Tind2),'ro')
%     stop
% Fraction of inhaled particles that are deposited after first
% breathing cycle
f(1) = Dt(Tind)/In1;
f(2) = Db(Tind)/In1;
f(3) = Dc(Tind)/In1;
f(4) = Dra(Tind)/In1;
f(5) = Dalv(Tind)/In1;
f(6) = Dtotal(Tind)/In1;

% Deposited at end of last breath
f(7) = (Dt(end)-Dt(Tind2))/In2;%max(Dt)/In;
f(8) = (Db(end)-Db(Tind2))/In2;%max(Db)/In;
f(9) = (Dc(end)-Dc(Tind2))/In2;%max(Dc)/In;
f(10) = (Dra(end)-Dra(Tind2))/In2;%max(Dra)/In;
f(11) = (Dalv(end)-Dalv(Tind2))/In2;
f(12) = (Dtotal(end)-Dtotal(Tind2))/In2;

% Suspended at end of first breath
g(1) = Nt(Tind)/In1;
g(2) = Nb(Tind)/In1;
g(3) = Nc(Tind)/In1;
g(4) = Nra(Tind)/In1;
g(5) = Nalv(Tind)/In1;
g(6) = Ntotal(Tind)/In1;

% Suspended at end of last breath
g(7) = Nt(end)/In2;
g(8) = Nb(end)/In2;
g(9) = Nc(end)/In2;
g(10) = Nra(end)/In2;
g(11) = Nalv(end)/In2;
g(12) = Ntotal(end)/In2;

end