% Yeh

dp_vals = [0.2 0.6 1 2 3 4 6];

load ResultsFromLit.mat

Yeh750 = results(:,:,1);
Yeh750total = sum(Yeh750);

% only have results up to 6 microns
Yeh1450 = results(:,:,2);
Yeh1450total = sum(Yeh1450(:,3:end-1));
Yeh2150 = results(:,:,3);

% particle sizes for 1986 paper are uniform across breathing pattterns
H_particles = [.2 .4 .7 1 2 3 4 5 6]; % intersecting particle sizes from Heyder

%RR = 15; TV = 1.5 (nasal breathing from Heyder 1986)
nasal = [.14 .14 .25 .45 .85 .96 .98 1 1];

% interpolate Yeh results
x = dp_vals; % exclude 10 microns
v = sum(Yeh1450(:,3:end-1)); % exclude the last column b/c no deposition for 10 microns
xq = H_particles;
vq1 = interp1(x,v,xq);

grayscale = [0.5 0.5 0.5];

figure(3)
hold on
% plot(x,v,'-ko')
plot(xq,vq1,'-.ko','DisplayName','(Yeh) RR=15 bpm; TV=1.45 L (nasal breathing)')
plot(H_particles,nasal,'o','Color',grayscale,'DisplayName','(data) RR=15 bpm; TV=1.5 L (nasal breathing)')
lgd = legend;
set(lgd,'box','off','location','northwest','FontSize',8)% 'box' 'off' 'location' 'northwest'
set(gca,'XScale','log')
set(gcf,'Color',[1 1 1])
xlabel('Particle Size ($\mu$m)')
ylabel('Deposition Efficiency, $\eta$')

rel_err = abs(nasal-vq1)./nasal

figure(4)
hold on
plot(H_particles,rel_err,'-.ko','DisplayName','Yeh')
sum(rel_err.^2)

my_result = f(5,:);

% Q = .5; RR = 15; TV=1
H_pattern7 = [.15 .12 .12 .15 .37 .56 .73 .82 .87];

figure(3)
plot(H_particles,my_result,'-ks','DisplayName','(model) RR=15 bpm; TV=1 L (oral breathing)')
plot(H_particles,H_pattern7,'s','Color',grayscale,'DisplayName','(data) RR=15 bpm; TV=1 L (oral breathing)')

my_rel_err = abs(H_pattern7-my_result)./H_pattern7

figure(4)
plot(H_particles,my_rel_err,'-ks','DisplayName','model')
lgd = legend;
set(lgd,'box','off','location','northeast','FontSize',8)% 'box' 'off' 'location' 'northwest'
set(gca,'XScale','log')
set(gcf,'Color',[1 1 1])
xlabel('Particle Size ($\mu$m)')
ylabel('Relative Error')
sum(my_rel_err.^2)
