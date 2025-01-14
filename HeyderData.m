% Yeh

dp_vals = [0.01 0.06 0.2 0.6 1 2 3 4 6 10];

load ResultsFromLit.mat

Yeh750 = results(:,:,1);
Yeh750total = sum(Yeh750);
Yeh1450 = results(:,:,2);
Yeh1450total = sum(Yeh1450);
Yeh2150 = results(:,:,3);
labels = ["(Yeh) NP+TB","(Yeh) P","(Asgharian) TB","(Asgharian) P"...
    "(Yeh) Total","(Asgharian) Total","(Heyder) 3.75 min$^{-1}$",...
    "(Heyder) 30 min$^{-1}$","(Heyder) 15 min$^{-1}$","(Heyder) 7.5 min$^{-1}$"];

% Asgharian

% TB
x1 = [0.014 0.027 0.057 0.107 0.218 0.444 0.888 1.75 3.56 7.011];
y1 = [0.504 0.307 0.185 0.116 0.08 0.062 0.065 0.098 0.203 0.484];

%P

x2 = [0.014 0.027 0.054 0.109 0.22 0.438 0.887 1.794 3.506 6.971];
y2 = [0.27 0.373 0.318 0.209 0.131 0.097 0.123 0.249 0.417 0.292];

% x2 = [0.014 0.028 0.055 0.224 0.438 0.857 1.733 3.506 6.971];
% y2 = [0.271 0.372 0.316 0.13 0.099 0.121 0.247 0.417 0.292];

% relative error
% rel_err = abs(y5-f(5,:))./y5
% figure(5)
% hold on
% plot(x5,rel_err,'-ks','DisplayName','model')
% 
% % relative error between Yeh/Schum and Heyder
% rel_err2 = abs(y8-Yeh750total(3:7))./Yeh750total(3:7)
% plot(x5,rel_err2,'-ko','DisplayName','Yeh/Schum')
% legend
% grayscale = [0.5 0.5 0.5];
% 
% % plot Yeh and Schum and Asgharian
% figure(1)
% hold on
% plot(dp_vals,Yeh750(1,:)+Yeh750(2,:),'LineStyle','-.','Color',grayscale,'Marker','o','MarkerFaceColor',grayscale,'DisplayName',labels(1))
% plot(dp_vals,Yeh750(3,:),'LineStyle','-.','Color',grayscale,'Marker','<','MarkerFaceColor',grayscale,'DisplayName',labels(2))
% plot(x1,y1,'LineStyle','--','Color',grayscale,'Marker','o','MarkerFaceColor',grayscale,'DisplayName',labels(3))
% plot(x2,y2,'LineStyle','--','Color',grayscale,'Marker','<','MarkerFaceColor',grayscale,'DisplayName',labels(4))
% set(gca,'XScale','log')
% legend 'box' 'off' 'location' 'northwest'
% grid on
% 
% % model predicted tracheobronchial and pulmonary (t+b and c+ra)
% h1 = gcf;
% hold on
% % plot(dp_list*10^6,100*(f(1,:)+f(2,:)),'-ko','DisplayName','t+b','MarkerFaceColor','k')
% % plot(dp_list*10^6,100*(f(3,:)+f(4,:)),'-k<','DisplayName','c+ra','MarkerFaceColor','k')
% plot(dp_list*10^6,(f(1,:)+f(2,:)+f(3,:)),'-ko','DisplayName','(model) t+b+c','MarkerFaceColor','k')
% plot(dp_list*10^6,(f(4,:)),'-k<','DisplayName','(model) ra','MarkerFaceColor','k')
% xlabel('Particle Size ($\mu$m)')
% ylabel('Deposition Efficiency, $\eta$')
% set(gca,'XScale','log')
% set(gcf,'Color',[1 1 1])
% saveas(h1,'CompareDepositionTheory.png')

H = struct();

%RR = 15; TV = 1.5 (nasal breathing from Heyder 1986)
nasal = [.91 .89 .86 .72 .61 .45 .36 .25 .14 .14 .25 .45 .85 .96 .98 1 1 1 1 1 1 1];

%RR=30;TV=1
H.y1975.tv1.bpm30.particles = [0.2 0.58 1 2 3];
H.y1975.tv1.bpm30.total = [0.08 0.07 0.08 0.24 0.44];

% RR=15 TV=1
H.y1975.tv1.bpm15.particles = [0.2 0.5 1 2 3];
H.y1975.tv1.bpm15.total = [0.14 0.11 0.15 0.34 0.58];

%RR=7.5 TV=1 - this same scenario is reported in their 1986 paper
% H.y1975.tv1.bpm7.particles = [0.2 0.49 1 2 3];
% H.y1975.tv1.bpm7.total = [0.22 0.17 .24 .52 .68];

% RR =3.75 TV=1
H.y1975.tv1.bpm4.particles = [0.2 0.46 1 2 3];
H.y1975.tv1.bpm4.total = [0.33 0.25 .4 .69 .79];

% particle sizes for 1986 paper are uniform across breathing pattterns
H.y1986.particles = [.005 .007 .01 .02 .03 .05 .07 .1 .2 .4 .7 1 2 3 4 5 6 7 8 9 10 12 15];

% interpolate Yeh results
x = dp_vals(1:end-1);
v = sum(Yeh1450(:,1:end-1));
xq = H.y1986.particles(3:end-6);
vq1 = interp1(x,v,xq);

% figure(3)
% % plot(x,v,'-ko')
% % plot(xq,vq1,'-ko');
% hold on
% particle_subset = H.y1986.particles(9:17);
% depo_subset = nasal(9:17);
% % plot(particle_subset,depo_subset,'rs')
% % plot(H.y1986.particles(1:end-1),nasal,'ko')
% lgd = legend;
% set(lgd,'box','off','location','northwest','FontSize',8)% 'box' 'off' 'location' 'northwest'
% % legend.FontSize = 8;
% set(gca,'XScale','log')
% set(gcf,'Color',[1 1 1])
% 
% rel_err = abs(depo_subset-vq1(7:end))./depo_subset;
% plot(particle_subset,rel_err,'bo')
% sum(rel_err.^2)


% stop
%Q=.25 L/s; RR=7.5 bpm; TV=1 L

H.y1986.pattern1.total = [0.82 0.81 0.8 0.74 0.67 0.52...
    0.43 0.34 0.21 0.18 0.19 0.25 0.53 0.67 0.76 0.81 0.85...
    0.89 0.91 0.92 0.93 0.94 0.95];

H.y1986.pattern1.L = [0 0 0 0 0 0 0 0 0.02 0.08 0.15 0.27...
    0.38 0.49 0.58 0.65 0.74 0.84];

H.y1986.pattern1.B = [0 0 0 0 0 0 0 0 0.03 0.07 0.14 0.18 ...
    0.22 0.23 0.22 0.22 0.18 0.11];

H.y1986.pattern1.A = [0.52 0.43 0.34 0.21 0.18 0.19 0.25 ...
    0.53 0.62 0.61 0.52 0.4 0.29 0.19 0.12 0.06 0.02 0];

%Q=.25 L/s; RR=15 bpm; TV=.5 L
H.y1986.pattern2.total = [.67 .65 .62 .52 .44 .33 .27 ...
    .21 .13 .11 .12 .15 .28 .44 .56 .65 .72 .78 .82 .84 .86 .87 .89];

H.y1986.pattern2.L = [0 0 0 0 0 0 0 .02 .08 .16 .24 ...
    .34 .43 .52 .59 .65 .74 .81];

H.y1986.pattern2.B = [0 0 0 0 0 0 0 .01 .04 .07 .11 ...
    .15 .18 .2 .19 .17 .12 .08];

H.y1986.pattern2.A = [.33 .27 .21 .13 .11 .12 .15 .25 ...
    .32 .33 .3 .23 .17 .1 .06 .04 .01 0];

%Q=.75 L/s; RR=15 bpm; TV=1.5 L
H.y1986.pattern3.total = [.87 .86 .84 .72 .61 .45 .36 ...
    .25 .14 .11 .12 .15 .39 .63 .77 .86 .9 .93 .95 .96 .97 .98];

H.y1986.pattern3.L = [0 0 0 0 0 0 0 .01 .08 .21 .4 .52 .61 .69 .77 .82 .89];

H.y1986.pattern3.B = [0 0 0 0 0 0 0 .01 .05 .11 .17 .2 .22 .21 .17 .14 .09];

H.y1986.pattern3.A = [.45 .36 .25 .14 .11 .12 .15 .37 .5 .45 .29 .18 .1 .05 .02 .01 0];

% Q=.125; RR=3.75; TV=1 - compared for 1975 paper
H.y1986.pattern4.total = [.82 .81 .79 .75 .71 .6 .52 .43 .33 .27 .31 .42 .64 .74 .79 .83 .86 .88 .89];

% Q=.25; RR=3.75; TV=2
H.y1986.pattern5.total = [.92 .91 .9 .85 .8 .7 .62 .51 .34 .29 .34 .46 .73 .83 .88 .91 .93 .94 .95];

% Q=.25; RR=5; TV=1.5
H.y1986.pattern6.total = [.89 .88 .86 .82 .75 .62 .53 .43 .28 .23 .26 .35 .66 .78 .85 .88 .9 .91 .96];

% Q = .5; RR = 15; TV=1
H.y1986.pattern7.total = [.76 .75 .72 .64 .57 .44 .35 .26 .15 .12 .12 .15 .37 .56 .73 .82 .87 .9 .93];

grayscale = [0.5 0.5 0.5];
figure(5)
hold on
% plot(H.y1986.particles,H.y1986.pattern1.total,'k+','Color',grayscale,'DisplayName','Q=.25 L/s; RR=7.5 bpm; TV=1 L')
% plot(H.y1986.particles,H.y1986.pattern2.total,'k^','Color',grayscale,'DisplayName','RR=15 bpm; TV=.5 L (oral breathing)')
% plot(H.y1986.particles(1:end-1),nasal,'ko','Color',grayscale,'DisplayName','RR=15 bpm; TV=1.5 L (nasal breathing)')
plot(H.y1986.particles(1:end-1),H.y1986.pattern3.total,'ks','Color',grayscale,'DisplayName','RR=15 bpm; TV=1.5 L (oral breathing)')
% plot(H.y1975.tv1.bpm30.particles,H.y1975.tv1.bpm30.total,'k*','Color',grayscale,'DisplayName','RR=30;TV=1')
% plot(H.y1975.tv1.bpm15.particles,H.y1975.tv1.bpm15.total,'ko','Color',grayscale,'DisplayName','RR=15;TV=1')
% % plot(H.y1975.tv1.bpm7.particles,H.y1975.tv1.bpm7.total,'k+','Color',grayscale,'DisplayName','RR=7.5;TV=1')
% plot(H.y1975.tv1.bpm4.particles,H.y1975.tv1.bpm4.total,'k>','Color',grayscale,'DisplayName','RR=3.75;TV=1')
% plot(H.y1986.particles(1:end-4),H.y1986.pattern4.total,'ks','Color',grayscale,'DisplayName','Q=.125 L/s; RR=3.75 bpm; TV=1 L')
% plot(H.y1986.particles(1:end-4),H.y1986.pattern5.total,'ks','Color',grayscale,'DisplayName','Q=.25 L/s; RR=3.75 bpm; TV=2 L')
% plot(H.y1986.particles(1:end-4),H.y1986.pattern6.total,'ks','Color',grayscale,'DisplayName','Q=.25 L/s; RR=5 bpm; TV=1.5 L')
% plot(H.y1986.particles(1:end-4),H.y1986.pattern7.total,'ks','Color',grayscale,'DisplayName','RR=15 bpm; TV=1 L')

H.y1986.particles = [.005 .007 .01 .02 .03 .05 .07 .1 .2 .4 .7 1 2 3 4 5 6 7 8 9 10 12 15];
my_result = f(5,:);
% my_rel_err = abs(H.y1986.pattern7.total(12:19)-my_result(12:19))./H.y1986.pattern7.total(12:19);
% plot(H.y1986.particles(12:19),H.y1986.pattern7.total(12:19),'ro')
% figure(3)
% plot(H.y1986.particles(12:19),my_rel_err,'rs')
% sum(my_rel_err.^2)
% % ylim = ([0 1.5]);
% xlabel('Particle Size ($\mu$m)')
% ylabel('Deposition Efficiency, $\eta$')
% % title('RR = 15 bpm')
% lgd = legend;
% set(lgd,'box','off','location','northwest','FontSize',8)% 'box' 'off' 'location' 'northwest'
% % legend.FontSize = 8;
% set(gca,'XScale','log')
% set(gcf,'Color',[1 1 1])
% grid on
% xlim([0 10])
% ylim([0 1])