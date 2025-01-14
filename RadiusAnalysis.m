load('Deposition_results')

f1 = f;
g1 = g;

dp_list1 = linspace(.01,.1,10);
dp_list2 = linspace(.1,1,10);
dp_list3 = linspace(1,10,10);
dp_list = [dp_list1 dp_list2 dp_list3]; % particle sizes in meters

sub = ["t","b","c","ra","total"];
colorvec = ["-bo","-rd","-k^","-m|","--ks"];
markervec = ["o","x","^","|","s"];

% model predicted compartmental deposition
h3 = figure;
hold on
plot(dp_list*10^6,100*f1(1,:),'-k','Marker',markervec(2)','DisplayName',sub(1));
plot(dp_list*10^6,100*f1(2,:),'-k','Marker',markervec(4),'DisplayName',sub(2));
plot(dp_list*10^6,100*f1(3,:),'-k','Marker',markervec(3),'DisplayName',sub(3));
plot(dp_list*10^6,100*f1(4,:),'-k','Marker',markervec(1),'DisplayName',sub(4));
plot(dp_list*10^6,100*f1(5,:),'-k','Marker',markervec(5),'DisplayName',sub(5));
set(gca,'XScale','log')
xlabel('Particle Size ($\mu$m)')
ylabel('Deposition Efficiency, $\eta$ (\%)')
legend location northwest box off
set(gcf,'Color',[1 1 1])
grid on
ylim([0 100])

% Radius and length scaled down 
load('Radius_scaled')

f2 = f;
g2 = g;

% model predicted compartmental deposition
h4 = figure;
hold on
plot(dp_list*10^6,100*f2(1,:),'-k','Marker',markervec(2)','DisplayName',sub(1));
plot(dp_list*10^6,100*f2(2,:),'-k','Marker',markervec(4),'DisplayName',sub(2));
plot(dp_list*10^6,100*f2(3,:),'-k','Marker',markervec(3),'DisplayName',sub(3));
plot(dp_list*10^6,100*f2(4,:),'-k','Marker',markervec(1),'DisplayName',sub(4));
plot(dp_list*10^6,100*f2(5,:),'-k','Marker',markervec(5),'DisplayName',sub(5));
set(gca,'XScale','log')
xlabel('Particle Size ($\mu$m)')
ylabel('Deposition Efficiency, $\eta$ (\%)')
legend location northwest box off
set(gcf,'Color',[1 1 1])
grid on
ylim([0 100])

h5 = figure;
hold on
plot(dp_list*10^6,abs(f2(1,:)-f1(1,:))./f2(1,:),'-k','Marker',markervec(2)','DisplayName',sub(1));
plot(dp_list*10^6,abs(f2(2,:)-f1(2,:))./f2(2,:),'-k','Marker',markervec(4),'DisplayName',sub(2));
plot(dp_list*10^6,abs(f2(3,:)-f1(3,:))./f2(3,:),'-k','Marker',markervec(3),'DisplayName',sub(3));
plot(dp_list*10^6,abs(f2(4,:)-f1(4,:))./f2(4,:),'-k','Marker',markervec(1),'DisplayName',sub(4));
plot(dp_list*10^6,abs(f2(5,:)-f1(5,:))./f2(5,:),'-k','Marker',markervec(5),'DisplayName',sub(5));
set(gca,'XScale','log')
xlabel('Particle Size ($\mu$m)')
ylabel('residual')
legend location northwest box off
set(gcf,'Color',[1 1 1])
grid on

h6 = figure;
hold on
plot(sum((f2(1,:)-f1(1,:)).^2),'-k','Marker',markervec(2),'DisplayName',sub(1));
plot(sum((f2(2,:)-f1(2,:)).^2),'-k','Marker',markervec(4),'DisplayName',sub(2));
plot(sum((f2(3,:)-f1(3,:)).^2),'-k','Marker',markervec(3),'DisplayName',sub(3));
plot(sum((f2(4,:)-f1(4,:)).^2),'-k','Marker',markervec(1),'DisplayName',sub(4));
plot(sum((f2(5,:)-f1(5,:)).^2),'-k','Marker',markervec(5),'DisplayName',sub(5));
% set(gca,'XScale','log')
% xlabel('Particle Size ($\mu$m)')
ylabel('sum of residuals')
legend location northwest box off
set(gcf,'Color',[1 1 1])
grid on