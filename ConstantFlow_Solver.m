function [f,g] = ConstantFlow_Solver(vol,res,Y,mech,breaths,geom_mode,ratio)

% FOR FIGURES
set(groot,'defaultLineLineWidth',1, 'defaulttextInterpreter', 'latex', ...
     'defaultLegendInterpreter', 'latex')
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultAxesFontSize', 18)
set(groot,'defaultFigureColor',[1 1 1])
set(groot,'defaultLegendFontSize',18,'DefaultLegendFontSizeMode','manual')

par = F_Parameters();
Ppl0 = (vol.FRC-(par.Vut+par.Vub+par.Vuc+par.Vura))/(-(par.Cb+par.Cc+par.Cra));

% GET CONSTANT FLOWS AND VOLUMES
IC = [Ppl0 0 0 0 0]; % [Ppl Pt Pb Pc Pra]
dt = 10^-3;
tvals = 0:dt:res.T;
options = odeset('RelTol',10^-10);
[tt,ss] = ode15s(@(tt,ss) D_LungMechanicsODEs(tt,ss,par,Y,res,vol,geom_mode,ratio),tvals,IC,options);
[Qin,Qex,Vin,Vex] = ConstantFlow(tt,ss,Y,vol,par,res,geom_mode,ratio);

% IC = [Ppl0 0 0 0 0 0 0 0 0 0 0 0 0];

% SET TIME STEP------------------------------------------------------------
dt = 10^-3;
% tvals = tinit(end):dt:breaths*res.T;
% tvals = 0:dt:breaths*res.TI;
% tvals = 0:dt:breaths*res.T;

% SOLVE ODE SYSTEM=========================================================
%==========================================================================

dp_list1 = linspace(.01,.1,10);
dp_list2 = linspace(.2,1,9);
dp_list3 = linspace(2,10,9);
dp_list = [dp_list1 dp_list2 dp_list3]*10^-6; % particle sizes in meters

% dp_list = 1*10^-6;

for i = 1:length(dp_list)
    dp = dp_list(i);
    options = odeset('RelTol',10^-10);
    tvals1 = 0:dt:res.T;
%     IC = [Ppl0 0 0 0 0 0 0 0 0 0 0 0 0];
    IC = [0 0 0 0 0 0 0 0 0 0];
    [t1,s1] = ode15s(@(t1,s1) ConstantFlow_DepositionODEs(t1,s1,par,Y,res,vol,dp,mech,...
                geom_mode,ratio,Qin,Qex,Vin,Vex),tvals1,IC,options);

%     figure(1)
%     TIind = find(abs(res.TI-t1)==min(abs(res.TI-t1)));
%     plot(t1(1:TIind),(s1(1:TIind,7)/Vin(4))/par.Cin)
%     hold on
%     plot(t1(TIind:end),(s1(TIind:end,7)/Vex(4))/par.Cin)
%     
%     stop

    breaths = 2;
    tvals2 = t1(end):dt:breaths*res.T;
    IC = s1(end,:);
    diff = 1;
    t = t1;
    s = s1;

    while diff>=.001
        [t2,s2] = ode15s(@(t2,s2) ConstantFlow_DepositionODEs(t2,s2,par,Y,res,vol,dp,mech,...
            geom_mode,ratio,Qin,Qex,Vin,Vex),tvals2,IC,options);
        t = [t; t2(2:end)];
        s = [s; s2(2:end,:)];
        
        Nt = s2(:,1);
        Nb = s2(:,3);
        Nc = s2(:,5);
        Nra = s2(:,7);
        Nalv = s2(:,9);
        Ntotal = Nt+Nb+Nc+Nra+Nalv;

        Tind = abs((breaths-1)*res.T-t)==min(abs((breaths-1)*res.T-t)); % index at end of first breath
%             plot(t,(Ntotal./Vtotal)/par.Cin)
%             hold on
%             plot(t(end),Ntotal(end)/Vtotal(end)/par.Cin,'ko')
%             plot(t(Tind),Ntotal(Tind)/Vtotal(Tind)/par.Cin,'ko')
%             pause
%         diff = abs(Ntotal(end)/Vtotal(end)/par.Cin-Ntotal(Tind)/Vtotal(Tind)/par.Cin);
%             diff=abs((s(end,7)/Vex(4))/par.Cin-(s(Tind,7)/Vex(4))/par.Cin);
          diff=abs((Ntotal(end)/Vex(6))/par.Cin-(Ntotal(Tind)/Vex(6))/par.Cin);
%             figure(1)
%             hold on
%             plot(t,s(:,7))
%             plot(t(end),s(end,7),'ro')
%             plot(t(Tind),s(Tind,7),'bo')
%             pause
        breaths = breaths+1;
        tvals2 = t2(end):dt:breaths*res.T;
        IC = s2(end,:);
    end
    
%     figure(11)
%     plot(t,((s(:,1)+s(:,2))/Vin(1)+(s(:,3)+s(:,4))/Vin(2)+(s(:,5)+s(:,6))/Vin(3))/par.Cin)
%     hold on
% 
%     figure(12)
%     hold on
%     plot(t,((s(:,7)+s(:,8))/Vin(4))/par.Cin)
% %     stop
% 
%     figure(13)
%     hold on
%     plot(t,(s(:,8)/Vin(4))/par.Cin)
%     stop
    saturation(i) = breaths-1;
    [f(:,i),g(:,i)] = ConstantFlow_CalculateDeposition(t,s,dp,par,res,Y,vol,breaths-1,Qin,Qex,Vin,Vex);
end

figure(100)
plot(dp_list(1:2:end)*10^6,saturation(1:2:end),'-kx','MarkerSize',10)
grid on
xlabel('Particle Size ($\mu$m)')
ylabel('Breathing Cycles')
set(gca,'XScale','log')
set(gca,'XTick',[10^-2 10^-1 10^0 10])
set(gcf,'Color',[1 1 1])

% plot total deposition
    figure(1)
    hold on
    % total
%     p1_1=plot(dp_list*10^6,f(5,:),'-k>','MarkerSize',6,'DisplayName','average flow');
    p1_1=plot(dp_list*10^6,f(6,:),'-k>','MarkerSize',6,'DisplayName','constant');
%     p3_2=plot(dp_list(1:2:end)*10^6,f(10,1:2:end),'k>','MarkerSize',6);
    p1_1.MarkerFaceColor='w';
    legend 'box' 'off' 'location' 'northwest'
    grid on
    xlabel('Particle Size ($\mu$m)')
    ylabel('Deposition Fraction, $\eta_{total}$')
    set(gca,'XScale','log')
    set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
    set(gcf,'Color',[1 1 1])
    ylim([0 1])

    % plot tracheobronchial deposition
    figure(2)
    hold on
    % total
    p1_1=plot(dp_list*10^6,f(1,:)+f(2,:)+f(3,:),'-ks','MarkerSize',6,'DisplayName','constant');
%     p3_2=plot(dp_list(1:2:end)*10^6,f(10,1:2:end),'k>','MarkerSize',6);
    p1_1.MarkerFaceColor='k';
    legend 'box' 'off' 'location' 'northwest'
    grid on
    xlabel('Particle Size ($\mu$m)')
    ylabel('$\eta_{tr}+\eta_{b}+\eta_{c}$')
    set(gca,'XScale','log')
    set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
    set(gcf,'Color',[1 1 1])
    ylim([0 1])

    % plot bronchi deposition
    figure(21)
    hold on
    % total
    p1_1=plot(dp_list*10^6,f(1,:)+f(2,:),'-ks','MarkerSize',6,'DisplayName','constant');
%     p3_2=plot(dp_list(1:2:end)*10^6,f(10,1:2:end),'k>','MarkerSize',6);
    p1_1.MarkerFaceColor='k';
    legend 'box' 'off' 'location' 'northwest'
    grid on
    xlabel('Particle Size ($\mu$m)')
    ylabel('$\eta_{tr}+\eta_{b}$')
    set(gca,'XScale','log')
    set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
    set(gcf,'Color',[1 1 1])
    ylim([0 1])

    % plot pulmonary deposition
    figure(3)
    hold on
    % total
    p1_1=plot(dp_list*10^6,f(4,:)+f(5,:),'-ko','MarkerSize',6,'DisplayName','constant');
%     p3_2=plot(dp_list(1:2:end)*10^6,f(10,1:2:end),'k>','MarkerSize',6);
    p1_1.MarkerFaceColor='k';
    legend 'box' 'off' 'location' 'northwest'
    grid on
    xlabel('Particle Size ($\mu$m)')
    ylabel('$\eta_{ra}$')
    set(gca,'XScale','log')
    set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
    set(gcf,'Color',[1 1 1])
    ylim([0 1])

    % plot total suspension
    figure(4)
    hold on
    % total
    p1_1=plot(dp_list*10^6,g(6,:),'-k>','MarkerSize',6,'DisplayName','constant');
%     p3_2=plot(dp_list(1:2:end)*10^6,f(10,1:2:end),'k>','MarkerSize',6);
    p1_1.MarkerFaceColor='k';
    legend 'box' 'off' 'location' 'northwest'
    grid on
    xlabel('Particle Size ($\mu$m)')
    ylabel('$\alpha_{total}$')
    set(gca,'XScale','log')
    set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
    set(gcf,'Color',[1 1 1])
    ylim([0 1])

    % plot tracheobronchial suspension
    figure(5)
    hold on
    % total
    p1_1=plot(dp_list*10^6,g(1,:)+g(2,:)+g(3,:),'-ks','MarkerSize',6,'DisplayName','constant');
%     p3_2=plot(dp_list(1:2:end)*10^6,f(10,1:2:end),'k>','MarkerSize',6);
    p1_1.MarkerFaceColor='k';
    legend 'box' 'off' 'location' 'northwest'
    grid on
    xlabel('Particle Size ($\mu$m)')
    ylabel('$\alpha_{tr}+\alpha_{b}+\alpha_{c}$')
    set(gca,'XScale','log')
    set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
    set(gcf,'Color',[1 1 1])
    ylim([0 1])

    % plot pulmonary suspension
    figure(6)
    hold on
    % total
    p1_1=plot(dp_list*10^6,g(4,:)+g(5,:),'-ko','MarkerSize',6,'DisplayName','constant');
%     p3_2=plot(dp_list(1:2:end)*10^6,f(10,1:2:end),'k>','MarkerSize',6);
    p1_1.MarkerFaceColor='k';
    legend 'box' 'off' 'location' 'northwest'
    grid on
    xlabel('Particle Size ($\mu$m)')
    ylabel('$\alpha_{ra}$')
    set(gca,'XScale','log')
    set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
    set(gcf,'Color',[1 1 1])
    ylim([0 1])

    % plot total retention
    figure(7)
    hold on
    % total
    p1_1=plot(dp_list*10^6,f(4,:)+f(5,:)+g(4,:)+g(5,:),'-ko','MarkerSize',6,'DisplayName','constant');
%     p3_2=plot(dp_list(1:2:end)*10^6,f(10,1:2:end),'k>','MarkerSize',6);
    p1_1.MarkerFaceColor='k';
    legend 'box' 'off' 'location' 'northeast'
    grid on
    xlabel('Particle Size ($\mu$m)')
    ylabel('$\eta_{ra}+\alpha_{ra}$')
    set(gca,'XScale','log')
    set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
    set(gcf,'Color',[1 1 1])
    ylim([0 1])

    % plot tracheobronchial retention
%     figure(8)
%     hold on
%     % total
%     p1_1=plot(dp_list*10^6,f(1,:)+f(2,:)+f(3,:)+g(1,:)+g(2,:)+g(3,:),'-ks','MarkerSize',6,'DisplayName','constant');
% %     p3_2=plot(dp_list(1:2:end)*10^6,f(10,1:2:end),'k>','MarkerSize',6);
%     p1_1.MarkerFaceColor='k';
%     legend 'box' 'off' 'location' 'northwest'
%     grid on
%     xlabel('Particle Size ($\mu$m)')
%     ylabel('$\eta_{tr}+\eta_{b}+\eta_{c}+\alpha_{tr}+\alpha_{b}+\alpha_{c}$')
%     set(gca,'XScale','log')
%     set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
%     set(gcf,'Color',[1 1 1])
%     ylim([0 1])

    % plot pulmonary suspension
%     figure(9)
%     hold on
%     % total
%     p1_1=plot(dp_list*10^6,f(4,:)+g(4,:),'-ko','MarkerSize',6,'DisplayName','constant');
% %     p3_2=plot(dp_list(1:2:end)*10^6,f(10,1:2:end),'k>','MarkerSize',6);
%     p1_1.MarkerFaceColor='k';
%     legend 'box' 'off' 'location' 'northwest'
%     grid on
%     xlabel('Particle Size ($\mu$m)')
%     ylabel('$\eta_{ra}+\alpha_{ra}$')
%     set(gca,'XScale','log')
%     set(gca, 'YTick',0:0.1:1,'XTick',[10^-2 10^-1 10^0 10])
%     set(gcf,'Color',[1 1 1])
%     ylim([0 1])
end
