Rissler = readtable('data.xlsx','Sheet','Rissler','VariableNamingRule','preserve');
Rissler.PD = Rissler.PD/1000; % convert from nanometers to microns

grayscale = [0.5 0.5 0.5];

figure(1)
riss = plot(Rissler.PD,Rissler.DF,'--k>','MarkerSize',6,'Color',grayscale,'MarkerFaceColor',grayscale,'DisplayName','Rissler');
hold on
% riss.MarkerFaceColor='k';

Heyder1986 = readtable('data.xlsx','Sheet','Heyder1986','VariableNamingRule','preserve');
x = Heyder1986.PD;
% y = Heyder1986.RR7_5TV1
% y = Heyder1986.RR15TV0_5;
% y = Heyder1986.RR15TV1_5;
% y = Heyder1986.RR3_75TV1;
% y = Heyder1986.RR3_75TV2;
% y = Heyder1986.RR12TV1_5;
y = Heyder1986.RR15TV1;

hey = plot(x,y,'-.k>','MarkerSize',6,'Color',grayscale,'MarkerFaceColor',grayscale,'DisplayName','Heyder');
% hey.MarkerFaceColor='k';

Dar = readtable('data.xlsx','Sheet','Darquenne2004','VariableNamingRule','preserve');
x = Dar.DP;
y = Dar.DF;
dar = plot(x,y,':k>','MarkerSize',6,'Color',grayscale,'MarkerFaceColor',grayscale,'DisplayName','Darquenne');
% dar.MarkerFaceColor='k';
% set(gca,'XScale','log')
% xlim([0.5 10])

