clear;%clc;%close all

%% MANUAL USER INPUTS

% choose respiratory settings
% TLC_liters = 6.7; % average TLC reported by Rissler (2017)
% TLC_liters = 5.19;%6.34;%5.6; % TLC of Yeh and Schum
TLC_liters = 5.6;
% TV_liters = 1; % to compare to Yeh
TV_liters = 0.75; % to compare to Heyder 
RR = 15; % enter in breaths/minute
IE = 1; % dimensionless (shifts the transition between inspiration and expiration in Pmus)
% FRC_liters = 2.5;%3.04;%3.3;
% FRC_liters = 3.04;%3.3;
FRC_liters = 3.3;%3.3;

% choose morphometry data
data_tog = 1; % 1 = Yeh, 2 = Weibel

% set maximum particle size
max_size = 10; % enter in microns

% choose number of breaths to simulate for
breaths = 1;

% choose deposition scenario
mech = 4; % 0 = 'no deposition', 1='impaction', 2='sedimentation', 3 = 'diffusion', 4 = all, 5 = 'no diffusion';

% choose constant geometry or dynamic
geom_mode = 2; % 1 = constant, 2 = dynamic, 3 = constant upper and dynamic respiratory

ratio = 1; % 1 = based on morphometric data, 2 = fit to data

% choose desired output
output = 0; % 0 = compare to data 
            % 1 = single breathing cycle dynamics 
            % 2 = efficiency of inhale versus exhale
            % 3 = steady state
            % 4 = particles inhaled on first breath only
            % 5 = dynamic curves
            % 6 = constant airflows and volumes
%% AUTOMATIC SET UP
% Make necessary conversions
% Calculate dependent variables 
% Store volumes in structure called 'vol' 
% Store respiratory parameters in structure called 'res'
% Call the solver

vol.TV = TV_liters/1e3; % convert to cubic meters
vol.TLC = TLC_liters/1e3; % m^3 - total lung capacity ~ 6 liters in healthy adult male
% vol.IC = .6*vol.TLC; % inspiratory capacity ~ 60% of TLC (Stahl 1967)?
% vol.FRC = .6*vol.TLC; % functional residual capacity ~ 40% of TLC
% vol.FRC = .4*vol.TLC; % functional residual capacity ~ 50% based on Table 2 in Rissler (2017)
vol.FRC = FRC_liters/1e3;

res.RR = RR; % breaths/min
res.IE = IE; % dimensionless (shifts the transition between inspiration and expiration in Pmus)
res.T = 60/res.RR; % respiratory period in seconds
res.TE = res.T/(1+res.IE);
res.TI = res.TE*res.IE;
res.tau = res.TE/5; 
slope = -0.107255226657383/1e3; % the relationship between TV and Pmusmin is linear, which I found by plotting TV for several values of Pmusmin
slope = 1.08*slope;
res.Pmusmin = vol.TV/slope;%-5; % cmH2O %*pcf; (this parameter changes the amplitude of Pmus)
% res.Pmusmin = -5; % cmH2O %*pcf; (this parameter changes the amplitude of Pmus)

% MORPHOMETRY

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
par = F_Parameters();

% MAKE CALL TO THE SOLVER--------------------------------------------------
% [f,g] = C_DepositionSolver(vol,res,data_tog,max_size,mech,breaths);
if output==6
    f = ConstantFlow_Solver(vol,res,Y,mech,breaths,geom_mode,ratio);
else
    [t,s,f,g] = C_DepositionSolver(vol,res,Y,max_size,mech,breaths,geom_mode,output,ratio);
end