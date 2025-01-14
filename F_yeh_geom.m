function geom = yeh_geom(par)

% morhphology study from Yeh (1980)

% number of tubes in generation
geom.N1 = 1;
geom.N2 = 2;
geom.N3 = 4;
geom.N4 = 8;
geom.N5 = 16;
geom.N6 = geom.N5*2;
geom.N7 = geom.N6*2;
geom.N8 = geom.N7*2;
geom.N9 = geom.N8*2;
geom.N10 = geom.N9*2;
geom.N11 = geom.N10*2;
geom.N12 = geom.N11*2;
geom.N13 = geom.N12*2;
geom.N14 = geom.N13*2;
geom.N15 = geom.N14*2;
geom.N16 = geom.N15*2;
geom.N17 = geom.N16*2;
geom.N18 = geom.N17*2;
geom.N19 = geom.N18*2;
geom.N20 = geom.N19*2;
geom.N21 = geom.N20*2;
geom.N22 = geom.N21*2;
geom.N23 = geom.N22*2;

% vessel radius in cm
geom.r1 = 2.01/2;
geom.r2 = 1.56/2;
geom.r3 = 1.13/2;
geom.r4 = 0.827/2;
geom.r5 = 0.651/2;
geom.r6 = 0.574/2;
geom.r7 = 0.435/2;
geom.r8 = 0.373/2;
geom.r9 = 0.322/2;
geom.r10 = 0.257/2;
geom.r11 = 0.198/2;
geom.r12 = 0.156/2;
geom.r13 = 0.118/2;
geom.r14 = 0.092/2;
geom.r15 = 0.073/2;
geom.r16 = 0.06/2;
geom.r17 = 0.054/2;
geom.r18 = 0.05/2;
geom.r19 = 0.047/2;
geom.r20 = 0.045/2;
geom.r21 = 0.044/2;
geom.r22 = 0.044/2;
geom.r23 = 0.043/2;

% vessel length in cm
geom.L1 = 10;
geom.L2 = 4.36;
geom.L3 = 1.78;
geom.L4 = 0.965;
geom.L5 = 0.995;
geom.L6 = 1.01;
geom.L7 = 0.89;
geom.L8 = 0.962;
geom.L9 = 0.867;
geom.L10 = 0.667;
geom.L11 = 0.556;
geom.L12 = 0.446;
geom.L13 = 0.359;
geom.L14 = 0.275;
geom.L15 = 0.212;
geom.L16 = 0.168;
geom.L17 = 0.134;
geom.L18 = 0.12;
geom.L19 = 0.092;
geom.L20 = 0.08;
geom.L21 = 0.07;
geom.L22 = 0.063;
geom.L23 = 0.057;

% parallel resistances in each generation
geom.R1 = (8/pi)*(par.mu/98.0665)*((geom.L1/geom.r1^4)*1000)/geom.N1; % generation 0
geom.R2 = (8/pi)*(par.mu/98.0665)*((geom.L2/geom.r2^4)*1000)/geom.N2; % generation 1
geom.R3 = (8/pi)*(par.mu/98.0665)*((geom.L3/geom.r3^4)*1000)/geom.N3;
geom.R4 = (8/pi)*(par.mu/98.0665)*((geom.L4/geom.r4^4)*1000)/geom.N4;
geom.R5 = (8/pi)*(par.mu/98.0665)*((geom.L5/geom.r5^4)*1000)/geom.N5;
geom.R6 = (8/pi)*(par.mu/98.0665)*((geom.L6/geom.r6^4)*1000)/geom.N6;
geom.R7 = (8/pi)*(par.mu/98.0665)*((geom.L7/geom.r7^4)*1000)/geom.N7;
geom.R8 = (8/pi)*(par.mu/98.0665)*((geom.L8/geom.r8^4)*1000)/geom.N8;
geom.R9 = (8/pi)*(par.mu/98.0665)*((geom.L9/geom.r9^4)*1000)/geom.N9;
geom.R10 = (8/pi)*(par.mu/98.0665)*((geom.L10/geom.r10^4)*1000)/geom.N10;
geom.R11 = (8/pi)*(par.mu/98.0665)*((geom.L11/geom.r11^4)*1000)/geom.N11;
geom.R12 = (8/pi)*(par.mu/98.0665)*((geom.L12/geom.r12^4)*1000)/geom.N12;
geom.R13 = (8/pi)*(par.mu/98.0665)*((geom.L13/geom.r13^4)*1000)/geom.N13;
geom.R14 = (8/pi)*(par.mu/98.0665)*((geom.L14/geom.r14^4)*1000)/geom.N14;
geom.R15 = (8/pi)*(par.mu/98.0665)*((geom.L15/geom.r15^4)*1000)/geom.N15;
geom.R16 = (8/pi)*(par.mu/98.0665)*((geom.L16/geom.r16^4)*1000)/geom.N16;
geom.R17 = (8/pi)*(par.mu/98.0665)*((geom.L17/geom.r17^4)*1000)/geom.N17;
geom.R18 = (8/pi)*(par.mu/98.0665)*((geom.L18/geom.r18^4)*1000)/geom.N18;
geom.R19 = (8/pi)*(par.mu/98.0665)*((geom.L19/geom.r19^4)*1000)/geom.N19;
geom.R20 = (8/pi)*(par.mu/98.0665)*((geom.L20/geom.r20^4)*1000)/geom.N20;
geom.R21 = (8/pi)*(par.mu/98.0665)*((geom.L21/geom.r21^4)*1000)/geom.N21;
geom.R22 = (8/pi)*(par.mu/98.0665)*((geom.L22/geom.r22^4)*1000)/geom.N22;
geom.R23 = (8/pi)*(par.mu/98.0665)*((geom.L23/geom.r23^4)*1000)/geom.N23;

% Lumped resistances (series resistances) - compartments based on Fukui (1972)

% traches
geom.Rt = geom.R1;

% main bronchus and bronchioles (generations 1-4)
geom.Rb = geom.R2+geom.R3+geom.R4+geom.R5;

% conductive and transitory airways (generations 5-15)
geom.Rc = geom.R6+geom.R7+geom.R8+geom.R9+geom.R10+geom.R11+geom.R12+...
    geom.R13+geom.R14+geom.R15+geom.R16;

% respiratory airways (generations 16-23)
geom.Rra = geom.R17+geom.R18+geom.R19+geom.R20+geom.R21+geom.R22+geom.R23;
end

