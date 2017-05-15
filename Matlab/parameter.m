%close all; clear; clc;

%% Parameter and Operating conditions

% Parameter p
p.A     =   238;        % [m^2] : Heat transfer area
p.cpg   =   35500;      % [J/(kmol K)] : Molar heat capacity of gas mixture
p.cpc   =   1100;       % [J/(kg K)] : Specific heat capacity of catalyst
p.Cp    =   138.4e6;    % [J/K] : Total heat capacity of catalyst
p.deltaH=   -92.4e6;    % [J/kmol] : Enthalpy of the reaction
p.E     =   1.98464e8;  % [J/kmol] : Activation energy of re- verse reaction
p.epsilon=  0.42;       % [1] : Void fraction of catalyst
p.f     =   4.75;       % [1] : Catalyst activity factor
p.k0    =   2.5714e16;  % [(kmol atm^(1/2))/(m^3 h)] Pre?exponential factor of reverse reaction
p.mC    =   125840;    % [kg] : Total mass of catalyst
p.MAr   =   39.95;      % [kg/kmol] : Molar mass of Ar atom
p.MH2   =   2.016;      % [kg/kmol] : Molar mass of H2 molecule
p.MN2   =   28.02;      % [kg/kmol] : Molar mass of N2 molecule
p.MNH3  =   17.034;     % [kg/kmol] : Molar mass of NH3 molecule
p.M     =   [p.MH2 p.MN2 p.MNH3 p.MAr];
p.N     =   150;        % [1] : Number of reactor compartments (Decided after a few trials)
p.n     =   4;          % [1] : Number of  species
p.eta   =   [-3 -1 2 0];  % [1] : Stoichiometric matrix [H2 N2 NH3], changed because inert gas does not take part in reaction
p.pAtmos=   1.01325e5;  % [Pa] : Atmospheric pressure
p.R     =   8314;       % [J/kmol K] : Universal gas constant
p.rhoC  =   2200;       % [kg/m^3] : Packing density of cat- alyst
p.U     =   1.9296e6;   % [J/(h m^2 K)] Overall heat transfer coefficient
p.V     =   57.2;       % [m^3] : Volume of the reactor
p.timeStep= 11;         % [1] : Number of timesteps
p.speciesNames={'H2','N2','NH3','Ar'};

% Operating conditions
oc.mDot =   67.6/3600;    % [kg/h] : Mass flow rate - reactor inlet
oc.p    =   178e5;      % [Pa] : Controlled reactor pressure
oc.Tin  =   350+273.15; % [�C] : Feed temperature (heat exchanger inlet)
oc.xH2  =   0.6972;     % [1] : Mole fraction of H2 at reactor inlet
oc.xN2  =   0.24;       % [1] : Mole fraction of N2 at reactor inlet
oc.xNH3 =   0.0212;     % [1] : Mole fraction of NH3 at reactor inlet
oc.xAr  =   1-(oc.xH2 + oc.xN2 + oc.xNH3); % [1] : Mole fraction of NH3 at reactor inlet
oc.x    =   [oc.xH2 oc.xN2 oc.xNH3 oc.xAr];


% molar inlet stream
oc.nDotFeed=oc.mDot/sum(oc.x.*p.M); % [kmol/s]

% dimensionless pressure
p.pStar=oc.p/p.pAtmos; % [10]

% % Initial condition as feed composition
% oc.IC=(oc.p*p.V*p.epsilon)/(p.N*oc.Tin*p.R).*oc.x;
% oc.nAr=[oc.IC(4)];

% Initial condition as Figure 4
% oc.IC(1:3)=0.48;
% oc.IC(4:7)=0.46;
% oc.IC(8:20)=0.43;
% oc.IC(21:24)=0.425;
% oc.IC(25:150)=0.42;

oc.IC(1:3)=linspace(0.49,0.48,3);
oc.IC(4:7)=linspace(0.48,0.46,4);
oc.IC(8:20)=linspace(0.46,0.43,13);
oc.IC(21:24)=linspace(0.43,0.425,4);
oc.IC(25:150)=linspace(0.425,0.42,126);

oc.IC=oc.x.*oc.IC';
%oc.Tic=[440 475 525 540]+273.15;
% amount of Ar at each compartment
% oc.nAr=(oc.p*p.V*p.epsilon)/(p.N.*oc.Tic*p.R).*oc.x(4);
