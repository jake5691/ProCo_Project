% Matlab Modell
close all; clear; clc;

%% Predifinitions

% 'order' of molecules
% [H2 N2 NH3 Ar]

% naming of variables
% -
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
p.mC    =   125840;     % [kg] : Total mass of catalyst
p.MAr   =   39.95;      % [kg/kmol] : Molar mass of Ar atom
p.MH2   =   2.016;      % [kg/kmol] : Molar mass of H2 molecule
p.MN2   =   28.02;      % [kg/kmol] : Molar mass of N2 molecule
p.MNH3  =   17.034;     % [kg/kmol] : Molar mass of NH3 molecule
p.M     =   [p.MH2 p.MN2 p.MNH3 p.MAr];
p.N     =   150;        % [1] : Number of reactor compartments (Decided after a few trials)
p.n     =   4;          % [1] : Number of  species
p.eta   =   [-3 -1 2 0]; % [1] : Stoichiometric matrix [H2 N2 NH3 Ar]
p.pAtmos=   1.01325e5;  % [Pa] : Atmospheric pressure
p.R     =   8314;       % [J/kmol K] : Universal gas constant
p.rhoC  =   2200;       % [kg/m^3] : Packing density of cat- alyst
p.U     =   1.9296e6;   % [J/(h m^2 K)] Overall heat transfer coefficient
p.V     =   57.2;       % [m^3] : Volume of the reactor
p.DeltaV=   p.V/p.N;    % [m^3] : Volume of compartement

% Operating conditions
oc.mDot =   67.6;       % [kg/s] : Mass flow rate - reactor inlet
oc.p    =   178e5;      % [Pa] : Controlled reactor pres- sure
oc.Tin  =   350;        % [°C] : Feed temperature (heat exchanger inlet)
oc.xH2  =   0.6972;     % [1] : Mole fraction of H2 at reactor inlet
oc.xN2  =   0.24;       % [1] : Mole fraction of N2 at reactor inlet
oc.xNH3 =   0.0212;     % [1] : Mole fraction of NH3 at reactor inlet
oc.xAr  =   1-(oc.xH2 + oc.xN2 + oc.xNH3);
oc.x    =   [oc.xH2 oc.xN2 oc.xNH3 oc.xAr];

% change of unit for Feed Temperature °C->K
oc.Tin=oc.Tin+273.15; % [K] : Feed temperature (heat exchanger inlet)

% avarege molar mass for inlet of reactor
p.Mtilde=0;
for i=1:4
    p.Mtilde=p.Mtilde+sum([oc.x(i)*p.M(i)]);
end % for i
%p.M(5)=p.Mtilde;
clear i;