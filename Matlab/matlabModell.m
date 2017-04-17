% Matlab Modell
close all; clear; clc;

%% Predifinitions

% 'order' of molecules
% [H2 N2 NH3 Ar]

% naming of variables
% -
%% Parameter and Operating conditions

% Parameter p
p.A=238; % [m^] : Heat transfer area
p.cp=35500; % [J/(kmol K)] : Molar heat capacity of gas mixture
p.cpc=1100; % [J/(kg K)] : Specific heat capacity of catalyst
p.Cp=138.4e6; % [J/K] : Total heat capacity of catalyst
p.deltaH=-92.4e6; % [J/kmol] : Enthalpy of the reaction
p.E=1.98464e8; % [J/kmol] : Activation energy of re- verse reaction
p.epsilon=0.42; % [1] : Void fraction of catalyst
p.f=4.75; % [1] : Catalyst activity factor
p.k_0=2.5714e16; % [(kmol atm^(1/2))/(m^3 h)] Pre?exponential factor of reverse reaction
p.mC=125840; % [kg] : Total mass of catalyst
p.M_Ar=39.95; % [kg/kmol] : Molar mass of Ar atom
p.M_H2=2.016; % [kg/kmol] : Molar mass of H2 molecule
p.M_N2=28.02; % [kg/kmol] : Molar mass of N2 molecule
p.M_NH3=17.034; % [kg/kmol] : Molar mass of NH3 molecule
p.M=[p.M_H2 p.M_N2 p.M_NH3 p.M_Ar];
p.N=150; % [1] : Number of reactor com- partments (Decided af- ter a few trials)
p.eta=[-3 -1 2 0]; % [1] : Stoichiometric matrix [H2 N2 NH3 Ar]
p.p_atmos=1.01325e5; % [Pa] : Atmospheric pressure
p.R=8314; % [J/kmol K] : Universal gas constant
p.rho_cata=2200; % [kg/m^3] : Packing density of cat- alyst
p.U=1.9296e6; % [J/(h m^2 K)] Overall heat transfer coefficient
p.V=57.2; % [m^3] : Volume of the reactor
p.DeltaV=p.V/p.N; % [m^3] : Volume of compartement
% Operating conditions
oc.m_dot=67.6; % [kg/s] : Mass flow rate - reactor inlet
oc.p=178e5; % [Pa] : Controlled reactor pres- sure
oc.T_in=350; % [°C] : Feed temperature (heat exchanger inlet)
oc.x_H2_in=0.6972; % [1] : Mole fraction of H2 at reactor inlet
oc.x_N2_in=0.24; % [1] : Mole fraction of N2 at reactor inlet
oc.x_NH3_in=0.0212; % [1] : Mole fraction of NH3 at reactor inlet
oc.x_Ar_in=1-(oc.x_H2_in + oc.x_N2_in + oc.x_NH3_in);
oc.x=[oc.x_H2_in oc.x_N2_in oc.x_NH3_in oc.x_Ar_in];

% change of unit for Feed Temperature °C->K
oc.T_in=oc.T_in+273.15; % [K] : Feed temperature (heat exchanger inlet)

% avarege molar mass for inlet of reactor
p.Mtilde=0;
for i=1:4
p.Mtilde=p.M_tilde+sum([oc.x(i)*p.M(i)]);
end % for i
p.M(5)=p.M_tilde;
