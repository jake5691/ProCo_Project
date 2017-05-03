%% Basics:
% References:
% main aper: Asanthi Jinasena, Bernt Lie, and Bjørn Glemmestad.
%            Dynamic model of an ammonia synthesis reactor based on open
%            information. University College Southeast Norway, 9th EUROSIM
%            Congress on Modelling and Simulation, September 2016.
% comments:
% -[i] after an equation in matlab file refers to an equation in main paper
% -'order' of molecules: [H2 N2 NH3 Ar]
function main

parameter(); % to load p and oc in workspaceç

tSpan=linspace(0,1,p.timeStep); % [h]
n0=reshape(oc.IC',p.n*p.N,1);

options = odeset();

[t,n] = ode23(@(t,n)ODE_System(t,n,p,oc),tSpan,n0,options);

make_Plots(p,oc,t,n);

end % function main