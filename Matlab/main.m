%% Basics:
% References:
% main aper: Asanthi Jinasena, Bernt Lie, and Bjørn Glemmestad.
%            Dynamic model of an ammonia synthesis reactor based on open
%            information. University College Southeast Norway, 9th EUROSIM
%            Congress on Modelling and Simulation, September 2016.
% comments:
% -[i] after an equation in matlab file refers to an equation in main paper
function main

parameter(); % to load p and oc in workspace

tSpan=linspace(0,1,11); % no idea what unit
n0=reshape(repmat([0.3836 0.1321 0.0117 0.0229],p.N,1)',p.n*p.N,1);
%n0=zeros(p.N*p.n,1);
%n0(1:4)=[0.3836 0.1321 0.0117 0.0229];
options = odeset();

[t,n] = ode45(@(t,n)ODE_System(t,n,p,oc),tSpan,n0,options);

plot(t,n)

end % function main