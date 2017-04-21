%% Basics:
% References:
% main aper: Asanthi Jinasena, Bernt Lie, and Bjørn Glemmestad.
%            Dynamic model of an ammonia synthesis reactor based on open
%            information. University College Southeast Norway, 9th EUROSIM
%            Congress on Modelling and Simulation, September 2016.
% comments:
% -[i] after an equation in matlab file refers to an equation in main paper
function main

parameter(); % to load p and oc in workspaceç

tSpan=[0:0.1:2]; % [h]

% loaad initial values
%load('/Users/PascalBock/Documents/ProCo_Project/Matlab/InitialValuesForThirdRun.mat');

IV=(oc.p*p.V*p.epsilon)/(p.N*p.R*oc.Tin)*oc.x;
%n0=reshape(repmat([[0.214261413392714 0.0769250916777039 0.209163693642374 0.0301320958933694]],p.N,1)',p.n*p.N,1);
%n0=reshape(repmat([0 0 0 0.5503],p.N,1)',p.n*p.N,1);

%n0(1:4)=IV;
n0=reshape(repmat(IV,p.N,1)',p.n*p.N,1);
%n0=zeros(p.N*p.n,1);
%n0(1:4)=[0.3836 0.1321 0.0117 0.0229];
options = odeset();

[t,n] = ode45(@(t,n)ODE_System(t,n,p,oc),tSpan,n0,options);

save('WS.mat','t','n')

figure(1)
subplot(2,3,1)
nt1=reshape(n(1,:),p.n,p.N)';
plot(nt1(:,[1 2]),fliplr([1:150]),'+')
title('t=0h')
yticklabels({'150','100','50','1'})


subplot(2,3,2)
nt2=reshape(n(11,:),p.n,p.N)';
plot(nt2(:,[1 2]),fliplr([1:150]),'+')
title('t=0.5h')
yticklabels({'150','100','50','1'})


subplot(2,3,3)
nt3=reshape(n(21,:),p.n,p.N)';
plot(nt3(:,[1 2]),fliplr([1:150]),'+')
title('t=1h')
yticklabels({'150','100','50','1'})


subplot(2,3,4)
nt1=reshape(n(1,:),p.n,p.N)';
plot(nt1(:,[3 4]),fliplr([1:150]),'+')
title('t=0h')
yticklabels({'150','100','50','1'})


subplot(2,3,5)
nt2=reshape(n(11,:),p.n,p.N)';
plot(nt2(:,[3 4]),fliplr([1:150]),'+')
title('t=0.5h')
yticklabels({'150','100','50','1'})


subplot(2,3,6)
nt3=reshape(n(21,:),p.n,p.N)';
plot(nt3(:,[3 4]),fliplr([1:150]),'+')
title('t=1h')
yticklabels({'150','100','50','1'})



end % function main