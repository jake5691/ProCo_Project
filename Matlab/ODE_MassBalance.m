%% ODE-System: Mass balances
% mole balances of reactor
function dndt=ODE_MassBalance(p,oc)

for k=1:N % over all volume compartements
    for l=1:3 % over all components (Ar is inert, so no reaction)
        if k==1
            dndt(l,k)=nInlet(p,oc)+nReac(p,oc,l,T_r);
        else
            dndt(k)=nReac(p,oc,l,T_r);
        end % if
    end % for l
end % for k

end %function ODE_MassBalance

function out=nInlet(p,oc)

out=(oc.x(l)*p.M(l)/p.M(5)*oc.m_dot)/p.M(l);

end %function


function out=nReac(p,oc,l)
%% Calculate Rate Constants kPlus/kMinus
[kPlus,KMinus]=CalcRateConst(p,oc);

%% Reaction Rates r
% calculating reaction rate r [4]
r=p.f/p.rho_cata * (kPlus * (p(1)*p(2)^1.5) / p(3) - kMinus * p(3)/p(2)^1.5);

%% Generation rate n_dot
% calculation generation rate of component l={1,2,3} [3]
out=p.eta(l)*r*p.m_c*1/p.N;

end % function


function [kPlus,KMinus]=CalcRateConst(p,oc)
%% Calculate Temperature
% calculate Temperature in compartements
T_r=CalcTemp(p,oc);

%% Rate Constants kPLus/kMinus
% rate constant backwards [5]
kMinus=p.k_0*exp(-p.E/(p.R*T)); % define T

% calculating rate constant forward reaction by calculaition of:
%       T_reactor->[alpha, K_GBStar, pStar]->K_GB->Kp->KPlus

% calculation alpha [8]
alpha=0.1191849/T_r+91.87212/T_r^2 + 25122730/T_r^4;
% calculating K_GBStar [9]
K_GBstar=10^(-2.69112*log(T_r) - 5.51926e-5*T_r + 1.84886e-7 T_r^2 + 2001.6/T_r + 2.6899);
% Calculating pStar [10]
pStar=oc.p/p.p_atmos;

% calculating K_GB [7]
K_GB=K_GBstar*10^(alpha*pStar);

% Calculating Kp [11]
Kp=K_GB^2;

% rate constant forward [6]
kPlus=kMinus*Kp;
end % function

function T_r=CalcTemp(p,oc)

% Calculating total number of moles in reactor volume
n_r=

% Calculating Temperature in compartement [12, 13]
T_r=(oc.p * p.DeltaV)/(n_r * p.R) * p.epsilon

end % function