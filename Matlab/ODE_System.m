function dndt=ODE_System(t,nVec,p,oc)

%% values needed for later computations
% making vector to (p.Nxp.n)-Matrix
nMat=reshape(nVec',p.n,p.N)';
% sum over all species in one compartment
nSum=sum(nMat,2)+oc.nAr;
% Temperature
T=(oc.p * p.V * p.epsilon)./( p.N * nSum * p.R); % [12]
% molar fraction x for each compartment
x=nMat./nSum;

%% Computing terms in ODE
% initializing nMatrix and nGenMatrix
nFlowMatrix=NaN(p.N,p.n); % molar stream between compartments
nGenMatrix=NaN(p.N,p.n); % generated amount of species due to reactions

% computing generated amount of species by reaction
[nGenMatrix, r]=nGen(p,oc,T,x);

% computing inlet/outlet stream of each compartment
nFlowMatrix=nFlow(p,oc,nSum,T,x,sum(nGenMatrix,2),r);

% calculating three ODEs (Ar calulated by 'closing condition': sum(x_j)=1)
% for each compartment (x_j: mole fraction of species j)
dndtMat=[oc.nDotFeed*oc.x(1:p.n); nFlowMatrix(1:end-1,:)]-nFlowMatrix+nGenMatrix;


dndt=reshape(dndtMat',p.N*p.n,1); % making matrix to vector

t
end % function ODE_MassBalance

function [nGen, r]=nGen(p,oc,T,x)
%% calculating the produced/reacted amount of species

%% some basics for further calulations
% partial pressure pPart(j) by using molar fraction x for each compartment
pPart=oc.p * x * 9.8692e-6; % CALCULATED IN !!!atm!!!

%% reaction constant backward reaction
% reaction constant for backwards reaction
kMinus=p.k0*exp(-p.E./(p.R*T)); % [5]

%% reaction constant forward reaction
% using Gillespie-Beattie coreelation
alpha= 0.1191849./T + 91.87212./T.^2 + 25122730./T.^4; % [8]
kGBStar= 10.^(-2.69112*log10(T) - 5.51926e-5*T +...
    1.84886e-7*T.^2 + 2001.6./T + 2.6899); % [9]
kGB=kGBStar.*10.^(alpha*p.pStar); % [7]
Kp=kGB.^2;
% finally: reaction constant for forward reaction
kPlus=kMinus.*Kp; % [6]

%% rate of reaction by Temkin-Pyzhev equation
r=p.f/p.rhoC * (kPlus  .* ( pPart(:,2).*pPart(:,1).^1.5 )./pPart(:,3) - ...
    kMinus .*   pPart(:,3)                   ./pPart(:,1).^1.5); % [4]

% rate of genrating moles inside one reactor compartment
nGen=p.eta.*r*p.mC*1/p.N; % [3]

end % function

function nFlow=nFlow(p,oc,nVec,T,x,nGenVec,r)

%% Generating matrix A and vector b
% computing heat capacity of reactor compartment Cp
Cp=nVec*p.cpg + p.Cp; % [22]
%% making vector b
% computing Temperature at inlet of reactor (after heat exchanger heated up gas)
Tinlet=TempInlet(p,oc,T);
% first element: special case because inlet/feed
bOne=Cp(1)*T(1)*nGenVec(1) - p.deltaH*p.mC*nVec(1)*r(1) + oc.nDotFeed *...
    (T(1)*Cp(1) + p.cpg*nVec(1) * (Tinlet - T(1)) ); % [26]
% rest of b
b=T(2:end).*Cp(2:end).*nGenVec(2:end) - p.deltaH*p.mC*nVec(2:end).*r(2:end); % [27]
% 'correcting' b by overwritting falsely computed first element
b=[bOne; b];

%% Making matrix A
% writting entries on main diagonal
A=diag(T.*Cp,0); % [28]
% writting entries on first minor diagonal below
a=-T(2:end).*Cp(2:end) + p.cpg*nVec(2:end) .* (T(1:end-1) - T(2:end)); % [29]
A=A+diag(a,-1);


%% Solving the Linear System of equation: A*nDot=b
%Solving linear system of equation (A*nDot=b)
nDot=linsolve(A,b); % options for linsolve (?maybe upper triangel?)

% molar stream leaving compartment V_i/ comes from V_{i-1}
nFlow=x.*nDot; % [30]
end % function


function Tinlet=TempInlet(p,oc,T)
% Temperatute at inlet of reactor (warmed up in heat exchanger)
Tinlet=(oc.Tin + (p.U*p.A)/(oc.mDot*p.cpg) * T(p.N)) / (1 + (p.U*p.A)/(oc.mDot*p.cpg));
end % function