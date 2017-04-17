function dndt=ODE_System(~,n,p,oc)


nMatrix=SolveLSE(p,oc,n);
nGenMatrix=nRG(p,oc,n);
% for both matrices:
% row: 1:N : volume compartment
% column: 1:3 : species

% special case: inlet
p.nDotFeed=oc.mDot/sum(oc.x.*p.M) * oc.x;

for i=1:p.N
    for j=1:3
        if i==1 % special case: inlet feed
            dndt(i,j)=p.nDotFeed(j)-nMatrix(i,j)+nGenMatrix(i,j);
        else % i>2
            dndt(i,j)=nMatrix(i-1,j)-nMatrix(i,j)+nGenMatrix(i,j);
        end % if
    end % for j
end % for i

dndt=reshape(dndt',1,1:p.N);

end % function ODE_MassBalance


function nDotrj=SolveLSE(p,oc,nVec)
%% Calculation of Linear System of Equation b=A* \dot{n}^r
b=NaN(p.N,1);
A=zeros(p.N,p.N);

nGen=sum(nRG(p,oc,nVec));

% making b vector
n=nVec(1:3);
T=TCalc(p,oc,n);
Tfeed=TfeedCalc(p,oc,nVec);
b(1)=T*p.Cp*nGen - p.deltaH*n*r(p,oc,n)*p.mC + p.nFeed *...
    (T*p.Cp + sum(p.nDotFeed)*p.Mtilde+sum(n)*(Tfeed-T));
for i=2:p.N
    n=nVec(3*i-2:3*i);
    b(i)=T(p,oc,n)*p.Cp*nGen - p.deltaH*sum(n)*r(p,oc,n)*p.mC;
end % for

for i=1:p.N
    n=nVec(3*i-2:3*i); % amount of species in V_i
    T=TCalc(p,oc,n); % Temperature at V_i
    A(i,i)=T*p.Cp; % main diagonal of matrix
    if i>1
        nMinus=nVec(3*(i-1)-2:3*(i-1)); % amount of species in V_{i-1}
        TMinus=TCalc(p,oc,nMinus); % Temperature at V_{i-1}
        A(i,i-1)=-T*p.Cp - p.cp*nGes(oc,n) * (TMinus - T);  % second diagonal band of matrix
    end % if
end % for

%% Solving LSE
nDotr=NaN(p.N,1);
nDotr=linsolve(A,b); % solving equation gives you amount of species at each volume
% compartment
% maybe using OPT=UT: Upper triangular makes sense?

% getting value for each species by using mole fraction x
nDotrj=NaN(p.N,3);
for i=1:p.N
    n=nVec(3*i-2:3*i);
    nDotrj(i,:)=nDotr(i)*x(n);
end % for i

end




function nRG=nRG(p,oc,nVec)
%% Calculating the generation term n^{r,g}
for i=1:p.N
    for j=1:3 % for all three species
        n=nVec(3*i-2:3*i);
        nRG(i,j)=p.eta(j)*r(p,oc,n)*p.mC*1/p.N;
    end % for j
end % for i
end % funtion CalcGenMoles

function r=r(p,oc,n)
%% Calculating reaction rate r

% one time calculation of T
T=T(p,oc,n);

% calculation of rate constant backwards reaction k-
kMinus=p.k0* exp(- p.E/(p.R*T);

% Calculating dimensionless pressure
pStar=oc.p/p.pAtmos;

% using Gillespie-Beattie correlation
alpha= 0.1191849/T + 91.87212/T^2 + 25122730/T^4;
KpGBStar= 10^(-2.69112*log(T) - 5.51926e-5*T + 1.84886e-7*T^2 + 2001.6/T + 2.6899);
KpGB=kpGBStar*10^(alpha*pStar);

% calculation of equilibrium constant
Kp=KpGB^2;

% calculation of rate konstant forward reaction k+
kPlus=kMinus*Kp;

% reaction rate by using Temkin-Pyzhev equation
r=p.f/p.rhoC * (kPlus * (pPart(2)*pPart(1)^1.5)/pPart(3) - kMinus * pPart(3)/pPart(1)^1.5 );
end % function r


function T=TCalc(p,oc,n)
%% Calulating Temperature in one volume compartment
T=(oc.p*p.DeltaV)/(nGes(n)*p.R)*p.epsilon;
end % function T

function Tfeed=TfeedCalc(p,oc,nVec)
%% Calculation the inlet temperature of reactor
n=nVec(3*N-2:3*N);
Tout=T(p,oc,n);
Tfeed=(oc.Tin + ( (p.U*p.A)/(oc.mDot*p.cpg) * Tout) ) / (1 + (p.U*p.A)/(oc.mDot*p.cpg) );
end % function Tfeed

function pPart=pPart(oc,n)
%% Calculating partial pressure
pPart=x(n)*oc.p;
end % function pPart

function x=x(n)
%% calculating mole fraction
x=n./nGes(n);
end % function x

function nGes=nGes(n)
%% calculating overall moles
nGes=sum(n);
end % function nGes




