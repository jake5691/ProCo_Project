function dndt=ODE_System(~,n,p,oc)

% special case: inlet
p.nDotFeed=oc.mDot/sum(oc.x.*p.M) * oc.x;

nMatrix=SolveLSE(p,oc,n);
nGenMatrix=nRG(p,oc,n);
% for both matrices:
% row: 1:N : volume compartment
% column: 1:3 : species



for i=1:p.N
    for j=1:p.n
        if i==1 % special case: inlet feed
            dndtMat(i,j)=p.nDotFeed(j)-nMatrix(i,j)+nGenMatrix(i,j);
        else % i>2
            dndtMat(i,j)=nMatrix(i-1,j)-nMatrix(i,j)+nGenMatrix(i,j);
        end % if
    end % for j
end % for i

dndt=reshape(dndtMat',p.N*p.n,1);

end % function ODE_MassBalance


function nDotrj=SolveLSE(p,oc,nVec)
%% Calculation of Linear System of Equation b=A* \dot{n}^r
b=NaN(p.N,1);
A=zeros(p.N,p.N);

% Calculating amount in volume compartment V_i in a vector, where row i is
% the value for V_i
nGen=sum(nRG(p,oc,nVec),2); % adding up rows

% making b vector
n=nVec(1:p.n);
T=TCalc(p,oc,n);
Tfeed=TfeedCalc(p,oc,nVec);
b(1)=T*p.Cp*nGen(1) - p.deltaH*nGes(n)*r(p,oc,n)*p.mC + sum(p.nDotFeed) *...
    (T*p.Cp + sum(p.nDotFeed)*p.Mtilde+nGes(n)*(Tfeed-T));
for i=2:p.N
    n=nVec(p.n*i-(p.n-1):p.n*i);
    b(i)=TCalc(p,oc,n)*p.Cp*nGen(i) - p.deltaH*nGes(n)*r(p,oc,n)*p.mC;
end % for

for i=1:p.N
    n=nVec(p.n*i-(p.n-1):p.n*i); % amount of species in V_i
    T=TCalc(p,oc,n); % Temperature at V_i
    A(i,i)=T*p.Cp; % main diagonal of matrix
    if i>1
        nMinus=nVec(p.n*(i-1)-(p.n-1):p.n*(i-1)); % amount of species in V_{i-1}
        TMinus=TCalc(p,oc,nMinus); % Temperature at V_{i-1}
        A(i,i-1)=-T*p.Cp - p.cpg*nGes(n) * (TMinus - T);  % second diagonal band of matrix
    end % if
end % for

%% Solving LSE
nDotr=NaN(p.N,1);
nDotr=linsolve(A,b); % solving equation gives you amount of species at each volume
% compartment
% maybe using OPT=UT: Upper triangular makes sense?

% getting value for each species by using mole fraction x
nDotrj=NaN(p.N,p.n);
for i=1:p.N
    n=nVec(p.n*i-(p.n-1):p.n*i);
    nDotrj(i,:)=nDotr(i)*x(n);
end % for i

end




function nRG=nRG(p,oc,nVec)
%% Calculating the generation term n^{r,g}
for i=1:p.N
    for j=1:p.n % for all three species
        n=nVec(p.n*i-(p.n-1):p.n*i);
        nRG(i,j)=p.eta(j)*r(p,oc,n)*p.mC*1/p.N;
    end % for j
end % for i
end % funtion CalcGenMoles

function r=r(p,oc,n)
%% Calculating reaction rate r

% one time calculation of T
T=TCalc(p,oc,n);

% calculation of rate constant backwards reaction k-
kMinus=p.k0* exp(- p.E/(p.R*T));

% Calculating dimensionless pressure
pStar=oc.p/p.pAtmos;

% using Gillespie-Beattie correlation
alpha= 0.1191849/T + 91.87212/T^2 + 25122730/T^4;
KpGBStar= 10^(-2.69112*log10(T) - 5.51926e-5*T + 1.84886e-7*T^2 + 2001.6/T + 2.6899);
KpGB=KpGBStar*10^(alpha*pStar);

% calculation of equilibrium constant
Kp=KpGB^2;

% calculation of rate konstant forward reaction k+
kPlus=kMinus*Kp;

% reaction rate by using Temkin-Pyzhev equation
pPart=pPartCalc(oc,n);
r=p.f/p.rhoC * (kPlus * (pPart(2)*pPart(1)^1.5)/pPart(3) - kMinus * pPart(3)/pPart(1)^1.5 );
end % function r


function T=TCalc(p,oc,n)
%% Calulating Temperature in one volume compartment
T=(oc.p*p.V)/(nGes(n)*p.R*p.N) *p.epsilon; % [K]
end % function T

function Tfeed=TfeedCalc(p,oc,nVec)
%% Calculation the inlet temperature of reactor
% Temperature at the end of reactor
n=nVec(p.n*p.N-(p.n-1):p.n*p.N);
Tout=TCalc(p,oc,n);
% feed temperature after heat exchanger
Tfeed=(oc.Tin + ( (p.U*p.A)/(oc.mDot*p.cpg) * Tout) ) / (1 + (p.U*p.A)/(oc.mDot*p.cpg) );
end % function Tfeed

function pPart=pPartCalc(oc,n)
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




