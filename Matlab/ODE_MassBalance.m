function ODE_MassBalance






end % function ODE_MassBalance


function CalcLSE

b=zeros(p.N,1);
A=zeros(p.N,p.N);

% making b vector
b(1)=T(p,oc,n)*p.Cp*nGen() - p.deltaH*n*r()*p.mC + p.nFeed *...
    (T(p,oc,n)*p.Cp + p.nDotFeed*p.Mtilde+n*(Tfeed-T(p,oc,n)));
for i=2:p.N
    b(i)=T(p,oc,n)*p.Cp*nGen() - p.deltaH*n*r()*p.mC;
end % for

for i=1:p.N
    A(i,i)=T()*p.Cp;
    if i>1
        A(i,i-1)=-T()*p.Cp - p.cp*nGes(oc,n) * ...
                 (T() - T());
    end % if
end % for
end



function T=T(p,oc,n)
% Calulating Temperature in one volume compartment
T=(oc.p*p.DeltaV)/(nGes(oc,n)*p.R)*p.epsilon;
end % function T

function pPart=pPart(oc,n)
% Calculating partial pressure
pPart=x(oc,n)*oc.p;
end % function pPart

function x=x(oc,n)
% calculating mole fraction
x=nGes(p,oc,n)./n;
end % function x

function nGes=nGes(oc,n)
% calculating overall moles
nGes=sum(n);
end % function nGes




