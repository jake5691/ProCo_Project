function main

parameter(); % to load p and oc in workspace

tSpan=[0 100]; % no idea what unit
%n0=reshape(repmat([0.2 0.2 0.05 0.05],p.N,1)',p.n*p.N,1);
n0=zeros(p.N*p.n,1);
n0(1:4)=[0.2 0.2 0.05 0.05];
options = odeset();

[t,n] = ode45(@(t,n)ODE_System(t,n,p,oc),tSpan,n0,options);

plot(t,n)

end % function main