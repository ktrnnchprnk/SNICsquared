function xdot=TM(~,x,p)

% initialise the system

xdot=zeros(2,1);


% PARAMETERS 
omega = p.omega;
tau=p.tau;
taur=p.taur;
I=p.I;
uu=p.uu;
r01=p.r01;
vmax1=p.vmax1;
v01=p.v01;
r02=p.r02;
vmax2=p.vmax2;
v02=p.v02;

%sigmoid
function out=phi(r,vmax, x,v0)
    out = vmax./(1+exp(r.*(v0-x)));
end

% define the system of ODEs

xdot(1)=(-x(1)+x(2).*uu.*omega.*phi(r01,vmax1,x(1),v01)+sqrt(tau)+I)./tau;
xdot(2)=(1-x(2))./taur-uu.*x(2).*phi(r02,vmax2,x(1),v02);
end