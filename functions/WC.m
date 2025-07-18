function xdot=WC(~,x,par)

% initialise the system

xdot=zeros(2,1);


% PARAMETERS 
cee = par.cee; % self-excitation coupling constant 
cie = par.cie; % 
cei = par.cei; % 
ai = par.ai; % 
ae = par.ae; % 
thetai = par.thetai; % 
thetae = par.thetae;  %
Kp = par.Kp;  %  input 
alpha = par.alpha; % non-dimensional ration constant 

% define the inputs into the sigmoid function 

Fe = cee*x(1)-cie*x(2)+alpha*Kp;
Fi = cei*x(1)+(1-alpha)*Kp;
% define the sigmoid function 

function out=phi(a,F,theta)
    out = 1./(1+exp(-a.*(F-theta)))-1./(1+exp(a.*theta));
end

% define the system of ODEs

xdot(1)=-x(1)+(1-x(1)).*phi(ae,Fe,thetae); %e
xdot(2)=-x(2)+(1-x(2)).*phi(ai,Fi,thetai); %i
end
