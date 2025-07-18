function xdot=ML(~,x,p)

% initialise the system

xdot=zeros(2,1);


% parameters
C = p.C; % Capacitance of the cell membrane
gL = p.gL; % leak conductance 
gCa = p.gCa; %calcium conductance 
gK = p.gK; % potassium conductance 
I = p.I; % applied current
VL = p.VL; % equilibium potential of leak channel
VCa = p.VCa; % equilibium potential of calcium channel
VK = p.VK; % equilibium potential of potassium channel
V1 = p.V1; % tuning parameter for steady state and time constant
V2 = p.V2; % tuning parameter for steady state and time constant
V3 = p.V3; % tuning parameter for steady state and time constant
V4 = p.V4; % tuning parameter for steady state and time constant

% auxiliary functions 

% Open state probability functions
function out=Mss(x)
    out = (1+tanh((x-V1)./V2))/2;
end

function out=Nss(x)
    out = (1+tanh((x-V3)./V4))/2;
end
% no name function 
function out=tau(x)
    out =cosh((x-V3)./V4); % shady function
end


xdot(1) = (-gL.*(x(1)-VL)-gCa.*Mss(x(1)).*(x(1)-VCa)-gK.*x(2).*(x(1)-VK)+I)./C;
xdot(2) = tau(x(1)).*(Nss(x(1))-x(2));

end