% Function for Janse-Rit model 
function xdot = full_JP(~, x, p)
%Parameters 
A = p.A; %  Maximum size of excitatory postsynaptic potential 
B = p.B; % Maximum size of inhibitory postsynaptic potential 
a = p.a; % Time constant of excitatory postsynaptic potential 
b = p.b; %Time constant of inhibitory postsynaptic potential

c1 = p.c1; % Connection strength parameter from pyramidal to excitatory population
c2 = p.c2; % Connection strength parameter from excitatory to pyramidal population
c3 = p.c3; % Connection strength parameter from pyramidal to inhibitory population 
c4 = p.c4; % Connection strength parameter from inhibitory to pyramidal population

vmax1 = p.vmax1; % Maximum threshold (pyramidal)
vmax2 = p.vmax2; %  Maximum threshold (excitatory)
vmax3 = p.vmax3; %   Maximum threshold (inhibitory)

v01 = p.v01; % Half-maximum activation (pyramidal)
v02 = p.v02; %  Half-maximum activation (excitatory)
v03 = p.v03; % Half-maximum activation (inhibitory)

r1 = p.r1; % Reciprocal of the activation slope (pyramidal)
r2 = p.r2; % Reciprocal of the activation slope (excitatory)
r3 = p.r3; % Reciprocal of the activation slope (inhibitory)

I = p.I; % External input

% Initialize output
xdot = zeros(6,1);


% sigmoid function inline
sgm = @(vmax,v0, scale, v) vmax ./ (1 + exp((v0 - v) .* scale));


xdot(1) = x(4);
xdot(2) = x(5);
xdot(3) = x(6);

xdot(4) = A.*a.*(I+sgm(vmax1,v01,r1,(c2.*x(2)-c4.*x(3))))-2.*a.*x(4)-a^2.*x(1) ;
xdot(5) = A.*a.*(sgm(vmax2,v02,r2,c1.*x(1)))-2.*a.*x(5)-a^2.*x(2);
xdot(6) = B.*b.*(sgm(vmax3,v03,r3,c3.*x(1)))-2.*b.*x(6)-b^2.*x(3);
end

