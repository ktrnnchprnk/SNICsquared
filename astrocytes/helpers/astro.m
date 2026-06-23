function xdot = astro(~,x,p)
% Wilson-Cowan model coupled to an astrocyte population

xdot=zeros(3,1);

% Recurrent and cross-population connection strengths
cee = p.cee;
cii = p.cii;
cei = p.cei;
cie = p.cie;

% Connections from neural populations to astrocytes
cea = p.cea;
cia = p.cia;

% Feedback from astrocytes to neural populations
cae = p.cae;
cai = p.cai;

% Response function gains
ae = p.ae;
ai = p.ai;
aa = p.aa;

% Response function thresholds
thetae = p.thetae;
thetai = p.thetai;
thetaa = p.thetaa;

% External drive and its allocation between populations
Kp = p.Kp;
alpha = p.alpha;

% Total inputs to each population
F_e = cee*x(1) - cie*x(2) + cae*x(3) + Kp.*alpha;
F_i = -cii*x(2) + cei*x(1) + cai*x(3) + Kp.*(1-alpha);
F_a = cea*x(1) - cia*x(2);

% Population response function
function out = phi(a,F,theta)
    out = 1./(1+exp(-a.*(F-theta))) ...
        - 1./(1+exp(a.*theta));
end

% Excitatory population
xdot(1) = -x(1) + (1-x(1)).*phi(ae,F_e,thetae);

% Inhibitory population
xdot(2) = -x(2) + (1-x(2)).*phi(ai,F_i,thetai);

% Astrocyte population
xdot(3) = -x(3) + (1-x(3)).*phi(aa,F_a,thetaa);

end