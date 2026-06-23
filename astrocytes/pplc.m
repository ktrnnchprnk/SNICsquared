clc; clearvars
addpath('helpers/')
load('colour.mat')
load("astro.mat")

% Set coupling strengths
p.cae=16.8;
p.cai=9;

% Generate initial guesses for equilibrium search
j.X1=-0.1:0.03:0.45;
j.X2=-0.1:0.03:0.45;
j.X3=-0.1:0.03:0.45;

[j.y1, j.y2, j.y3] = ndgrid(j.X1, j.X2, j.X3);
j.gridpts = [j.y1(:), j.y2(:), j.y3(:)];

% Preallocate arrays for fixed points and solver status
FixedPTS=NaN(length(j.gridpts),3);
FixedPTS2=NaN(length(j.gridpts),3);
Exit=NaN(length(j.gridpts),1);

options = optimset('Display', 'off', 'TolFun', 1e-13, 'TolX', 1e-11);

% Search for equilibria from all grid points
for i = 1:length(j.gridpts)

    [FixedPTS(i,:), ~, Exit(i)] = ...
        fsolve(@(x) astro([], x,p), j.gridpts(i,:), options);

    if Exit(i)==1
        FixedPTS2(i,:)=FixedPTS(i,:);
        % disp(FixedPTS2(i,:))
    end

end

out = NaN(1,1);

% Remove unsuccessful solutions
FixedPTS3=FixedPTS2(~any(isnan(FixedPTS2), 2),:);

% Merge duplicate equilibria
UniqueEig = unique(round(FixedPTS3,4), 'rows');

% Compute Jacobian, eigenvalues, and eigenvectors at each equilibrium
for i=1:size(UniqueEig,1)

    U(i,1)= UniqueEig(i,1);
    U(i,2)= UniqueEig(i,2);
    U(i,3)= UniqueEig(i,3);

    A = jacobianest(@(x) astro([], x,p), UniqueEig(i,:));

    [v, evalue] = eig(A);

    eval(i,:) = [evalue(1,1); evalue(2,2); evalue(3,3)];

    v1(i,:) = [v(:,1)];
    v2(i,:) = [v(:,2)];
    v3(i,:) = [v(:,3)];

end

disp(U)
disp(eval)

%%

% Integrate a trajectory to approach the limit cycle
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);

epsi=[0.01,0.01,0.01];
Tspanf=0:0.01:2500;

% Initial condition
X0u2 = [0;0;0];

% X0u2 = [0.0;0.0;0.0];

[~, Y1]=ode45(@astro,Tspanf,X0u2,opts,p);

%%

close all

f2=figure(2);
f2.Units="centimeters";
f2.OuterPosition = [2 30 12.2 12];

set(gca,'FontSize',16,'fontname','times')

hold on
box on
grid off

ylabel('$I$', Interpreter='latex')
xlabel('$E$', Interpreter='latex')
zlabel('$A$', Interpreter='latex')

% Remove transient and plot the limit cycle in the E-I projection
n=500;
plot(Y1(n:end,1),Y1(n:end,2), ...
    'Color',colour.orange,'LineWidth',2,'LineStyle','-')

% Equilibrium point
n=1;
plot(U(n,1),U(n,2),'o', ...
    'MarkerFaceColor','w', ...
    'MarkerEdgeColor','black', ...
    'MarkerSize',9,'LineWidth',1)

xlim([-0.05 0.55])
ylim([-0.05 0.55])

% Export figure
saveas(f2,'ppLC.svg','svg')

% saveas(f2,'down.fig','fig')