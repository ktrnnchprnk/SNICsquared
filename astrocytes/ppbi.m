clc; clearvars
addpath('helpers/')
load('colour.mat')
load("astro.mat")

% Parameters in the bistable regime
p.cae = 17.2;
p.cai = 20;

% Initial grid for equilibrium search
j.X1 = -0.1:0.03:0.45;
j.X2 = -0.1:0.03:0.45;
j.X3 = -0.1:0.03:0.45;

[j.y1, j.y2, j.y3] = ndgrid(j.X1, j.X2, j.X3);
j.gridpts = [j.y1(:), j.y2(:), j.y3(:)];

FixedPTS = NaN(length(j.gridpts),3);
FixedPTS2 = NaN(length(j.gridpts),3);
Exit = NaN(length(j.gridpts),1);

options = optimset('Display','off','TolFun',1e-13,'TolX',1e-11);

% Locate equilibria from multiple initial conditions
for i = 1:length(j.gridpts)
    [FixedPTS(i,:),~,Exit(i)] = ...
        fsolve(@(x) astro([],x,p), j.gridpts(i,:), options);

    if Exit(i) == 1
        FixedPTS2(i,:) = FixedPTS(i,:);
    end
end

% Remove duplicates
FixedPTS3 = FixedPTS2(~any(isnan(FixedPTS2),2),:);
UniqueEig = unique(round(FixedPTS3,4),'rows');

% Compute eigenvalues and eigenvectors at each equilibrium
for i = 1:size(UniqueEig,1)
    U(i,1) = UniqueEig(i,1);
    U(i,2) = UniqueEig(i,2);
    U(i,3) = UniqueEig(i,3);

    A = jacobianest(@(x) astro([],x,p), UniqueEig(i,:));
    [v,evalue] = eig(A);

    eval(i,:) = [evalue(1,1);evalue(2,2);evalue(3,3)];
    v1(i,:) = [v(:,1)];
    v2(i,:) = [v(:,2)];
    v3(i,:) = [v(:,3)];
end

disp(U)
disp(eval)

%% Integrate trajectories along unstable directions

opts = odeset('RelTol',1e-6,'AbsTol',1e-8);

epsi = [0.01,0.01,0.01];
Tspanf = 0:0.01:2500;
X0u2 = U(2,:) + epsi.*v1(2,:);
[~,Y1] = ode45(@astro,Tspanf,X0u2,opts,p);

epsi = [0.01,0.01,0.01];
X0u2 = U(2,:) - epsi.*v1(2,:);
[~,Y2] = ode45(@astro,Tspanf,X0u2,opts,p);

epsi = [0.01,0.01,0.01];
X0u2 = U(4,:) + epsi.*v1(4,:);
[~,Y3] = ode45(@astro,Tspanf,X0u2,opts,p);

epsi = [0.01,0.01,0.01];
X0u2 = U(4,:) - epsi.*v1(4,:);
[~,Y4] = ode45(@astro,Tspanf,X0u2,opts,p);

%% Plot phase portrait

close all
f2 = figure(2);
f2.Units = "centimeters";
f2.OuterPosition = [60 30 12.2 12];

set(gca,'FontSize',16,'fontname','times');
hold on; box on; grid off

ylabel('$I$', Interpreter='latex')
xlabel('$E$', Interpreter='latex')
zlabel('$A$', Interpreter='latex')

% Trajectories
plot(Y1(:,1),Y1(:,2),'Color',colour.blue,'LineWidth',2)
plot(Y2(:,1),Y2(:,2),'Color',colour.blue,'LineWidth',2)
plot(Y3(:,1),Y3(:,2),'Color',colour.blue,'LineWidth',2)
plot(Y4(:,1),Y4(:,2),'Color',colour.blue,'LineWidth',2)

% Equilibria
n = 1;
plot(U(n,1),U(n,2),'o','MarkerFaceColor', ...
    'k','MarkerEdgeColor','black','MarkerSize',9,'LineWidth',1)

n = 2;
plot(U(n,1),U(n,2),'^','MarkerFaceColor', ...
    'w','MarkerEdgeColor','black','MarkerSize',8,'LineWidth',1)

n = 3;
plot(U(n,1),U(n,2),'o','MarkerFaceColor', ...
    'w','MarkerEdgeColor','black','MarkerSize',9,'LineWidth',1)

n = 5;
plot(U(n,1),U(n,2),'o','MarkerFaceColor', ...
    'k','MarkerEdgeColor','black','MarkerSize',9,'LineWidth',1)

n = 4;
plot(U(n,1),U(n,2),'^','MarkerFaceColor', ...
    'w','MarkerEdgeColor','black','MarkerSize',8,'LineWidth',1)

xlim([-0.05 0.55])
ylim([-0.05 0.55])

saveas(f2,'ppbi.svg','svg')