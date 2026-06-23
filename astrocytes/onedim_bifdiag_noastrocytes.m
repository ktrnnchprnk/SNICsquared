clc; clearvars
load("colour.mat")
addpath('auto')

% Load equilibrium continuation data
load('r1.dat')

% Separate stable and unstable equilibrium branches
for i = 1:length(r1)
    Stable(i,1) = r1(i,1);
    Unstable(i,1) = r1(i,1);

    if r1(i,6) == 3
        Stable(i,2) = r1(i,3);
        Stable(i,3) = r1(i,4);
        Unstable(i,2) = NaN;
        Unstable(i,3) = NaN;
    else
        Unstable(i,2) = r1(i,3);
        Unstable(i,3) = r1(i,4);
        Stable(i,2) = NaN;
        Stable(i,3) = NaN;
    end
end

% Saddle-node bifurcation points
LP = r1([142,246],[1,3]);

close all
f2 = figure(2);
f2.Units = "centimeters";
f2.OuterPosition = [60 30 15 12];
set(gca,'FontSize',16,'fontname','times');
hold on;

ylabel('Excitatory population activity', Interpreter='latex')
xlabel('External input to the excitatory population', Interpreter='latex')

% Equilibrium branches
plot(Stable(:,1),Stable(:,2), 'Color','k', 'LineWidth',2, 'LineStyle','-')
plot(Unstable(:,1),Unstable(:,2), 'Color','k', 'LineWidth',2, 'LineStyle',':')

% Saddle-node points
plot(LP(:,1),LP(:,2),'^','MarkerFaceColor', ...
    'black','MarkerEdgeColor','black', 'MarkerSize',8,'LineWidth',1)

legend('stable', 'unstable', 'saddle node', ...
    'Location', 'best')

xlim([0 2])

saveas(f2,'bif1.svg', 'svg')