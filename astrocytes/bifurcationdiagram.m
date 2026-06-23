clc; clearvars
load("colour.mat")
addpath('auto')

% Load bifurcation curves
load('HB1.dat')
load('HB2.dat')
load('SN11.dat')
load('SN12.dat')
load('SN21.dat')
load('SN22.dat')
load('HC12.dat')

close all
f1 = figure(1);
f1.Units = "centimeters";
f1.OuterPosition = [2 30 12 16];
set(gca, 'FontSize', 16, 'FontName', 'Times');
hold on; box on

% Saddle-node curves
plot(SN11(:,1),SN11(:,6), 'Color',colour.green, 'LineWidth',2, 'LineStyle','-')
plot(SN12(:,1),SN12(:,6), 'Color',colour.green, 'LineWidth',2, 'LineStyle','-')
plot(SN21(:,1),SN21(:,6), 'Color',colour.green, 'LineWidth',2, 'LineStyle','-')
plot(SN22(:,1),SN22(:,6), 'Color',colour.green, 'LineWidth',2, 'LineStyle','-')

% Homoclinic curve and SNIC segment
n = 1600;
plot(HC12(1:n,1),HC12(1:n,6), 'Color',colour.grey, 'LineWidth',2, 'LineStyle','-')

n2 = 37250;
plot(HC12(n:n2,1),HC12(n:n2,6), 'Color',colour.purple, 'LineWidth',2, 'LineStyle','-')

% Hopf curves
plot(HB1(:,1),HB1(:,6), 'Color',colour.red, 'LineWidth',2, 'LineStyle','-')
plot(HB2(:,1),HB2(:,6), 'Color',colour.red, 'LineWidth',2, 'LineStyle','-')

% Codimension-two points
nstar = 8650;
plot(HC12(nstar,1),HC12(nstar,6),'p','MarkerFaceColor', ...
    'w','MarkerEdgeColor','k', 'MarkerSize', 12,'LineWidth',1)

nstar2 = 1600;
plot(HC12(nstar2,1),HC12(nstar2,6),'p','MarkerFaceColor', ...
    'k','MarkerEdgeColor','k', 'MarkerSize', 12,'LineWidth',1)

ylabel('Connection strength from $A$ to $I$', 'Interpreter','latex')
xlabel('Connection strength from $A$ to $E$', 'Interpreter','latex')

legend('SN', '', '', '', ...
    'HC', 'SNIC','', ...
    '','SNIC^2','SNSL', 'Location', 'northwest')

xlim([16.3 17.5]);
ylim([3 20]);

saveas(f1,'bif1.svg', 'svg')
saveas(f1,'bif1.fig', 'fig')