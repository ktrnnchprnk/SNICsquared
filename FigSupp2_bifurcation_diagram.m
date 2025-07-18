clc; clearvars

load("colour.mat")                  % Custom colour scheme
addpath('auto/ML')                  % Path to AUTO output

% Load AUTO-generated bifurcation data
load('HB1.dat')
load('HB2.dat')
load('SN11.dat')
load('SN12.dat')
load('SN21.dat')
load('SN22.dat')
load('HC11.dat')
load('HC12.dat')
load('GH1.dat')

close all
f1 = figure(1);
f1.Units = "centimeters";
f1.OuterPosition = [2 10 12 16];

set(gca, 'FontSize', 16, 'FontName', 'times');
hold on; box on; grid off

% Plot Hopf bifurcation curves
plot(HB1(:,1), HB1(:,5), 'Color', colour.blood, 'LineWidth', 2)
plot(HB2(:,1), HB2(:,5), 'Color', colour.blood, 'LineWidth', 2)

% Plot Saddle-Node bifurcation curves
plot(SN11(:,1), SN11(:,5), 'Color', colour.green, 'LineWidth', 2)
plot(SN12(:,1), SN12(:,5), 'Color', colour.green, 'LineWidth', 2)
plot(SN21(:,1), SN21(:,5), 'Color', colour.green, 'LineWidth', 2)
plot(SN22(:,1), SN22(:,5), 'Color', colour.green, 'LineWidth', 2)

% Plot Homoclinic (HC) branches
n1 = 1800;
plot(HC11(n1:end,1), HC11(n1:end,5), 'Color', colour.purple, 'LineWidth', 2)
plot(HC11(1:n1,1), HC11(1:n1,5), 'Color', colour.grey, 'LineWidth', 2)

% Special points
nr = 5469;
plot(SN21(nr,1), SN21(nr,5), 'd', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'black', 'MarkerSize', 15, 'LineWidth', 1)

plot(HB1(end,1), HB1(end,5), 's', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 15, 'LineWidth', 1)

n = 1800;
plot(HC11(n,1), HC11(n,5), 'p', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'black', 'MarkerSize', 15, 'LineWidth', 1)

n = 10260;
plot(HC11(n,1), HC11(n,5), 'p', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'black', 'MarkerSize', 15, 'LineWidth', 1)

legend('HB', '', '', 'SN', '', '', ...
    'SNIC', 'HC', 'CP', 'BT', 'SNSL', 'SNIC^2', 'Location', 'northwest')

ylim([0 40])
xlim([-85 60])

ylabel('Maximum potassium conductance')
xlabel('External input')

saveas(f1, 'figure_components/supp_2/bif_diagram.svg', 'svg')
