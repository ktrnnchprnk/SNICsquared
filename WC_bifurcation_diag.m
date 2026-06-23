clc; clearvars

%% Load colour map and bifurcation data
load("colour.mat")              % Custom colour scheme
addpath('auto/WC')              % Path to AUTO bifurcation data

% Load bifurcation data files
load('HB11.dat')   % Hopf bifurcation branch 1
load('HB12.dat')   % Hopf bifurcation branch 2

load('SN11.dat')   % Saddle-node bifurcation branches
load('SN12.dat')
load('SN21.dat')
load('SN22.dat')

load('HC11.dat')   % Homoclinic bifurcation data

%% Plot bifurcation diagram

close all
f1 = figure(1);
f1.Units = "centimeters";
f1.OuterPosition = [2 10 12 16];

set(gca, 'FontSize', 16, 'FontName', 'Times');
hold on; box on; grid off

% Plot saddle-node branches in green
plot(SN11(:,1), SN11(:,5), 'Color', colour.green, 'LineWidth', 2);
plot(SN12(:,1), SN12(:,5), 'Color', colour.green, 'LineWidth', 2);
plot(SN21(:,1), SN21(:,5), 'Color', colour.green, 'LineWidth', 2);
plot(SN22(:,1), SN22(:,5), 'Color', colour.green, 'LineWidth', 2);

% Plot Hopf bifurcations in blood red
plot(HB12(:,1), HB12(:,5), 'Color', colour.blood, 'LineWidth', 2);
plot(HB11(:,1), HB11(:,5), 'Color', colour.blood, 'LineWidth', 2);

% Split homoclinic branch into two color segments
n1 = 68200;
plot(HC11(1:n1,1), HC11(1:n1,5), 'Color', colour.purple, 'LineWidth', 2);
plot(HC11(n1:end,1), HC11(n1:end,5), 'Color', colour.grey, 'LineWidth', 2);

% Highlight special points with markers
plot(HB11(end,1), HB11(end,5), 's', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 15);

plot(SN11(982,1), SN11(982,5), 'd', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 15);
plot(SN11(4701,1), SN11(4701,5), 'd', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 15);

plot(HC11(n1,1), HC11(n1,5), 'p', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 15);
n2 = 45200;
plot(HC11(n2,1), HC11(n2,5), 'p', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 15);

% Axes labels and limits
xlabel('External input')
ylabel('Inhibitory connection strength')
xlim([0.5 1.6])
ylim([3 34])

% Legend entries for plotted branches and markers
legend('SN', '', '', '', 'HB', '', 'HC', 'SNIC', 'BT', '', '', 'SNHL', 'SNIC^2', 'Location', 'northwest')

%% Save figure
saveas(f1, 'figure_components/f2/bif_diagram.svg', 'svg')
