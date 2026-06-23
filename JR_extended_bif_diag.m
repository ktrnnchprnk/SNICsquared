clc; clearvars; close all

% Load custom color scheme
load("colour.mat")

% Add path to AUTO-generated bifurcation data
addpath('auto/JR')

% Load Hopf bifurcation branches
load('HB1.dat')   % Hopf branch 1 (not plotted here)
load('HB2.dat')   % Hopf branch 2

% Load saddle-node bifurcation branches
load('SN11.dat')  % Saddle-node 1
load('SN12.dat')  % Saddle-node 2
load('SN21.dat')  % Saddle-node 3
load('SN22.dat')  % Saddle-node 4

% Load heteroclinic connection branches
load('HC11.dat')  % Heteroclinic branch 1 (main)
load('HC12.dat')  % Heteroclinic branch 2 (not used)

% Create figure window
f1 = figure(1);
f1.Units = "centimeters";
f1.OuterPosition = [2 10 20 17];
set(gca, 'FontSize', 16, 'FontName', 'Times');
hold on; box on; grid off;

% Plot Hopf bifurcation branch 2
plot(HB2(:,1), HB2(:,8), ...
    'Color', colour.blood, 'LineWidth', 2)

% Plot all Saddle-Node bifurcation branches
plot(SN11(:,1), SN11(:,8), 'Color', colour.green, 'LineWidth', 2)
plot(SN12(:,1), SN12(:,8), 'Color', colour.green, 'LineWidth', 2)
plot(SN21(:,1), SN21(:,8), 'Color', colour.green, 'LineWidth', 2)
plot(SN22(:,1), SN22(:,8), 'Color', colour.green, 'LineWidth', 2)

% Plot Heteroclinic bifurcation branch in two segments:
% Main segment in purple, followed by transition to grey
nt = 2100;
nt2 = 2390;
plot(HC11(1:nt,1), HC11(1:nt,8), ...
    'Color', colour.purple, 'LineWidth', 2)
plot(HC11(nt:nt2,1), HC11(nt:nt2,8), ...
    'Color', colour.grey, 'LineWidth', 2)

% Mark key bifurcation points
plot(HB2(end,1), HB2(end,8), 's', ...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w', ...
    'MarkerSize', 13, 'LineWidth', 1)       

plot(HC11(nt,1), HC11(nt,8), 'p', ...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', ...
    'MarkerSize', 15, 'LineWidth', 1)         

ns = 1200;
plot(HC11(ns,1), HC11(ns,8), 'p', ...
    'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', ...
    'MarkerSize', 15, 'LineWidth', 1)         
ns = 1620;
plot(SN22(ns,1), SN22(ns,8), 'd', ...
    'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', ...
    'MarkerSize', 13, 'LineWidth', 1)         

% Axis labels
xlabel('External input')
ylabel('Inhibitory connection strength')

% Legend (placeholder labels – update as needed)
legend('HB', '', '', 'SN', '', ...
       'SNIC', 'HC', 'BT', 'SNSL', 'SNIC^2', 'CP', ...
       'Location', 'northwest')

% Axis limits
xlim([-0.05 0.01])
ylim([-0.1 17.5])

% Save figure
saveas(f1, 'figure_components/supp_1/bif_diagram.svg', 'svg')
