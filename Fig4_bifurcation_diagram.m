clc; clearvars

% Load colour map and AUTO-generated bifurcation data
load("colour.mat")              % Custom colour scheme
addpath('auto/JR')              % Path to AUTO output for Jansen-Rit model

load('HB1.dat')                 % Hopf bifurcation branch 1
load('HB2.dat')                 % Hopf bifurcation branch 2

load('SN11.dat')                % Saddle-node branch 1
load('SN12.dat')                % Saddle-node branch 2
load('SN21.dat')                % Saddle-node branch 3
load('SN22.dat')                % Saddle-node branch 4

load('HC11.dat')                % Heteroclinic connection branch 1
load('HC12.dat')                % Heteroclinic connection branch 2

close all
f1 = figure(1);
f1.Units = "centimeters";
f1.OuterPosition = [2 10 12 16];

set(gca, 'FontSize', 16, 'FontName', 'Times');
hold on; box on; grid off

% Plot HB2 branch (Hopf)
plot(HB2(:,1), HB2(:,8), 'Color', colour.blood, 'LineWidth', 2)

% Plot SN branches (Saddle-nodes)
plot(SN11(:,1), SN11(:,8), 'Color', colour.green, 'LineWidth', 2)
plot(SN12(:,1), SN12(:,8), 'Color', colour.green, 'LineWidth', 2)
plot(SN21(:,1), SN21(:,8), 'Color', colour.green, 'LineWidth', 2)
plot(SN22(:,1), SN22(:,8), 'Color', colour.green, 'LineWidth', 2)

% Plot HC11 as two segments: oscillatory (purple) and unstable (grey)
nt = 2100;
nt2 = 2390;
plot(HC11(1:nt,1), HC11(1:nt,8), 'Color', colour.purple, 'LineWidth', 2)
plot(HC11(nt:nt2,1), HC11(nt:nt2,8), 'Color', colour.grey, 'LineWidth', 2)

% Plot marker for endpoint of HB2
plot(HB2(end,1), HB2(end,8), 's', 'MarkerFaceColor', ...
    'k', 'MarkerEdgeColor', 'w', 'MarkerSize', 13, 'LineWidth', 1)

% Marker at HC11 transition
plot(HC11(nt,1), HC11(nt,8), 'p', 'MarkerFaceColor', ...
    'k', 'MarkerEdgeColor', 'black', 'MarkerSize', 15, 'LineWidth', 1)

% Marker on SNIC segment
ns = 1200;
plot(HC11(ns,1), HC11(ns,8), 'p', 'MarkerFaceColor', ...
    'w', 'MarkerEdgeColor', 'black', 'MarkerSize', 15, 'LineWidth', 1)

% Marker on SN22 branch
ns = 1620;
plot(SN22(ns,1), SN22(ns,8), 'd', 'MarkerFaceColor', ...
    'k', 'MarkerEdgeColor', 'black', 'MarkerSize', 13, 'LineWidth', 1)

ylabel('Inhibitory connection strength')
xlabel('External input')

legend('HB','','','SN','', ...
    'SNIC','', '','','SNIC^2', '', 'Location', 'northwest')

xlim([-0.05 0.05])
ylim([-0.01 1])

saveas(f1, 'figure_components/f4/bif_diagram.svg', 'svg')
