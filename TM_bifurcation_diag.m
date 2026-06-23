clc; clearvars; close all

%% Load colour map and AUTO-generated bifurcation data
load("functions/colour.mat")      % Custom colour scheme
addpath('auto/TM')                % Path to AUTO output

% Load bifurcation curve data
load('HB1.dat')
load('HB2.dat')
load('SN11.dat')
load('SN12.dat')
load('SN21.dat')
load('SN22.dat')
load('HC11.dat')
load('HC12.dat')

%% Prepare figure
f1 = figure(1);
f1.Units = "centimeters";
f1.OuterPosition = [2 10 12 16];
set(gca, 'FontSize', 16, 'FontName', 'Times');
hold on; box on

%% Plot Hopf bifurcations
nd = 227;
plot(HB1(:,1), HB1(:,5), 'Color', colour.blood, 'LineWidth', 2)              % Supercritical HB
plot(HB2(1:nd,1), HB2(1:nd,5), 'Color', colour.blood, 'LineWidth', 2)        % Supercritical HB
plot(HB2(nd:end,1), HB2(nd:end,5), 'Color', colour.blood, 'LineWidth', 2, 'LineStyle', '--') % Subcritical HB

%% Plot Saddle-node bifurcations
plot(SN11(:,1), SN11(:,5), 'Color', colour.green, 'LineWidth', 2)
plot(SN12(:,1), SN12(:,5), 'Color', colour.green, 'LineWidth', 2)
plot(SN21(:,1), SN21(:,5), 'Color', colour.green, 'LineWidth', 2)
plot(SN22(:,1), SN22(:,5), 'Color', colour.green, 'LineWidth', 2)

%% Plot homoclinic bifurcations and SNICs
ns = 370; n1 = 480;
plot(HC11(1:ns,1), HC11(1:ns,5), 'Color', colour.purple, 'LineWidth', 2)      % Homoclinic
plot(HC11(ns:n1,1), HC11(ns:n1,5), 'Color', colour.grey, 'LineWidth', 2)      % SNIC

nt = 12200;
plot(HC12(1:nt,1), HC12(1:nt,5), 'Color', colour.purple, 'LineWidth', 2)
plot(HC12(nt:end,1), HC12(nt:end,5), 'Color', colour.grey, 'LineWidth', 2)

%% Mark special bifurcation points
plot(SN11(61,1), SN11(61,5), 'd', 'MarkerFaceColor', 'k', 'MarkerSize', 15)   % Cusp
plot(SN21(65,1), SN21(65,5), 'd', 'MarkerFaceColor', 'k', 'MarkerSize', 15)
plot(HB1(end,1), HB1(end,5), 's', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w', 'MarkerSize', 15) % BT
plot(HB2(end,1), HB2(end,5), 's', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w', 'MarkerSize', 15)
plot(HC11(ns,1), HC11(ns,5), 'p', 'MarkerFaceColor', 'k', 'MarkerSize', 15)  % SNSL
plot(HC12(nt,1), HC12(nt,5), 'p', 'MarkerFaceColor', 'k', 'MarkerSize', 15)
plot(HC12(10260,1), HC12(10260,5), 'p', 'MarkerFaceColor', 'w', 'MarkerSize', 15)  % SNIC²
plot(HB2(nd,1), HB2(nd,5), 'h', 'MarkerFaceColor', 'k', 'MarkerSize', 15)     % GH

%% Labels and legend
xlabel('External Input', 'Interpreter', 'latex')
ylabel('Average synaptic strength', 'Interpreter', 'latex')
legend({'super-HB','','sub-HB','','SN','','','','HC','SNIC','','','BT','SNSL','','SNIC$^2$','GH'}, ...
    'Interpreter','latex','Location','north')

xlim([1 3.5])
ylim([14 40])

%% Save figure
saveas(f1, 'figure_components/f3/bif_diagram.svg', 'svg')
