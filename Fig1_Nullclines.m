clc; clearvars; close all

%% Plot nullclines for SNIC² and bistable regimes

%% SNIC² regime

% Define E and I variable ranges
E1 = linspace(0, 0.5, 10000);
I2 = linspace(0, 0.5, 2000);

% E-nullcline parameters
ae = 4;
thetae = 1.25;
cee = 22;
cie = 17.6;
P = 0;

% Compute I-nullcline from E equation
I1 = (cee * E1 + (1 / ae) .* log((1 - 2 * E1) ./ E1) + P - thetae) ./ cie;

% I-nullcline parameters
ai = 1;
thetai = 5.5;
cei = 21;
cii = 0;

% Compute E-nullcline from I equation
E2 = (cii .* I2 + thetai - (1 / ai) .* log((1 - 2 .* I2) ./ I2)) ./ cei;

% Plot nullclines
f1 = figure(1); f1.Units = "centimeters"; f1.OuterPosition = [5 15 11 10]; 
hold on; box on; set(gca, 'FontSize', 16, 'FontName', 'Times');

plot(E1, I1, 'b', 'LineWidth', 1.5);
plot([E1(2), 0], [I1(2), 0.55], 'b', 'LineWidth', 1.5);
plot([E1(end-1), 0.5], [I1(end-1), -0.05], 'b', 'LineWidth', 1.5);

plot(E2, I2, 'r', 'LineWidth', 1.5);

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
xlim([-0.05, 0.55]);
ylim([-0.05, 0.55]);

%% Bistable regime

% Update parameters
cie = 16;
P = -0.3;

% Recompute nullclines
I1 = (cee * E1 + (1 / ae) .* log((1 - 2 * E1) ./ E1) + P - thetae) ./ cie;
E2 = (cii .* I2 + thetai - (1 / ai) .* log((1 - 2 .* I2) ./ I2)) ./ cei;

% Plot nullclines
f2 = figure(2); f2.Units = "centimeters"; f2.OuterPosition = [5 15 11 10]; 
hold on; box on; set(gca, 'FontSize', 16, 'FontName', 'Times');

plot(E1, I1, 'b', 'LineWidth', 1.5);
plot([E1(2), 0], [I1(2), 0.55], 'b', 'LineWidth', 1.5);
plot([E1(end-1), 0.5], [I1(end-1), -0.05], 'b', 'LineWidth', 1.5);

plot(E2, I2, 'r', 'LineWidth', 1.5);

xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
xlim([-0.05, 0.55]);
ylim([-0.05, 0.55]);

%% Save figures
if ~exist('figure_components/f1', 'dir')
    mkdir('figure_components/f1');
end

saveas(f1, 'figure_components/f1/null_SNIC2.svg', 'svg');
saveas(f2, 'figure_components/f1/null_bistable.svg', 'svg');
