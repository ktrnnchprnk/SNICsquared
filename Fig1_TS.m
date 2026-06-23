clc; clearvars; close all

% Add function path
addpath(genpath('\home\kn356\Desktop\SNIC2\functions'))

% Load model parameters
load("wc_par_SNIC2_noise.mat")  % Loads struct `p`

% Simulation settings
dt = 0.01;           % Time step
T = 500;            % Total simulation time
X0 = [0.0; 0.0];     % Initial condition
p.noise_std = 0.01;

%% Bistable regime
p.Kp = 1.1;
p.cie =7.8;
[t, X] = simulate_wc_sde(p, X0, dt, T);

if exist('f2','var') && isgraphics(f2)
    close(f2)
end
f2 = figure(2); f2.Units = "centimeters"; f2.OuterPosition = [10 30 16 12];
set(gca, 'FontSize', 16, 'FontName', 'Times'); hold on; box on;
xlabel('Time', 'Interpreter', 'latex')
ylabel('Population activity', 'Interpreter', 'latex')
plot(t, X(1,:), 'k', 'LineWidth', 2)
ylim([-0.1 0.6])
% saveas(f2, 'figure_components/f1/ts_bi.svg', 'svg')

%% Down state
p.Kp = 1.075;
p.cie =7.9;
X0 = [0.1; 0.5];  % Initial condition targeting down state
if exist('f3','var') && isgraphics(f3)
    close(f3)
end

[t, X] = simulate_wc_sde(p, X0, dt, T);

f3 = figure(3); f3.Units = "centimeters"; f3.OuterPosition = [10 30 16 12];
set(gca, 'FontSize', 16, 'FontName', 'Times'); hold on; box on;
xlabel('Time', 'Interpreter', 'latex')
ylabel('Population activity', 'Interpreter', 'latex')
plot(t, X(1,:), 'k', 'LineWidth', 2)
ylim([-0.1 0.6])
% saveas(f3, 'figure_components/f1/ts_down.svg', 'svg')

%% Up state
p.Kp = 1.44;
p.cie =8.35;
p.noise_std = 0.02;
X0 = [0.5; 0.5];
if exist('f4','var') && isgraphics(f4)
    close(f4)
end

[t, X] = simulate_wc_sde(p, X0, dt, T);

f4 = figure(4); f4.Units = "centimeters"; f4.OuterPosition = [10 30 16 12];
set(gca, 'FontSize', 16, 'FontName', 'Times'); hold on; box on;
xlabel('Time', 'Interpreter', 'latex')
ylabel('Population activity', 'Interpreter', 'latex')
plot(t, X(1,:), 'k', 'LineWidth', 2)
ylim([-0.1 0.6])
%
saveas(f4, 'figure_components/f1/ts_up.svg', 'svg')

%% Limit cycle regime
p.Kp = 1.25;
p.cie =9.3;
X0 = [0.0; 0.0];
p.noise_std = 0.01;
T = 100;  % Shorter time for oscillations
if exist('f5','var') && isgraphics(f5)
    close(f5)
end
[t, X] = simulate_wc_sde(p, X0, dt, T);

f5 = figure(5); f5.Units = "centimeters"; f5.OuterPosition = [50 40 16 12];
set(gca, 'FontSize', 16, 'FontName', 'Times'); hold on; box on;
xlabel('Time', 'Interpreter', 'latex')
ylabel('Population activity', 'Interpreter', 'latex')
plot(t, X(1,:), 'k', 'LineWidth', 2)
ylim([-0.1 0.6])
saveas(f5, 'figure_components/f1/ts_LC.svg', 'svg')