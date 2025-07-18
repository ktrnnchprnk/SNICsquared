clc; clearvars; close all

%% Add function path
addpath(genpath('/home/kn356/Desktop/SNIC2/functions'))

%% Load model parameters
load("wc_par_SNIC2_noise.mat")  % Loads struct `p`

% Simulation settings
dt = 0.01;           % Time step
T = 1000;            % Total simulation time
X0 = [0.0; 0.0];     % Initial condition

%% SNIC² regime
p.Kp = 1.09208;
p.beta1 = 0.776736;
p.noise_std = 0.001;

[t, X] = simulate_wc_sde(p, X0, dt, T);

f1 = figure(1); f1.Units = "centimeters"; f1.OuterPosition = [2 10 16 12];
set(gca, 'FontSize', 16, 'FontName', 'Times'); hold on; box on;
xlabel('Time', 'Interpreter', 'latex')
ylabel('Population activity', 'Interpreter', 'latex')
plot(t, X(2,:), 'k', 'LineWidth', 2)
ylim([-0.1 0.6])
saveas(f1, 'figure_components/f1/ts_hetero.svg', 'svg')

%% Bistable regime
p.Kp = 1.06;
p.beta1 = 0.78;
p.noise_std = 0.012;

[t, X] = simulate_wc_sde(p, X0, dt, T);

f2 = figure(2); f2.Units = "centimeters"; f2.OuterPosition = [2 10 16 12];
set(gca, 'FontSize', 16, 'FontName', 'Times'); hold on; box on;
xlabel('Time', 'Interpreter', 'latex')
ylabel('Population activity', 'Interpreter', 'latex')
plot(t, X(2,:), 'k', 'LineWidth', 2)
ylim([-0.1 0.6])
saveas(f2, 'figure_components/f1/ts_bi.svg', 'svg')

%% Down state
p.Kp = 1.075;
p.beta1 = 0.65;
p.noise_std = 0.01;
X0 = [0.1; 0.5];  % Initial condition targeting down state

[t, X] = simulate_wc_sde(p, X0, dt, T);

f3 = figure(3); f3.Units = "centimeters"; f3.OuterPosition = [2 10 16 12];
set(gca, 'FontSize', 16, 'FontName', 'Times'); hold on; box on;
xlabel('Time', 'Interpreter', 'latex')
ylabel('Population activity', 'Interpreter', 'latex')
plot(t, X(2,:), 'k', 'LineWidth', 2)
ylim([-0.1 0.6])
saveas(f3, 'figure_components/f1/ts_down.svg', 'svg')

%% Up state
p.Kp = 1.483;
p.beta1 = 0.76354;
p.noise_std = 0.01;

[t, X] = simulate_wc_sde(p, X0, dt, T);

f4 = figure(4); f4.Units = "centimeters"; f4.OuterPosition = [2 10 16 12];
set(gca, 'FontSize', 16, 'FontName', 'Times'); hold on; box on;
xlabel('Time', 'Interpreter', 'latex')
ylabel('Population activity', 'Interpreter', 'latex')
plot(t, X(2,:), 'k', 'LineWidth', 2)
ylim([-0.1 0.6])
saveas(f4, 'figure_components/f1/ts_up.svg', 'svg')

%% Limit cycle regime
p.Kp = 1.3;
p.beta1 = 0.68;
p.noise_std = 0.03;
T = 100;  % Shorter time for oscillations

[t, X] = simulate_wc_sde(p, X0, dt, T);

f5 = figure(5); f5.Units = "centimeters"; f5.OuterPosition = [2 10 16 12];
set(gca, 'FontSize', 16, 'FontName', 'Times'); hold on; box on;
xlabel('Time', 'Interpreter', 'latex')
ylabel('Population activity', 'Interpreter', 'latex')
plot(t, X(2,:), 'k', 'LineWidth', 2)
ylim([-0.1 0.6])
saveas(f5, 'figure_components/f1/ts_LC.svg', 'svg')
