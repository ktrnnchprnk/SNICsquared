clc; clearvars
load("colour.mat")  % Load custom colour scheme
load("functions/parameters/wc_par.mat")  % Load model parameters
addpath("functions")  % Add path to functions

% SNIC^2 phase portrait
fprintf('SNIC^2\n');
p.Kp = 1.09206;
p.cie = (1 - 0.776615) * 35;

[u, e, v1, v2] = compute_fp(@(x) WC([], x, p), [0, 5], [0 5], 0.05, 0.05, 5);
for i = 1:size(u, 1)
    fprintf('Steady State: (%.3f, %.3f), Eigenvalues: (%.3f, %.3f)\n', u(i, 1), u(i, 2), e(i, 1), e(i, 2));
end

% Simulate trajectory from saddle unstable direction
X0u2 = u(1,:) + [0.005, 0.005] .* v1(1,:);
[T, Y1] = ode45(@WC, 0:0.1:600, X0u2, p.opts, p);

close all
f1 = figure(1);
f1.Units = "centimeters";
f1.OuterPosition = [2 10 12.2 12.3];
set(gca, 'FontSize', 16, 'FontName', 'Times'); hold on; box on;
xlabel('Excitatory population activity');
ylabel('Inhibitory population activity');
plot(Y1(:,1), Y1(:,2), 'Color', colour.pink, 'LineWidth', 2);
plot(u([1,3],1), u([1,3],2), '^', 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black', 'MarkerSize', 9);
plot(u(2,1), u(2,2), 'o', 'MarkerFaceColor', 'white', 'MarkerEdgeColor', 'black', 'MarkerSize', 8);
xlim([0 0.5]); ylim([-0.05 0.55]);
saveas(f1, 'figure_components/f2/pp_SNICeroclinic.svg', 'svg');

% Bistability phase portrait
fprintf('Bistability\n');
p.Kp = 1.06;
p.cie = (1 - 0.78) * 35;

[u, e, v1, v2] = compute_fp(@(x) WC([], x, p), [0, 5], [0 5], 0.05, 0.05, 5);
for i = 1:size(u, 1)
    fprintf('Steady State: (%.3f, %.3f), Eigenvalues: (%.3f, %.3f)\n', u(i, 1), u(i, 2), e(i, 1), e(i, 2));
end

% Saddle separatrices (2 saddles at u(2) and u(4))
epsi = [0.005, 0.005];
[T1, Y1] = ode45(@WC, 0:0.01:100, u(2,:) + epsi .* v1(2,:), p.opts, p);
[T2, Y2] = ode45(@WC, 0:0.01:100, u(2,:) - epsi .* v1(2,:), p.opts, p);
[T3, Y3] = ode45(@WC, 10:-0.01:0, u(2,:) - epsi .* v1(2,:), p.opts, p);
[T4, Y4] = ode45(@WC, 10:-0.01:0, u(2,:) + epsi .* v2(2,:), p.opts, p);
[~, Y21] = ode45(@WC, 0:0.01:100, u(4,:) + epsi .* v1(4,:), p.opts, p);
[~, Y22] = ode45(@WC, 0:0.01:100, u(4,:) - epsi .* v1(4,:), p.opts, p);
[~, Y23] = ode45(@WC, 10:-0.01:0, u(4,:) - epsi .* v1(4,:), p.opts, p);
[~, Y24] = ode45(@WC, 10:-0.01:0, u(4,:) - epsi .* v2(4,:), p.opts, p);

f2 = figure(2);
f2.Units = "centimeters";
f2.OuterPosition = [2 10 12 12];
set(gca, 'FontSize', 16, 'FontName', 'Times'); hold on; box on;
xlabel('Excitatory population activity');
ylabel('Inhibitory population activity');
plot(Y1(:,1), Y1(:,2), 'b', 'LineWidth', 2);
plot(Y2(:,1), Y2(:,2), 'b', 'LineWidth', 2);
plot(Y3(:,1), Y3(:,2), 'r', 'LineWidth', 2);
plot(Y4(:,1), Y4(:,2), 'r', 'LineWidth', 2);
plot(Y21(:,1), Y21(:,2), 'b', 'LineWidth', 2);
plot(Y22(:,1), Y22(:,2), 'b', 'LineWidth', 2);
plot(Y23(:,1), Y23(:,2), 'r', 'LineWidth', 2);
plot(Y24(:,1), Y24(:,2), 'r', 'LineWidth', 2);
plot(u([1,5],1), u([1,5],2), 'o', 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black', 'MarkerSize', 9);
plot(u([2,4],1), u([2,4],2), '^', 'MarkerFaceColor', 'white', 'MarkerEdgeColor', 'black', 'MarkerSize', 9);
plot(u(3,1), u(3,2), 'o', 'MarkerFaceColor', 'white', 'MarkerEdgeColor', 'black', 'MarkerSize', 8);
xlim([0 0.5]); ylim([-0.05 0.55]);
saveas(f2, 'figure_components/f2/pp_bistable.svg', 'svg');

% Up state
fprintf('Up State\n');
p.Kp = 1.25;
p.cie = (1 - 0.8) * 35;

[u, e, v1, v2] = compute_fp(@(x) WC([], x, p), [0, 5], [0 5], 0.05, 0.05, 5);
for i = 1:size(u, 1)
    fprintf('Steady State: (%.3f, %.3f), Eigenvalues: (%.3f, %.3f)\n', u(i, 1), u(i, 2), e(i, 1), e(i, 2));
end

[~, Y1] = ode45(@WC, 0:0.01:100, u(2,:) + epsi .* v1(2,:), p.opts, p);
[~, Y2] = ode45(@WC, 0:0.01:100, u(2,:) - epsi .* v1(2,:), p.opts, p);
[~, Y3] = ode45(@WC, 10:-0.01:0, u(2,:) - epsi .* v1(2,:), p.opts, p);
[~, Y4] = ode45(@WC, 10:-0.01:0, u(2,:) - epsi .* v2(2,:), p.opts, p);

f3 = figure(3);
f3.Units = "centimeters";
f3.OuterPosition = [2 10 12 12];
set(gca, 'FontSize', 16, 'FontName', 'Times'); hold on; box on;
xlabel('Excitatory population activity');
ylabel('Inhibitory population activity');
plot(Y1(:,1), Y1(:,2), 'b', 'LineWidth', 2);
plot(Y2(:,1), Y2(:,2), 'b', 'LineWidth', 2);
plot(Y3(:,1), Y3(:,2), 'r', 'LineWidth', 2);
plot(Y4(:,1), Y4(:,2), 'r', 'LineWidth', 2);
plot(u(3,1), u(3,2), 'o', 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black', 'MarkerSize', 9);
plot(u(2,1), u(2,2), '^', 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black', 'MarkerSize', 9);
plot(u(1,1), u(1,2), 'o', 'MarkerFaceColor', 'white', 'MarkerEdgeColor', 'black', 'MarkerSize', 8);
xlim([0 0.5]); ylim([-0.05 0.55]);
saveas(f3, 'figure_components/f2/pp_UP.svg', 'svg');

% Down state
fprintf('Down State\n');
p.Kp = 1;
p.cie = (1 - 0.68) * 35;

[u, e, v1, v2] = compute_fp(@(x) WC([], x, p), [0, 5], [0 5], 0.05, 0.05, 5);
for i = 1:size(u, 1)
    fprintf('Steady State: (%.3f, %.3f), Eigenvalues: (%.3f, %.3f)\n', u(i, 1), u(i, 2), e(i, 1), e(i, 2));
end

[~, Y1] = ode45(@WC, 0:0.01:100, u(2,:) + epsi .* v1(2,:), p.opts, p);
[~, Y2] = ode45(@WC, 0:0.01:100, u(2,:) - epsi .* v1(2,:), p.opts, p);
[~, Y3] = ode45(@WC, 10:-0.01:0, u(2,:) - epsi .* v1(2,:), p.opts, p);
[~, Y4] = ode45(@WC, 10:-0.01:0, u(2,:) + epsi .* v2(2,:), p.opts, p);

f4 = figure(4);
f4.Units = "centimeters";
f4.OuterPosition = [2 10 12.2 12.26];
set(gca, 'FontSize', 16, 'FontName', 'Times'); hold on; box on;
xlabel('Excitatory population activity');
ylabel('Inhibitory population activity');
plot(Y1(:,1), Y1(:,2), 'b', 'LineWidth', 2);
plot(Y2(:,1), Y2(:,2), 'b', 'LineWidth', 2);
plot(Y3(:,1), Y3(:,2), 'r', 'LineWidth', 2);
plot(Y4(:,1), Y4(:,2), 'r', 'LineWidth', 2);
plot(u(1,1), u(1,2), 'o', 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black', 'MarkerSize', 9);
plot(u(2,1), u(2,2), '^', 'MarkerFaceColor', 'black', 'MarkerEdgeColor', 'black', 'MarkerSize', 9);
plot(u(3,1), u(3,2), 'o', 'MarkerFaceColor', 'white', 'MarkerEdgeColor', 'black', 'MarkerSize', 8);
xlim([0 0.5]); ylim([-0.05 0.55]);
saveas(f4, 'figure_components/f2/pp_DOWN.svg', 'svg');

% Limit cycle
fprintf('Limit Cycle \n');
p.Kp = 1.3;
p.cie = (1 - 0.68) * 35;

[u, e, v1, v2] = compute_fp(@(x) WC([], x, p), [0, 5], [0 5], 0.05, 0.05, 5);
for i = 1:size(u, 1)
    fprintf('Steady State: (%.3f, %.3f), Eigenvalues: (%.3f, %.3f)\n', u(i, 1), u(i, 2), e(i, 1), e(i, 2));
end

X0 = [0; 0];
[~, Y1] = ode45(@WC, 0:0.01:1000, X0, p.opts, p);

f5 = figure(5);
f5.Units = "centimeters";
f5.OuterPosition = [2 10 12 12];
set(gca, 'FontSize', 16, 'FontName', 'Times'); hold on; box on;
xlabel('Excitatory population activity');
ylabel('Inhibitory population activity');
plot(Y1(500:end,1), Y1(500:end,2), 'Color', colour.purple, 'LineWidth', 2);
plot(u(1,1), u(1,2), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'black', 'MarkerSize', 9);
xlim([0 0.5]); ylim([-0.05 0.55]);
saveas(f5, 'figure_components/f2/pp_LC.svg', 'svg');
