clc; clearvars; close all

% Load color scheme and model parameters
load("colour.mat")
load("functions/parameters/ML_par.mat")
addpath("functions")

%% Phase portrait for SNIC^2
fprintf('SNIC^2\n');

% Parameter setup for SNIC^2
p.VCa = 127;
p.V1 = -1.23;
p.I   = 0.079;
p.gK  = 11.9639;

% Compute fixed points and their stability
[u, e, v1, v2] = compute_fp(@(x) ML([], x, p), [-60, 25], [-0.1, 1.1], 0.5, 0.05, 4);

% Print fixed points and corresponding eigenvalues
for i = 1:size(u, 1)
    fprintf('Steady State: (%.3f, %.3f), Eigenvalues: (%.3f, %.3f)\n', ...
        u(i, 1), u(i, 2), e(i, 1), e(i, 2));
end

%% Simulate separatrix from saddle node
epsi   = [0.005, 0.005];                     % Small offset to perturb initial condition
X0u2   = u(3,:) - epsi .* v1(3,:);           % Initial condition slightly off saddle
Tspan  = 0:0.005:2000;                       % Time span for simulation
[~, Y] = ode45(@ML, Tspan, X0u2, p.opts, p); % Integrate system

%% Plot phase portrait
f1 = figure(1);
f1.Units = "centimeters";
f1.OuterPosition = [2 10 12.2 12.3];
set(gca, 'FontSize', 16, 'FontName', 'Times');
hold on; box on;

xlabel('Membrane potential')
ylabel('Gating variable')

% Plot trajectory
plot(Y(:,1), Y(:,2), 'Color', colour.pink, 'LineWidth', 2)

% Plot fixed points: stable nodes (filled triangle), saddle (circle), and unstable (triangle)
plot(u(1,1), u(1,2), '^', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 9, 'LineWidth', 1)
plot(u(2,1), u(2,2), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'LineWidth', 1)
plot(u(3,1), u(3,2), '^', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'LineWidth', 1)

% Set plot limits
xlim([-65 25])
ylim([-0.05 1.05])

% Save figure
saveas(f1, 'figure_components/supp_2/pp_SNIC2.svg', 'svg')
%% Phase Portrait Bistable
fprintf('Bistability\n');

% Parameter setup for bistability
p.VCa = 127;
p.V1  = -1.23;
p.I   = -1;
p.gK  = 11.5;

% Compute fixed points and their eigenvalues
[u, e, v1, v2] = compute_fp(@(x) ML([], x, p), [-60, 30], [-0.1, 1.1], 0.5, 0.05, 4);

% Display fixed point information
for i = 1:size(u, 1)
    fprintf('Steady State: (%.3f, %.3f), Eigenvalues: (%.3f, %.3f)\n', ...
        u(i, 1), u(i, 2), e(i, 1), e(i, 2));
end

%% Compute separatrices around saddle nodes

% Small perturbation for stable/unstable manifold direction
epsi = [0.005, 0.005];

% Trajectories from fixed point 2 (saddle)
[~, Y1] = ode45(@ML, 0:0.005:2000, u(2,:) - epsi .* v1(2,:), p.opts, p);
[~, Y2] = ode45(@ML, 0:0.005:2000, u(2,:) + epsi .* v1(2,:), p.opts, p);
[~, Y3] = ode45(@ML, 10:-0.005:0, u(2,:) + epsi .* v2(2,:), p.opts, p);
[~, Y4] = ode45(@ML, 0.3:-0.005:0, u(2,:) - epsi .* v2(2,:), p.opts, p);

% Trajectories from fixed point 4 (saddle)
[~, Y5] = ode45(@ML, 0:0.005:2000, u(4,:) - epsi .* v1(4,:), p.opts, p);
[~, Y6] = ode45(@ML, 0:0.005:2000, u(4,:) + epsi .* v1(4,:), p.opts, p);
[~, Y7] = ode45(@ML, 10:-0.005:0, u(4,:) + epsi .* v2(4,:), p.opts, p);
[~, Y8] = ode45(@ML, 1.5:-0.005:0, u(4,:) - epsi .* v2(4,:), p.opts, p);

%% Plot phase portrait
f2 = figure(2);
f2.Units = "centimeters";
f2.OuterPosition = [2 10 12.2 12.3];
set(gca, 'FontSize', 16, 'FontName', 'Times');
hold on; box on;

xlabel('Membrane potential')
ylabel('Gating variable')

% Plot stable manifolds (blue) and unstable manifolds (red)
plot(Y1(:,1), Y1(:,2), 'b', 'LineWidth', 2)
plot(Y2(:,1), Y2(:,2), 'b', 'LineWidth', 2)
plot(Y3(:,1), Y3(:,2), 'r', 'LineWidth', 2)
plot(Y4(:,1), Y4(:,2), 'r', 'LineWidth', 2)

plot(Y5(:,1), Y5(:,2), 'b', 'LineWidth', 2)
plot(Y6(:,1), Y6(:,2), 'b', 'LineWidth', 2)
plot(Y7(:,1), Y7(:,2), 'r', 'LineWidth', 2)
plot(Y8(:,1), Y8(:,2), 'r', 'LineWidth', 2)

% Plot fixed points
plot(u(1,1), u(1,2), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 9, 'LineWidth', 1)  % stable node
plot(u(2,1), u(2,2), '^', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'LineWidth', 1) % saddle
plot(u(3,1), u(3,2), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'LineWidth', 1) % saddle/focus
plot(u(4,1), u(4,2), '^', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'LineWidth', 1) % saddle
plot(u(5,1), u(5,2), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 9, 'LineWidth', 1)  % stable node

% Set axis limits
xlim([-60 25])
ylim([-0.05 1.05])

% Save figure
saveas(f2, 'figure_components/supp_2/pp_bistability.svg', 'svg')


%% Phase Portrait Down
fprintf('Down State\n');

% Parameter setup for Down state
p.VCa = 127;
p.V1  = -1.23;
p.I   = -5;
p.gK  = 16;

% Compute fixed points and their stability
[u, e, v1, v2] = compute_fp(@(x) ML([], x, p), [-60, 20], [-0.1, 1.1], 0.5, 0.05, 4);

% Display fixed points and their eigenvalues
for i = 1:size(u, 1)
    fprintf('Steady State: (%.3f, %.3f), Eigenvalues: (%.3f, %.3f)\n', ...
        u(i, 1), u(i, 2), e(i, 1), e(i, 2));
end

% Integrate stable and unstable manifolds from saddle (point 2)
epsi = [0.005, 0.005];
[~, Y1] = ode45(@ML, 0:0.005:2000, u(2,:) - epsi .* v1(2,:), p.opts, p);
[~, Y2] = ode45(@ML, 0:0.005:2000, u(2,:) + epsi .* v1(2,:), p.opts, p);
[~, Y3] = ode45(@ML, 10:-0.005:0,    u(2,:) + epsi .* v2(2,:), p.opts, p);
[~, Y4] = ode45(@ML, 0.6:-0.005:0,   u(2,:) - epsi .* v2(2,:), p.opts, p);

% Plot Down state phase portrait
f3 = figure(3);
f3.Units = "centimeters";
f3.OuterPosition = [2 10 12.2 12.3];
set(gca, 'FontSize', 16, 'FontName', 'Times');
hold on; box on;

xlabel('Membrane potential')
ylabel('Gating variable')

plot(Y1(:,1), Y1(:,2), 'b', 'LineWidth', 2)
plot(Y2(:,1), Y2(:,2), 'b', 'LineWidth', 2)
plot(Y3(:,1), Y3(:,2), 'r', 'LineWidth', 2)
plot(Y4(:,1), Y4(:,2), 'r', 'LineWidth', 2)

% Plot fixed points
plot(u(1,1), u(1,2), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 9, 'LineWidth', 1)  % stable
plot(u(2,1), u(2,2), '^', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'LineWidth', 1)  % saddle
plot(u(3,1), u(3,2), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'LineWidth', 1)  % unstable

xlim([-65 25])
ylim([-0.05 1.05])

saveas(f3, 'figure_components/supp_2/pp_Down.svg', 'svg')

%% Phase Portrait Up
fprintf('Up State\n');

% Parameter setup for Up state
p.VCa = 127;
p.V1  = -1.23;
p.I   = 0.2;
p.gK  = 11.4;

% Compute fixed points and their stability
[u, e, v1, v2] = compute_fp(@(x) ML([], x, p), [-60, 20], [-0.1, 1.1], 0.5, 0.05, 4);

% Display fixed points and their eigenvalues
for i = 1:size(u, 1)
    fprintf('Steady State: (%.3f, %.3f), Eigenvalues: (%.3f, %.3f)\n', ...
        u(i, 1), u(i, 2), e(i, 1), e(i, 2));
end

% Integrate stable and unstable manifolds from saddle (point 2)
epsi = [0.005, 0.005];
[~, Y1] = ode45(@ML, 0:0.005:2000, u(2,:) - epsi .* v1(2,:), p.opts, p);
[~, Y2] = ode45(@ML, 0:0.005:2000, u(2,:) + epsi .* v1(2,:), p.opts, p);
[~, Y3] = ode45(@ML, 2:-0.005:0,     u(2,:) + epsi .* v2(2,:), p.opts, p);
[~, Y4] = ode45(@ML, 10:-0.005:0,    u(2,:) - epsi .* v2(2,:), p.opts, p);

% Plot Up state phase portrait
f4 = figure(4);
f4.Units = "centimeters";
f4.OuterPosition = [2 10 12.2 12.3];
set(gca, 'FontSize', 16, 'FontName', 'Times');
hold on; box on;

xlabel('Membrane potential')
ylabel('Gating variable')

plot(Y1(:,1), Y1(:,2), 'b', 'LineWidth', 2)
plot(Y2(:,1), Y2(:,2), 'b', 'LineWidth', 2)
plot(Y3(:,1), Y3(:,2), 'r', 'LineWidth', 2)
plot(Y4(:,1), Y4(:,2), 'r', 'LineWidth', 2)

% Plot fixed points
plot(u(1,1), u(1,2), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'LineWidth', 1)  % unstable
plot(u(2,1), u(2,2), '^', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 8, 'LineWidth', 1)  % saddle
plot(u(3,1), u(3,2), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 9, 'LineWidth', 1)  % stable

xlim([-65 25])
ylim([-0.05 1.05])

saveas(f4, 'figure_components/supp_2/pp_Up.svg', 'svg')


%% Phase Portrait Limit Cycle
fprintf('Limit Cycle\n');

% Parameter setup to produce a stable limit cycle
p.VCa = 127;
p.V1  = -1.23;
p.I   = 50;
p.gK  = 20;

% Compute fixed points and their stability
[u, e, v1, v2] = compute_fp(@(x) ML([], x, p), [-70, -30], [0, 1.2], 0.5, 0.05, 4);

% Display fixed point information
for i = 1:size(u, 1)
    fprintf('Steady State: (%.3f, %.3f), Eigenvalues: (%.3f, %.3f)\n', ...
        u(i, 1), u(i, 2), e(i, 1), e(i, 2));
end

% Simulate trajectory from arbitrary initial condition to trace the limit cycle
X0     = [0, 0];                     % Initial condition
Tspan  = 0:0.001:200;                % Long time span to reach limit cycle
[~, Y] = ode45(@ML, Tspan, X0, p.opts, p);

% Plot Limit Cycle phase portrait
f5 = figure(5);
f5.Units = "centimeters";
f5.OuterPosition = [2 10 12.2 12.3];
set(gca, 'FontSize', 16, 'FontName', 'Times');
hold on; box on;

xlabel('Membrane potential')
ylabel('Gating variable')

% Plot last part of trajectory (assumed to be on the limit cycle)
plot(Y(end-1000:end, 1), Y(end-1000:end, 2), ...
     'Color', colour.orange, 'LineWidth', 2)

% Plot one of the fixed points (e.g., saddle or unstable focus)
plot(u(1,1), u(1,2), 'o', 'MarkerFaceColor', 'w', ...
     'MarkerEdgeColor', 'k', 'MarkerSize', 9, 'LineWidth', 1)

xlim([-65 25])
ylim([-0.05 1.05])

% Save figure
saveas(f5, 'figure_components/supp_2/pp_LC.svg', 'svg')



