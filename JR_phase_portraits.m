clc; clearvars; close all
load("colour.mat")
load("functions/parameters/full_JRpar.mat")
addpath("functions")

%% Phase portrait Limit Cycle
fprintf('Limit Cycle\n');
p.A = 0.334;
p.c4 = 0.3;
p.I = 0.02;

% Create 6D grid points
j.X1 = 0:0.3:1;
j.X2 = 0:0.3:1;
j.X3 = 0:0.3:1;
j.X4 = 0:0.3:1;
j.X5 = 0:0.3:1;
j.X6 = 0:0.3:1;
[j.y1, j.y2, j.y3, j.y4, j.y5, j.y6] = ndgrid(j.X1, j.X2, j.X3, j.X4, j.X5, j.X6);
j.gridpts = [j.y1(:), j.y2(:), j.y3(:), j.y4(:), j.y5(:), j.y6(:)];

% Initialize arrays for fixed points
FixedPTS = NaN(length(j.gridpts),6);
FixedPTS2 = NaN(length(j.gridpts),6);
Exit = NaN(length(j.gridpts),1);

options = optimset('Display', 'off', 'TolFun', 1e-12, 'TolX', 1e-12);

% Find fixed points using fsolve
for i = 1:length(j.gridpts)
    [FixedPTS(i,:), ~, Exit(i)] = fsolve(@(x) full_JP([], x,p), j.gridpts(i,:), options);
    if Exit(i) == 1
        FixedPTS2(i,:) = FixedPTS(i,:);
    end
end

% Remove NaNs and find unique fixed points
FixedPTS3 = FixedPTS2(~any(isnan(FixedPTS2), 2), :);
UniqueEig = unique(round(FixedPTS3, 6), 'rows');

% Calculate Jacobian and eigenvalues at each unique fixed point
for i = 1:size(UniqueEig,1)
    U(i,:) = UniqueEig(i,:);
    A = jacobianest(@(x) full_JP([], x,p), UniqueEig(i,:));
    [v, evalue] = eig(A);
    eval(i,:) = diag(evalue)';
    v1(i,:) = v(:,1)';
    v2(i,:) = v(:,2)';
    v3(i,:) = v(:,3)';
    v4(i,:) = v(:,4)';
    v5(i,:) = v(:,5)';
    v6(i,:) = v(:,6)';
    fprintf('Steady State: (%.3f, %.3f, %.3f, %.3f, %.3f, %.3f), Eigenvalues: (%.3f, %.3f, %.3f, %.3f, %.3f, %.3f)\n', ...
        U(i, :), eval(i, :));
end

% ODE integration setup
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
epsi = 0.005 * ones(1,6);

% Initial condition near second eigenvector direction of first fixed point
X0u2 = U(1,:) + epsi .* v2(1,:);
Tspanf2 = 0:0.01:2000;
[~, Y1] = ode15s(@full_JP, Tspanf2, X0u2, opts, p);
Tspanf2 = 0:0.01:200;
[T, Y2] = ode15s(@full_JP, Tspanf2, Y1(end,:), opts, p);

% Plot phase portrait of Limit Cycle
f1 = figure(1);
f1.Units = "centimeters";
f1.OuterPosition = [2 10 12.2 12.3];
hold on; box on; grid off
set(gca, 'FontSize', 16, 'FontName', 'times');
ylabel('$y_3$', Interpreter='latex')
xlabel('$y_1$', Interpreter='latex')
i = 1; j = 3;

plot(Y1(160000:end, i), Y1(160000:end, j), 'Color', colour.orange, 'LineWidth', 2, 'LineStyle', '-')

plot(U(1,i), U(1,j), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'black', 'MarkerSize', 8, 'LineWidth', 1)
xlim([-0.03 0.4])
ylim([-0.1 1.6])
saveas(f1, 'figure_components/f4/pp_LC.svg', 'svg')

%% Phase portrait SNIC^2
fprintf('SNIC^2\n');
clearvars -except colour

load("functions/parameters/full_JRpar.mat")
addpath("functions")

p.A = 0.334;
p.c4 = 0.2104 + 0.0000905;
p.I = 2.9095e-4 + 0.000001;

j.X1 = -0.2:0.3:1;
j.X2 = -0.20:0.3:1;
j.X3 = -0.20:0.3:1;
j.X4 = 0:0.1:0.2;
j.X5 = 0:0.1:0.2;
j.X6 = 0:0.1:0.2;
[j.y1, j.y2, j.y3, j.y4, j.y5, j.y6] = ndgrid(j.X1, j.X2, j.X3, j.X4, j.X5, j.X6);
j.gridpts = [j.y1(:), j.y2(:), j.y3(:), j.y4(:), j.y5(:), j.y6(:)];

FixedPTS = NaN(length(j.gridpts),6);
FixedPTS2 = NaN(length(j.gridpts),6);
Exit = NaN(length(j.gridpts),1);

options = optimset('Display', 'off', 'TolFun', 1e-9, 'TolX', 1e-10);

% Find fixed points for SNIC^2 parameters
for i = 1:length(j.gridpts)
    [FixedPTS(i,:), ~, Exit(i)] = fsolve(@(x) full_JP([], x,p), j.gridpts(i,:), options);
    if Exit(i) == 1
        FixedPTS2(i,:) = FixedPTS(i,:);
    end
end

FixedPTS3 = FixedPTS2(~any(isnan(FixedPTS2), 2), :);
UniqueEig = unique(round(FixedPTS3, 5), 'rows');

for i = 1:size(UniqueEig,1)
    U(i,:) = UniqueEig(i,:);
    A = jacobianest(@(x) full_JP([], x,p), UniqueEig(i,:));
    [v, evalue] = eig(A);
    eval(i,:) = diag(evalue)';
    v1(i,:) = v(:,1)';
    v2(i,:) = v(:,2)';
    v3(i,:) = v(:,3)';
    v4(i,:) = v(:,4)';
    v5(i,:) = v(:,5)';
    v6(i,:) = v(:,6)';
    fprintf('Steady State: (%.3f, %.3f, %.3f, %.3f, %.3f, %.3f), Eigenvalues: (%.3f, %.3f, %.3f, %.3f, %.3f, %.3f)\n', ...
        U(i,:), eval(i,:));
end

% ODE integration near 3rd fixed point along first eigenvector
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
epsi = 0.005 * ones(1,6);
X0u2 = U(3,:) - epsi .* v1(3,:);
Tspanf2 = 0:0.1:15000;
[~, Y1] = ode15s(@full_JP, Tspanf2, X0u2, opts, p);

[~, Y2] = ode15s(@full_JP, Tspanf2, Y1(end,:), opts, p);

% Plot phase portrait SNIC^2
f2 = figure(2);
f2.Units = "centimeters";
f2.OuterPosition = [2 10 12.2 12.3];
hold on; box on; grid off
set(gca, 'FontSize', 16, 'FontName', 'times');
ylabel('$y_3$', Interpreter='latex')
xlabel('$y_1$', Interpreter='latex')
i = 1; j = 3;

plot(Y1(5600:end,i), Y1(5600:end,j), 'Color', colour.pink, 'LineWidth', 2, 'LineStyle', '-')

plot(U(1,i), U(1,j), '^', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'black', 'MarkerSize', 8, 'LineWidth', 1)
plot(U(2,i), U(2,j), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'black', 'MarkerSize', 8, 'LineWidth', 1)
plot(U(3,i), U(3,j), '^', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'black', 'MarkerSize', 8, 'LineWidth', 1)

xlim([-0.03 0.4])
ylim([-0.1 1.6])
saveas(f2, 'figure_components/f4/pp_SNIC2.svg', 'svg')

%% Phase portrait Up State
fprintf('Up State\n');
clearvars -except colour

load("functions/parameters/full_JRpar.mat")
addpath("functions")

p.A = 0.334;
p.c4 = 0.15;
p.I = 0.02;

j.X1 = -0.2:0.3:1;
j.X2 = -0.20:0.3:1;
j.X3 = -0.20:0.3:1;
j.X4 = 0:0.1:0.2;
j.X5 = 0:0.1:0.2;
j.X6 = 0:0.1:0.2;
[j.y1, j.y2, j.y3, j.y4, j.y5, j.y6] = ndgrid(j.X1, j.X2, j.X3, j.X4, j.X5, j.X6);
j.gridpts = [j.y1(:), j.y2(:), j.y3(:), j.y4(:), j.y5(:), j.y6(:)];

FixedPTS = NaN(length(j.gridpts),6);
FixedPTS2 = NaN(length(j.gridpts),6);
Exit = NaN(length(j.gridpts),1);

options = optimset('Display', 'off', 'TolFun', 1e-9, 'TolX', 1e-10);

% Find fixed points for Up state
for i = 1:length(j.gridpts)
    [FixedPTS(i,:), ~, Exit(i)] = fsolve(@(x) full_JP([], x,p), j.gridpts(i,:), options);
    if Exit(i) == 1
        FixedPTS2(i,:) = FixedPTS(i,:);
    end
end

FixedPTS3 = FixedPTS2(~any(isnan(FixedPTS2), 2), :);
UniqueEig = unique(round(FixedPTS3, 5), 'rows');

for i = 1:size(UniqueEig,1)
    U(i,:) = UniqueEig(i,:);
    A = jacobianest(@(x) full_JP([], x,p), UniqueEig(i,:));
    [v, evalue] = eig(A);
    eval(i,:) = diag(evalue)';
    v1(i,:) = v(:,1)';
    v2(i,:) = v(:,2)';
    v3(i,:) = v(:,3)';
    v4(i,:) = v(:,4)';
    v5(i,:) = v(:,5)';
    v6(i,:) = v(:,6)';
    fprintf('Steady State: (%.3f, %.3f, %.3f, %.3f, %.3f, %.3f), Eigenvalues: (%.3f, %.3f, %.3f, %.3f, %.3f, %.3f)\n', ...
        U(i,:), eval(i,:));
end

% ODE integration near second fixed point along first eigenvector
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
epsi = 0.005 * ones(1,6);
X0u2 = U(2,:) - epsi .* v1(2,:);
Tspanf2 = 0:0.1:15000;
[~, Y1] = ode15s(@full_JP, Tspanf2, X0u2, opts, p);

% Another integration near first fixed point
X0u2 = U(2,:) - epsi .* v1(1,:);
[~, Y2] = ode15s(@full_JP, Tspanf2, X0u2, opts, p);

% Plot phase portrait Up State
f3 = figure(3);
f3.Units = "centimeters";
f3.OuterPosition = [2 10 12.2 12.3];
hold on; box on; grid off
set(gca, 'FontSize', 16, 'FontName', 'times');
ylabel('$y_3$', Interpreter='latex')
xlabel('$y_1$', Interpreter='latex')
i = 1; j = 3;

plot(Y1(:, i), Y1(:, j), 'b-', 'LineWidth', 2)
plot(Y2(:, i), Y2(:, j), 'b-', 'LineWidth', 2)

plot(U(1,i), U(1,j), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'black', 'MarkerSize', 8, 'LineWidth', 1)
plot(U(2,i), U(2,j), '^', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'black', 'MarkerSize', 8, 'LineWidth', 1)
plot(U(3,i), U(3,j), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'black', 'MarkerSize', 8, 'LineWidth', 1)

xlim([-0.03 0.4])
ylim([-0.1 1.6])
saveas(f3, 'figure_components/f4/pp_UP.svg', 'svg')

%% Phase portrait Down State
fprintf('Down State\n');
clearvars -except colour

load("functions/parameters/full_JRpar.mat")
addpath("functions")

p.A = 0.334;
p.c4 = 0.4;
p.I = -0.02;

j.X1 = -0.2:0.3:1;
j.X2 = -0.20:0.3:1;
j.X3 = -0.20:0.3:1;
j.X4 = 0:0.1:0.2;
j.X5 = 0:0.1:0.2;
j.X6 = 0:0.1:0.2;
[j.y1, j.y2, j.y3, j.y4, j.y5, j.y6] = ndgrid(j.X1, j.X2, j.X3, j.X4, j.X5, j.X6);
j.gridpts = [j.y1(:), j.y2(:), j.y3(:), j.y4(:), j.y5(:), j.y6(:)];

FixedPTS = NaN(length(j.gridpts),6);
FixedPTS2 = NaN(length(j.gridpts),6);
Exit = NaN(length(j.gridpts),1);

options = optimset('Display', 'off', 'TolFun', 1e-9, 'TolX', 1e-10);

% Find fixed points for Down state
for i = 1:length(j.gridpts)
    [FixedPTS(i,:), ~, Exit(i)] = fsolve(@(x) full_JP([], x,p), j.gridpts(i,:), options);
    if Exit(i) == 1
        FixedPTS2(i,:) = FixedPTS(i,:);
    end
end

FixedPTS3 = FixedPTS2(~any(isnan(FixedPTS2), 2), :);
UniqueEig = unique(round(FixedPTS3, 5), 'rows');

for i = 1:size(UniqueEig,1)
    U(i,:) = UniqueEig(i,:);
    A = jacobianest(@(x) full_JP([], x,p), UniqueEig(i,:));
    [v, evalue] = eig(A);
    eval(i,:) = diag(evalue)';
    v1(i,:) = v(:,1)';
    v2(i,:) = v(:,2)';
    v3(i,:) = v(:,3)';
    v4(i,:) = v(:,4)';
    v5(i,:) = v(:,5)';
    v6(i,:) = v(:,6)';
    fprintf('Steady State: (%.3f, %.3f, %.3f, %.3f, %.3f, %.3f), Eigenvalues: (%.3f, %.3f, %.3f, %.3f, %.3f, %.3f)\n', ...
        U(i,:), eval(i,:));
end

% ODE integration near fixed points
opts = odeset('RelTol',1e-6,'AbsTol',1e-8);
epsi = 0.005 * ones(1,6);

X0u2 = U(2,:) - epsi .* v1(2,:);
Tspanf2 = 0:0.1:15000;
[~, Y1] = ode15s(@full_JP, Tspanf2, X0u2, opts, p);

X0u2 = U(2,:) - epsi .* v1(1,:);
[~, Y2] = ode15s(@full_JP, Tspanf2, X0u2, opts, p);

% Plot phase portrait Down State
close all
f4 = figure(4);
f4.Units = "centimeters";
f4.OuterPosition = [2 10 12.2 12.3];
hold on; box on; grid off
set(gca, 'FontSize', 16, 'FontName', 'times');
ylabel('$y_3$', Interpreter='latex')
xlabel('$y_1$', Interpreter='latex')
i = 1; j = 3;

plot(Y1(:, i), Y1(:, j), 'b-', 'LineWidth', 2)
plot(Y2(:, i), Y2(:, j), 'b-', 'LineWidth', 2)
plot(U(1,i), U(1,j), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'black', 'MarkerSize', 8, 'LineWidth', 1)
plot(U(2,i), U(2,j), '^', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'black', 'MarkerSize', 8, 'LineWidth', 1)
plot(U(3,i), U(3,j), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'black', 'MarkerSize', 8, 'LineWidth', 1)

xlim([-0.03 0.4])
ylim([-0.1 1.6])
saveas(f4, 'figure_components/f4/pp_DOWN.svg', 'svg')
%% Phase portrait Bistable state
fprintf('Bistable State\n');
clearvars -except colour

load("functions/parameters/full_JRpar.mat")
addpath("functions")
p.A = 0.334;
p.c4 = 0.15;
p.I = -0.02;

% Define 6D grid for initial guesses
j.X1 = -0.2:0.3:1;
j.X2 = -0.20:0.3:1;
j.X3 = -0.20:0.3:1;
j.X4 = 0:0.1:0.2;
j.X5 = 0:0.1:0.2;
j.X6 = 0:0.1:0.2;
[j.y1, j.y2, j.y3, j.y4, j.y5, j.y6] = ndgrid(j.X1, j.X2, j.X3, j.X4, j.X5, j.X6);
j.gridpts = [j.y1(:), j.y2(:), j.y3(:), j.y4(:), j.y5(:), j.y6(:)];

FixedPTS = NaN(length(j.gridpts),6);
FixedPTS2 = NaN(length(j.gridpts),6);
Exit = NaN(length(j.gridpts),1);

options = optimset('Display', 'off', 'TolFun', 1e-9, 'TolX', 1e-10);

% Find fixed points using fsolve
for i = 1:length(j.gridpts)
    [FixedPTS(i,:), ~, Exit(i)] = fsolve(@(x) full_JP([], x, p), j.gridpts(i,:), options);
    if Exit(i) == 1
        FixedPTS2(i,:) = FixedPTS(i,:);
    end
end

% Extract unique fixed points
FixedPTS3 = FixedPTS2(~any(isnan(FixedPTS2), 2), :);
UniqueEig = unique(round(FixedPTS3, 5), 'rows');

% Compute Jacobian, eigenvalues and eigenvectors at fixed points
for i = 1:size(UniqueEig,1)
    U(i,:) = UniqueEig(i,:);
    A = jacobianest(@(x) full_JP([], x, p), UniqueEig(i,:));
    [v, evalue] = eig(A);
    eval(i,:) = diag(evalue)';
    v1(i,:) = v(:,1);
    v2(i,:) = v(:,2);
    v3(i,:) = v(:,3);
    v4(i,:) = v(:,4);
    v5(i,:) = v(:,5);
    v6(i,:) = v(:,6);
    fprintf('Steady State: (%.3f, %.3f, %.3f, %.3f, %.3f, %.3f), Eigenvalues: (%.3f, %.3f, %.3f, %.3f, %.3f, %.3f)\n', ...
        U(i,:), eval(i,:));
end

opts = odeset('RelTol',1e-6, 'AbsTol',1e-8);
epsi=[0.005,0.005,0.05,0.005,0.005,0.005];
% Integrate trajectories near fixed points in positive and negative directions along eigenvector v1
X0u2 = U(2,:)-epsi.*v1(1,:);
Tspanf2 = 0:0.1:15000;
[~, Y1]=ode15s(@full_JP,Tspanf2,X0u2,opts, p);


epsi=[0.005,0.005,0.005,0.005,0.005,0.005];
X0u2 = U(2,:)+epsi.*v1(1,:);
Tspanf2 = 0:0.1:15000;
[~, Y2]=ode15s(@full_JP,Tspanf2,X0u2,opts, p);


X0u2 = U(3,:)+epsi.*v1(1,:);
Tspanf2 = 0:0.1:15000;
[~, Y21]=ode15s(@full_JP,Tspanf2,X0u2,opts, p);


X0u2 = U(3,:)-epsi.*v1(1,:);
Tspanf2 = 0:0.1:15000;
[~, Y22]=ode15s(@full_JP,Tspanf2,X0u2,opts, p);

X0u2 = U(4,:)+epsi.*v1(1,:);
Tspanf2 = 0:0.1:15000;
[~, Y3]=ode15s(@full_JP,Tspanf2,X0u2,opts, p);

X0u2 = U(4,:)-epsi.*v1(1,:);
Tspanf2 = 0:0.1:15000;
[~, Y4]=ode15s(@full_JP,Tspanf2,X0u2,opts, p);

% Plot phase portrait of Bistable state
f5 = figure(5);
f5.Units = "centimeters";
f5.OuterPosition = [2 10 12.2 12.3];
hold on; box on; grid off
set(gca, 'FontSize', 16, 'FontName', 'times');
ylabel('$y_3$', Interpreter='latex')
xlabel('$y_1$', Interpreter='latex')
i = 1; j = 3;

% Plot separatrices (trajectories near saddles)
plot(Y1(:, i), Y1(:, j), 'b-', 'LineWidth', 2)
plot(Y2(:, i), Y2(:, j), 'b-', 'LineWidth', 2)
% Uncomment below if you want to plot additional trajectories
% plot(Y21(:, i), Y21(:, j), 'b-', 'LineWidth', 2)
% plot(Y22(:, i), Y22(:, j), 'b--', 'LineWidth', 2)
plot(Y3(:, i), Y3(:, j), 'b-', 'LineWidth', 2)
plot(Y4(:, i), Y4(:, j), 'b-', 'LineWidth', 2)

% Plot fixed points with different markers
plot(U(1,i), U(1,j), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'black', 'MarkerSize', 8, 'LineWidth', 1)
plot(U(2,i), U(2,j), '^', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'black', 'MarkerSize', 8, 'LineWidth', 1)
plot(U(3,i), U(3,j), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'black', 'MarkerSize', 8, 'LineWidth', 1)
plot(U(5,i), U(5,j), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'black', 'MarkerSize', 8, 'LineWidth', 1)
plot(U(4,i), U(4,j), '^', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'black', 'MarkerSize', 8, 'LineWidth', 1)

xlim([-0.03 0.4])
ylim([-0.1 1.6])
saveas(f5, 'figure_components/f4/pp_bi.svg', 'svg')