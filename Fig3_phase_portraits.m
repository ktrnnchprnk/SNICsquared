clc; clearvars
load("colour.mat")

load("functions/parameters/TM_par.mat")
addpath("functions")

%% Phase portrait for SNIC^2
fprintf('SNIC^2\n');
% define the parameters
p.I=1.53562;
p.omega=26.8+0.005;
% find equilibria, associated eigenvalues and eigenvectors
[u,e,v1,v2] = compute_fp(@(x) TM([], x,p),[0, 22],[0.2, 0.6],0.5,0.05,4);
for i = 1:size(u, 1)
    fprintf('Steady State: (%.3f, %.3f), Eigenvalues: (%.3f, %.3f)\n', ...
        u(i, 1), u(i, 2), ...
        e(i, 1), e(i, 2));
end
% Find the heteroclinic loop
epsi=[0.005,0.005];
X0u2 = u(1,:)+epsi.*v1(1,:);
Tspanf2 = 0:0.0001:74;
[~, Y1]=ode45(@TM,Tspanf2,X0u2,p.opts, p);
close all
f1=figure(1);
f1.Units="centimeters";
f1.OuterPosition = [2 10 12.2 12.3];
set ( gca , 'FontSize' , 16 , 'fontname' , 'times');
hold on; box on; grid off
xlabel('Mean synaptic depression', Interpreter='latex')
ylabel('Mean activity', Interpreter='latex')



plot(Y1(1:end,2),Y1(1:end,1), 'Color',colour.pink, 'LineWidth',2, 'LineStyle','-')



n=1;
plot(u(n,2),u(n,1),'^','MarkerFaceColor', ...
'k','MarkerEdgeColor','black', 'MarkerSize', 9,'LineWidth',1)

n=2;
plot(u(n,2),u(n,1),'o','MarkerFaceColor', ...
'white','MarkerEdgeColor','black', 'MarkerSize', 8,'LineWidth',1)
n=3;
plot(u(n,2),u(n,1),'^','MarkerFaceColor', ...
'black','MarkerEdgeColor','black', 'MarkerSize', 9,'LineWidth',1)
ylim([3 33])
xlim([0.2, 1.1])


saveas(f1,'figure_components/f3/pp_SNIC2.svg', 'svg')
%% Phase portrait for bistability 
fprintf('Bistability\n');
p.I=1.2;
p.omega=28;

[u,e,v1,v2] = compute_fp(@(x) TM([], x,p),[0, 22],[0.2, 0.6],0.5,0.05,4);
for i = 1:size(u, 1)
    fprintf('Steady State: (%.3f, %.3f), Eigenvalues: (%.3f, %.3f)\n', ...
        u(i, 1), u(i, 2), ...
        e(i, 1), e(i, 2));
end
% Unstable manifold 
epsi=[0.005,0.005];
X0u2 = u(2,:)+epsi.*v1(2,:);
Tspanf2 = 0:0.0001:74;
[~, Y1]=ode45(@TM,Tspanf2,X0u2,p.opts, p);

X0u2 = u(2,:)-epsi.*v1(2,:);
Tspanf2 = 0:0.0001:74;
[~, Y2]=ode45(@TM,Tspanf2,X0u2,p.opts, p);

X0u2 = u(4,:)-epsi.*v1(4,:);
Tspanf2 = 0:0.0001:10;
[~, Y7]=ode45(@TM,Tspanf2,X0u2,p.opts, p);

X0u2 = u(4,:)+epsi.*v2(4,:);
Tspanf2 = 0:0.0001:10;
[~, Y8]=ode45(@TM,Tspanf2,X0u2,p.opts, p);

% Stable manifold 

X0u2 = u(2,:)-epsi.*v1(2,:);
Tspanf2 = 10:-0.0001:0;
[~, Y3]=ode45(@TM,Tspanf2,X0u2,p.opts, p);

epsi=[0.05,0.05];
X0u2 = u(2,:)-epsi.*v2(2,:);
Tspanf2 = 10:-0.0001:0;
[~, Y4]=ode45(@TM,Tspanf2,X0u2,p.opts, p);

epsi=[0.05,0.05];
X0u2 = u(4,:)-epsi.*v2(4,:);
Tspanf2 = 10:-0.0001:0;
[~, Y5]=ode45(@TM,Tspanf2,X0u2,p.opts, p);

epsi=[0.05,0.05];
X0u2 = u(4,:)-epsi.*v1(4,:);
Tspanf2 = 10:-0.0001:0;
[~, Y6]=ode45(@TM,Tspanf2,X0u2,p.opts, p);

f2=figure(2);
f2.Units="centimeters";
f2.OuterPosition = [2 10 12.2 12.3];
set ( gca , 'FontSize' , 16 , 'fontname' , 'times');
hold on; box on; grid off
xlabel('Mean synaptic depression', Interpreter='latex')
ylabel('Mean activity', Interpreter='latex')




plot(Y1(1:end,2),Y1(1:end,1), 'Color','b', 'LineWidth',2, 'LineStyle','-')
plot(Y2(1:end,2),Y2(1:end,1), 'Color','b', 'LineWidth',2, 'LineStyle','-')

plot(Y3(1:end,2),Y3(1:end,1), 'Color','r', 'LineWidth',2, 'LineStyle','-')
plot(Y4(1:end,2),Y4(1:end,1), 'Color','r', 'LineWidth',2, 'LineStyle','-')

plot(Y5(1:end,2),Y5(1:end,1), 'Color','r', 'LineWidth',2, 'LineStyle','-')
plot(Y6(1:end,2),Y6(1:end,1), 'Color','r', 'LineWidth',2, 'LineStyle','-')

plot(Y7(1:end,2),Y7(1:end,1), 'Color','b', 'LineWidth',2, 'LineStyle','-')
plot(Y8(1:end,2),Y8(1:end,1), 'Color','b', 'LineWidth',2, 'LineStyle','-')

n=1;
plot(u(n,2),u(n,1),'o','MarkerFaceColor', ...
'k','MarkerEdgeColor','black', 'MarkerSize', 9,'LineWidth',1)

n=2;
plot(u(n,2),u(n,1),'^','MarkerFaceColor', ...
'white','MarkerEdgeColor','black', 'MarkerSize', 8,'LineWidth',1)
n=3;
plot(u(n,2),u(n,1),'o','MarkerFaceColor', ...
'w','MarkerEdgeColor','black', 'MarkerSize', 9,'LineWidth',1)

n=4;
plot(u(n,2),u(n,1),'^','MarkerFaceColor', ...
'white','MarkerEdgeColor','black', 'MarkerSize', 8,'LineWidth',1)
n=5;
plot(u(n,2),u(n,1),'o','MarkerFaceColor', ...
'k','MarkerEdgeColor','black', 'MarkerSize', 9,'LineWidth',1)
ylim([3 33])
xlim([0.2, 1.1])

saveas(f2,'figure_components/f3/pp_Bi.svg', 'svg')

%% Phase Portrait DOWN
fprintf('Down State\n');
p.I=1.35;
p.omega=25;

[u,e,v1,v2] = compute_fp(@(x) TM([], x,p),[0, 22],[0.2, 0.6],0.5,0.05,4);
for i = 1:size(u, 1)
    fprintf('Steady State: (%.3f, %.3f), Eigenvalues: (%.3f, %.3f)\n', ...
        u(i, 1), u(i, 2), ...
        e(i, 1), e(i, 2));
end

% Unstable manifold
epsi=[0.005,0.005];
X0u2 = u(2,:)+epsi.*v1(2,:);
Tspanf2 = 0:0.0001:74;
[~, Y1]=ode45(@TM,Tspanf2,X0u2,p.opts, p);

X0u2 = u(2,:)-epsi.*v1(2,:);
Tspanf2 = 0:0.0001:74;
[~, Y2]=ode45(@TM,Tspanf2,X0u2,p.opts, p);

% Stable manifold
X0u2 = u(2,:)-epsi.*v1(2,:);
Tspanf2 = 10:-0.0001:0;
[~, Y3]=ode45(@TM,Tspanf2,X0u2,p.opts, p);

X0u2 = u(2,:)+epsi.*v2(2,:);
Tspanf2 = 10:-0.0001:0;
[~, Y4]=ode45(@TM,Tspanf2,X0u2,p.opts, p);

f3=figure(3);
f3.Units="centimeters";
f3.OuterPosition = [2 10 12.2 12.3];
set ( gca , 'FontSize' , 16 , 'fontname' , 'times');
hold on; box on; grid off
xlabel('Mean synaptic depression', Interpreter='latex')
ylabel('Mean activity', Interpreter='latex')




plot(Y1(1:end,2),Y1(1:end,1), 'Color','b', 'LineWidth',2, 'LineStyle','-')
plot(Y2(1:end,2),Y2(1:end,1), 'Color','b', 'LineWidth',2, 'LineStyle','-')

plot(Y3(1:end,2),Y3(1:end,1), 'Color','r', 'LineWidth',2, 'LineStyle','-')
plot(Y4(1:end,2),Y4(1:end,1), 'Color','r', 'LineWidth',2, 'LineStyle','-')
n=1;
plot(u(n,2),u(n,1),'o','MarkerFaceColor', ...
'k','MarkerEdgeColor','black', 'MarkerSize', 9,'LineWidth',1)

n=2;
plot(u(n,2),u(n,1),'^','MarkerFaceColor', ...
'white','MarkerEdgeColor','black', 'MarkerSize', 8,'LineWidth',1)
n=3;
plot(u(n,2),u(n,1),'o','MarkerFaceColor', ...
'w','MarkerEdgeColor','black', 'MarkerSize', 9,'LineWidth',1)
ylim([3 33])
xlim([0.2, 1.1])

saveas(f3,'figure_components/f3/pp_down.svg', 'svg')
%% Phase portrait UP
fprintf('Up State\n');
p.I=1.5;
p.omega=29;

[u,e,v1,v2] = compute_fp(@(x) TM([], x,p),[0, 22],[0.2, 0.6],0.5,0.05,4);
for i = 1:size(u, 1)
    fprintf('Steady State: (%.3f, %.3f), Eigenvalues: (%.3f, %.3f)\n', ...
        u(i, 1), u(i, 2), ...
        e(i, 1), e(i, 2));
end

% Unstable manifold 
epsi=[0.005,0.005];
X0u2 = u(2,:)+epsi.*v1(2,:);
Tspanf2 = 0:0.0001:74;
[~, Y1]=ode45(@TM,Tspanf2,X0u2,p.opts, p);

X0u2 = u(2,:)-epsi.*v1(2,:);
Tspanf2 = 0:0.0001:74;
[~, Y2]=ode45(@TM,Tspanf2,X0u2,p.opts, p);

% Stable manifold
X0u2 = u(2,:)-epsi.*v1(2,:);
Tspanf2 = 10:-0.0001:0;
[~, Y3]=ode45(@TM,Tspanf2,X0u2,p.opts, p);

X0u2 = u(2,:)-epsi.*v2(2,:);
Tspanf2 = 10:-0.0001:0;
[~, Y4]=ode45(@TM,Tspanf2,X0u2,p.opts, p);

close all
f4=figure(4);
f4.Units="centimeters";
f4.OuterPosition = [2 10 12.2 12.3];
set ( gca , 'FontSize' , 16 , 'fontname' , 'times');
hold on; box on; grid off
xlabel('Mean synaptic depression', Interpreter='latex')
ylabel('Mean activity', Interpreter='latex')


% contour(X1, X2, nullcline_x1, [0 0], 'b', 'LineWidth', 2); 
% contour(X1, X2, nullcline_x2, [0 0], 'r', 'LineWidth', 2);

plot(Y1(1:end,2),Y1(1:end,1), 'Color','b', 'LineWidth',2, 'LineStyle','-')
plot(Y2(1:end,2),Y2(1:end,1), 'Color','b', 'LineWidth',2, 'LineStyle','-')

plot(Y3(1:end,2),Y3(1:end,1), 'Color','r', 'LineWidth',2, 'LineStyle','-')
plot(Y4(1:end,2),Y4(1:end,1), 'Color','r', 'LineWidth',2, 'LineStyle','-')
% plot(u2(1,1),u2(1,2),'o','MarkerFaceColor', ...
% 'black','MarkerEdgeColor','black', 'MarkerSize', 8,'LineWidth',1)
n=1;
plot(u(n,2),u(n,1),'o','MarkerFaceColor', ...
'w','MarkerEdgeColor','black', 'MarkerSize', 9,'LineWidth',1)

n=2;
plot(u(n,2),u(n,1),'^','MarkerFaceColor', ...
'white','MarkerEdgeColor','black', 'MarkerSize', 8,'LineWidth',1)
n=3;
plot(u(n,2),u(n,1),'o','MarkerFaceColor', ...
'k','MarkerEdgeColor','black', 'MarkerSize', 9,'LineWidth',1)
ylim([3 33])
xlim([0.2, 1.1])


saveas(f4,'figure_components/f3/pp_up.svg', 'svg')
%% Phase Portrait Limit Cycle
fprintf('Limit Cycle\n');
p.I=1.8;
p.omega=25.9;

[u,e,v1,v2] = compute_fp(@(x) TM([], x,p),[0, 22],[0.2, 0.6],0.5,0.05,4);
for i = 1:size(u, 1)
    fprintf('Steady State: (%.3f, %.3f), Eigenvalues: (%.3f, %.3f)\n', ...
        u(i, 1), u(i, 2), ...
        e(i, 1), e(i, 2));
end


X0u2=[0;0];
Tspanf2 = 0:0.001:150;
[~, Y1]=ode45(@TM,Tspanf2,X0u2,p.opts, p);

f5=figure(5);
f5.Units="centimeters";
f5.OuterPosition = [2 10 12.2 12.3];
set ( gca , 'FontSize' , 16 , 'fontname' , 'times');
hold on; box on; grid off
xlabel('Mean synaptic depression', Interpreter='latex')
ylabel('Mean activity', Interpreter='latex')




plot(Y1(5000:end,2),Y1(5000:end,1), 'Color',colour.orange, 'LineWidth',2, 'LineStyle','-')

n=1;
plot(u(n,2),u(n,1),'o','MarkerFaceColor', ...
'w','MarkerEdgeColor','black', 'MarkerSize', 9,'LineWidth',1)

ylim([3 33])
xlim([0.2, 1.1])

saveas(f5,'figure_components/f3/pp_LC.svg', 'svg')