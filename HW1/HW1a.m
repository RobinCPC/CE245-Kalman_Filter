%% CMPE 245 HW1 a
% -------------------------
% CMPE 245 - Spring, 2016
% Chien-Pin Chen
% 04/05/2016
% MATLAB script for HW1 a
% --------------------------

% The original equation.
% mLx'' = -mgsin(x)
% x'' = -(g/L)sin(x)
% rewrite :
% x1' = x2
% x2' = x1'' = -(g/L)*sin(x1)   , where (g/L)=0.1
% I.C.: x1 = pi/6;   x2 = 0

% clear workspace
clear all;

% Simulate numerical equation.
T = 600;
fun = @(t,x) [x(2); -0.1*sin(x(1))];
[t,x] = ode23(fun, [0 T], [pi/6 0]);

% Plot the results
h1 = figure;
axes('Parent', h1, 'FontWeight', 'bold', 'FontSize', 12);
plot(t, x(:,1),'r');
title('Angular Position X_1(t)  v.s.  Time t','FontWeight','bold', 'FontSize',12);
xlabel('Time t, unit: sec', 'FontWeight','bold', 'FontSize',12);
ylabel('Angular Postion X_1(t), unit: rad', 'FontWeight','bold', 'FontSize',12);

h2 = figure;
axes('Parent', h2, 'FontWeight', 'bold', 'FontSize', 12);
plot(t, x(:,2),'k');
title('Angular Velocity X_2(t)  v.s.  Time t','FontWeight','bold', 'FontSize',12);
xlabel('Time t, unit: sec', 'FontWeight','bold', 'FontSize',12);
ylabel('Angular Velocity X_2(t), unit: rad', 'FontWeight','bold', 'FontSize',12);

h3 = figure;
axes('Parent', h3, 'FontWeight', 'bold', 'FontSize', 12);
plot(x(:,1), x(:,2));
axis equal
title('Angular Position X_1(t)  v.s.  Angular Velovity X_2(t)','FontWeight','bold', 'FontSize',12);
xlabel('Angular Postion X_1(t), unit: rad', 'FontWeight','bold', 'FontSize',12);
ylabel('Angular Velocity X_2(t), unit: rad', 'FontWeight','bold', 'FontSize',12);
%plot(t, x(:,1),'r-', t, x(:,2),'b', t,x(:,2));