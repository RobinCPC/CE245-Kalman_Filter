%% CMPE 245 HW1 b
% -------------------------
% CMPE 245 - Spring, 2016
% Chien-Pin Chen
% 04/05/2016
% MATLAB script for HW1 b
% --------------------------

% The original equation.
% mLx'' = -mgsin(x)
% x'' = -(g/L)sin(x) + f(t)
% rewrite :
% x1' = x2
% x2' = x1'' = -(g/L)*sin(x1) + f(t)  , where (g/L)=0.1

% clear workspace
clear all;

% Set time interval and parameter of 
T = 500;    % Total time interval
la = 1/0.02;  % parameter of exponential distribution

% Set ode function
fun = @(t,x) [x(2); -0.1*sin(x(1))];

% get # of events happenned by p(t_k - t_k-1)
t_p = (0);  % array to when external force hitting
k = 1;
while (t_p(end) < T)
    k = k + 1;
    t_p(end + 1) = t_p(end) + random('exp', la);
end
n_event = k;
t_p(end) = T;

% Simulate numerical equation.
x0 = [0 0];     % Initial Condition x1=0 x2=0
t=[];
x=[];

for k=1:(n_event-1)
    [ts, xs] = ode45(fun, [t_p(k) t_p(k+1)], x0);
    t = cat(1, t, ts);
    x = cat(1, x, xs);
    ts=0;xs=0;
    f = random('norm', 0, 0.2);
    x0 = x(end) + [0 f];
end

% Plot the results
h1 = figure;
axes('Parent', h1, 'FontWeight', 'bold', 'FontSize', 14);
plot(t, x(:,1), 'r', 'LineWidth',2);
title('Angular Position X_1(t)  v.s.  Time t','FontWeight','bold', 'FontSize',12);
xlabel('Time t, unit: sec', 'FontWeight','bold', 'FontSize',12);
ylabel('Angular Postion X_1(t), unit: rad', 'FontWeight','bold', 'FontSize',12);

h2 = figure;
axes('Parent', h2, 'FontWeight', 'bold', 'FontSize', 14);
plot(t, x(:,2), 'k', 'LineWidth',2);
title('Angular Velocity X_2(t)  v.s.  Time t','FontWeight','bold', 'FontSize',12);
xlabel('Time t, unit: sec', 'FontWeight','bold', 'FontSize',12);
ylabel('Angular Velocity X_2(t), unit: rad', 'FontWeight','bold', 'FontSize',12);

h3 = figure;
axes('Parent', h3,'FontWeight','bold', 'FontSize',14);
plot(x(:,1), x(:,2),'LineWidth',2);
axis equal
title('Angular Position X_1(t)  v.s.  Angular Velovity X_2(t)','FontWeight','bold', 'FontSize',12);
xlabel('Angular Postion X_1(t), unit: rad', 'FontWeight','bold', 'FontSize',12);
ylabel('Angular Velocity X_2(t), unit: rad', 'FontWeight','bold', 'FontSize',12);
