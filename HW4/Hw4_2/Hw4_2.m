%% HW4-2
% -------------------------
% CMPE 245 - Spring, 2016
% Provided by Prof. Dejan Milutinovic
% Edited by Chien-Pin Chen
% 05/08/2016
% MATLAB script for HW4-2
% --------------------------
clear all;

% setting initial value
vars0 = 100;
varv0 = 100;    % steady state value = 400/4 =100
covsv = 0;      % because s and v are independant for initial condition.

% Simulate numerical equation.
[t solX]=ode45(@dCovXfun,[0 5],[vars0; covsv; varv0]);

% plot the evolution of the variacne of s(t) and v(t)
P11=solX(:,1);  % extract the variance of s(t)
P12=solX(:,2);  % extract the covariance of [s(t) v(t)]
P22=solX(:,3);  % extract the variance of v(t)

% initalize figure and its handles for ploting line later
f1 = figure;        % delcare an empty canvas 
h_f1 = gca;         % get handle value (pointer) of above canvs
% set(h_f1, 'fontsize', 16, 'FontName','Times');
set(h_f1,'NextPlot', 'add');
xlabel('Time [s]','FontSize',16,'FontName','Times');
ylabel('var \{s(t)\} = P_{11}(t)','FontSize',16,'FontName','Times');

plot(h_f1, t,P11, 'LineWidth',2);
plot(h_f1, t,2*vars0*ones(length(t)),'r-');   %  find at what time var(s(t*) = 2*var(s(0))



% initalize figure and its handles for ploting line later
f2 = figure;        % delcare an empty canvas 
h_f2 = gca;         % get handle value (pointer) of above canvs
% set(h_f2, 'fontsize', 16, 'FontName','Times');
set(h_f2,'NextPlot', 'add');
xlabel('Time [s]','FontSize',16,'FontName','Times');
ylabel('var \{v(t)\} = P_{22}(t)','FontSize',16,'FontName','Times');

plot(h_f2, t,P22, 'LineWidth',2);


% initalize figure and its handles for ploting line later
f3 = figure;        % delcare an empty canvas 
h_f3 = gca;         % get handle value (pointer) of above canvs
%set(h_f3, 'fontsize', 16, 'FontName','Times');
set(h_f3,'NextPlot', 'add');
xlabel('Time [s]','FontSize',16,'FontName','Times');
ylabel('cov \{s(t)v(t)\} = P_{12}(t) = P_{21}(t)','FontSize',16,'FontName','Times');

plot(t,P12, 'LineWidth',2)
ylim([0 55]);
