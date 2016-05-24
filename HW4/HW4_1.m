%% HW4-1
% -------------------------
% CMPE 245 - Spring, 2016
% Chien-Pin Chen
% 05/08/2016
% MATLAB script for HW4-1
% --------------------------
clear all;

% setting intial value and constant matrix
A = [0.2 0.4; -0.4 1];
b = [1 0; 0 1];
eX_0 = [10; 20];
D = [0;1];


%% compute the expected value of X(k+1)
% vector to storage X(k+1) 
eX_k = zeros(2,15);

eX_k(:,1) = eX_0;%A*eX_0 + D;
for i = 2:15
    eX_k(:,i) = A*eX_k(:,i-1) + D;
end
t = (0:1:14);  % set the range of time step.

% initalize figure and its handles for ploting line later
f1 = figure;        % delcare an empty canvas fb: figure for black line 
h_f1 = gca;         % get handle value (pointer) of above canvs
set(h_f1, 'fontsize', 12, 'FontName','Times New Roman');
set(h_f1,'NextPlot', 'add');
xlabel ('Time step: k');
ylabel ('E\{ x(k) \}');
ylim([2 21]);
title('The evolution of expected value of x(k)');

plot(h_f1, t, eX_k(1,:), 'LineWidth', 2, 'DisplayName','E\{x_1 (k)\}');
plot(h_f1, t, eX_k(2,:), 'r', 'LineWidth', 2,'DisplayName','E\{x_2 (k)\}');
legend show;

%% compute the variances of X(k+1)
% setting intial variance and variable vector
P_0 = eye(2)*40;
Q = eye(2);
P_k = zeros(2, 2, 15);

P_k(:, :, 1) = P_0; %A*P_0*A' + b*Q*b';
for i = 2:length(P_k(1,1,:))
    P_k(:, :, i) = A*P_k(:, :, i-1)*A' + b*Q*b';
end

t = (0:1:14);  % set the range of time step.

% initalize figure and its handles for ploting line later
f3 = figure;        % delcare an empty canvas fb: figure for black line 
h_f3 = gca;         % get handle value (pointer) of above canvs
set(h_f3, 'fontsize', 12, 'FontName','Times New Roman');
set(h_f3,'NextPlot', 'add');
xlabel ('Time step: k');
ylabel ('var( x(k) )');
ylim([0 50]);
title('The evolution of variance of x(k)');

var_x1 = reshape( P_k(1,1,:), 1, []);  % extract the variance of x1 from P_k
var_x2 = reshape( P_k(2,2,:), 1, []);  % extract the variance of x2 from P_k
plot(h_f3, t, var_x1, 'LineWidth', 2, 'DisplayName','var\{x_1 (k)\}');
plot(h_f3, t, var_x2, 'r', 'LineWidth', 2, 'DisplayName','var\{x_2 (k)\}');
legend show;
