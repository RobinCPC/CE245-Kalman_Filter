%% HW7
% -------------------------
% CMPE 245 - Spring, 2016
% Chien-Pin Chen
% 05/27/2016
% MATLAB script for HW7
% --------------------------
clear all;
% measurement data and plot them on one figure
load robotMes.mat;

% plotting the measurement 
plot(y(1,:), 'LineWidth', 2, 'DisplayName','y_1 (k) The position of the robot');
hold on
plot(y(2,:), 'r', 'LineWidth', 2, 'DisplayName','y_2 (k) The position of the wall');
set(gca, 'fontsize', 12, 'FontName','Times New Roman');
title('The measurement with noise');
xlabel ('Time step: k (second)');
ylabel ('The measurement positions from the corner (cm)');
legend('show', 'Location', 'SouthEast');

% setting initial x_0, P_0 
%x_0 = [40; 0; 200]; % my first guess. (integer according y1(0) and y2(0))
x_0 = [y(1,1); 0.2; y(1,1)+y(2,1)]; % more reasoning guess (after course)
% P_0 = [10 0 0;
%         0 0 0;    % my first guess. (set as variance of the measurement 
%         0 0 10];  % for P11 P33. and set others as zeros)
P_0 = [10 0 10;
        0 1 0;
        10 0 10*2];
    
% setting some constant matrix in kalman filter (phi, gamma, H, R, Q)
phi = [1 1 0;
       0 1 0;
       0 0 1];
ga = [0 0.1 0]';
%u = [ 0 0.8 0]';
Q = 1;

H = [1 0 0;
    -1 0 1];
R = [10 0;
     0 10];

% array to compute E{x_k} and P_k 
Ex_k = zeros(3, length(y)+1);
P_k = zeros(3, 3, length(y)+1);
 
% array to storge E{x_k}(-) and P_k(-) before update (prior)
Ex_kn = zeros(3, length(y)+1);
P_kn = zeros(3, 3, length(y)+1);
 
% array to storge E{x_k}(+) and P_k(+) after update (posterior) 
Ex_kp = zeros(3, length(y)+1);
P_kp = zeros(3, 3, length(y)+1);

% array to storage E{x_k} and P_k before and after update
Ex_arr = zeros(3, 2*length(y)+1);
P_arr = zeros(3, 3, 2*length(y)+1);
K_arr = zeros(3,2,length(y));   % storage KF Gain (for trends)

% Do forward kalman filter to compute estimation and the varinace
Ex_k(:,1) = x_0;
P_k(:,:,1) = P_0;
Ex_kn(:,1) = x_0;
P_kn(:,:,1) = P_0;
Ex_kp(:,1) = x_0;
P_kp(:,:,1) = P_0;
Ex_arr(:,1) = x_0;
P_arr(:,:,1) = P_0;
for k=1: length(y)
    % prediction (prior)
    Ex_k(:,k+1) = phi* Ex_k(:,k);
    P_k(:,:,k+1) = phi*P_k(:,:,k)*phi' + ga*Q*ga';
    
    Ex_arr(:,k*2) = Ex_k(:,k+1);
    P_arr(:,:,k*2) = P_k(:,:,k+1);
    Ex_kn(:,k+1) = Ex_k(:,k+1);
    P_kn(:,:,k+1) = P_k(:,:,k+1);
    
    % observation

    % update and compute K and P(+) (posterior)
    K = P_k(:,:,k+1) * H'* inv(H*P_k(:,:,k+1)*H'+R);
    P_k(:,:,k+1) = (eye(3) - K*H)*P_k(:,:,k+1);
    
    Ex_k(:,k+1) = Ex_k(:,k+1) + K*(y(:,k) - H*Ex_k(:,k+1));
    
    K_arr(:,:,k) = K;
    Ex_arr(:,k*2+1) = Ex_k(:,k+1);
    P_arr(:,:,k*2+1) = P_k(:,:,k+1);
    Ex_kp(:,k+1) = Ex_k(:,k+1);
    P_kp(:,:,k+1) = P_k(:,:,k+1);
end


% Do kalman smoother
% array to compute E{x_k} and P_k 
Ex_ks = zeros(3, length(y)+1);
P_ks = zeros(3, 3, length(y)+1);

Ex_ks(:,1) = x_0;
P_ks(:,:,1) = P_0;

Ex_ks(:, end) = Ex_kp(:, end);
P_ks(:,:, end) = P_kp(:, :, end);

for k = length(y):-1:1
    A_k = P_kp(:, :, k)*phi'*inv(P_kn(:, :, k+1));
    
    Ex_ks(:,k) = Ex_kp(:,k) + A_k*( Ex_ks(:,k+1) - Ex_kn(:,k+1) );
    P_ks(:,:,k) = P_kp(:,:,k) + A_k*( P_ks(:,:,k+1) - P_kn(:,:,k+1) )*A_k';
end


% plot the result of the forward Kalman fillter and RTS backward 
% iterations for the robot position
t = (0:1:100);
f1 = figure;        % delcare an empty canvas
h_f1 = gca;         % get handle value (pointer) of above canvs
set(h_f1, 'fontsize', 12, 'FontName','Times');
set(h_f1,'NextPlot', 'add');
title('The estimated robot position x(t)');
xlabel('Time step: k (second)');
ylabel ('The esitmated position from the corner (cm)');

plot(h_f1, t, Ex_kp(1,:), 'DisplayName','Forward: E\{x1\}');
plot(h_f1, t, Ex_ks(1,:), 'r', 'DisplayName','RTS: E\{x1\}');
legend('show', 'Location', 'SouthEast');

% plot the result of the forward Kalman fillter and RTS backward 
% iterations for the robot velocity
f2 = figure;        % delcare an empty canvas
h_f2 = gca;         % get handle value (pointer) of above canvs
set(h_f2, 'fontsize', 12, 'FontName','Times');
set(h_f2,'NextPlot', 'add');
title('The estimated robot position dx(t)');
xlabel('Time step: k (second)');
ylabel ('The estimated velocity of the robot (cm/s)');

plot(h_f2, t, Ex_kp(2,:), 'DisplayName','Forward: E\{x2\}');
plot(h_f2, t, Ex_ks(2,:), 'r', 'DisplayName','RTS: E\{x2\}');
legend('show', 'Location', 'SouthEast');

% plot the result of the forward Kalman fillter and RTS backward 
% iterations for the wall position
f3 = figure;        % delcare an empty canvas
h_f3 = gca;         % get handle value (pointer) of above canvs
set(h_f3, 'fontsize', 12, 'FontName','Times');
set(h_f3,'NextPlot', 'add');
title('The estimated distance from the corner to the wall');
xlabel('Time step: k (second)');
ylabel ('The estimated position of the wall (cm)');
plot(h_f3, t, Ex_kp(3,:), 'DisplayName','Forward: E\{x3\}');
plot(h_f3, t, Ex_ks(3,:), 'r', 'DisplayName','RTS: E\{x3\}');
legend('show', 'Location', 'SouthEast');

% plot the variance of the robot position from Kalman filter and Kalman smoother.
var_x1p = reshape(P_kp(1,1,:), 1, []);
var_x1s = reshape(P_ks(1,1,:), 1, []);
f4 = figure;        % delcare an empty canvas
h_f4 = gca;         % get handle value (pointer) of above canvs
set(h_f4, 'fontsize', 12, 'FontName','Times');
set(h_f4,'NextPlot', 'add');
title('The vairance of the robot position from Kalman filter and smoother');
xlabel('Time step: k (second)');
ylabel ('The vairance of the robot position (cm)');
xlim([-1, 100]);
plot(h_f4, t, var_x1p, 'DisplayName','Forward: var\{x1\}');
plot(h_f4, t, var_x1s, 'r', 'DisplayName','RTS: var\{x1\}');
legend('show', 'Location', 'NorthEast');
