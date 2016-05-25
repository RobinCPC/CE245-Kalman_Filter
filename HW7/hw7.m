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
x_0 = [y(1,1); 1.6; y(1,1)+y(2,1)]; % more reasoning guess (after course)
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

% plot expexted X1 and X3 (expect and update)
t = (0:0.5:100);
figure;
plot(t, Ex_arr(1,:), 'LineWidth', 1, 'DisplayName','E\{ Position of the robot (x_1)\} ');
hold on
plot(t, Ex_arr(3,:), 'r', 'LineWidth', 1, 'DisplayName','E\{ Position of the wall (x_3)\}');
set(gca, 'fontsize', 12, 'FontName','Times New Roman');
title('The estimation of positions of the robot and the wall');
xlabel ('Time step: k (second)');
ylabel ('The esitmated positions from the corner (cm)');
legend('show', 'Location', 'SouthEast');


% plot variance X1 and X3 (expect and update)
var_x1a = reshape(P_arr(1,1,:), 1, []);
var_x3a = reshape(P_arr(3,3,:), 1, []);
figure;
plot(t, var_x1a, 'LineWidth', 1, 'DisplayName','var\{ position of the robot (x_1)\}');
hold on
plot(t, var_x3a, 'r', 'LineWidth', 1, 'DisplayName','var\{ position of the wall (x_3)\}');
set(gca, 'fontsize', 12, 'FontName','Times New Roman');
title('The vairance of positions of the robot and the wall');
xlabel ('Time step: k (second)');
ylabel ('The variance of the positions (cm)');
legend('show');


% plot expexted X1 and X3 (only update)
figure;
plot(Ex_kp(1,:), 'DisplayName','Ex1');
hold on
plot(Ex_kp(3,:), 'r', 'DisplayName','Ex3');
legend('show', 'Location', 'SouthEast');

% plot variance X1 and X3
var_x1 = reshape(P_kp(1,1,:), 1, []);
var_x3 = reshape(P_kp(3,3,:), 1, []);
figure;
plot(var_x1, 'DisplayName','VARx1');
hold on
plot(var_x3, 'r', 'DisplayName','VARx3');
legend('show');
