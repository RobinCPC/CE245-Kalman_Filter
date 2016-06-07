%% Project 2  Part1 and 2
% -------------------------
% CMPE 245 - Spring, 2016
% Chien-Pin Chen
% 05/30/2016
% MATLAB script for Project 2
% --------------------------
clear all;

load data3D.mat;



%T = 0.1;
%c = 1e-16;
g = 9.81;

R = [.1^2 0 0 0;
     0 .1^2 0 0;
     0 0 .1^2 0;
     0 0 0 1];
%eye(4)*1.225*1; 
% [12.25 0;
%      0 12.25]; % do we get this ?? or guess by ourself

% ga = [0 0 0 0;
%       0 0 0 0;
%       0 0 0 0;
%       sqrt(T) 0 0 0;
%       0 sqrt(T) 0 0;
%       0 0 sqrt(T) 0;
%       0 0 0 sqrt(T)];
Q = eye(4);

% setting initial x_0, P_0
x0 = xm(45);
y0 = ym(45);
z0 = zm(45);
dx0 = (xm(45)-xm(44))/dtrec(45);
dy0 = (ym(45)-ym(44))/dtrec(45);
dz0 = (zm(45)-zm(44))/dtrec(45);
yaw0 = yawm(45);
Ex0 = [ x0 y0 z0 dx0 dy0 dz0 yaw0]';

P11 = 0.1^2%1.225*.1;    % var{x}
P22 = P11;      % var{y}
P33 = P11;      % var{z}
P44 = 2*P11/(0.0085^2);
P55 = 2*P11/(0.0085^2);
P66 = 2*P11/(0.0085^2);
P14 = P11;
P25 = P11;
P36 = P11;
P12 = 0;        % temp set to zero
P13 = 0;        % temp set to zero
P23 = 0;        % temp set to zero
P45 = 0;        % temp set to zero
P46 = 0;        % temp set to zero
P56 = 0;        % temp set to zero
Px0 = eye(length(Ex0)) * P11;
% Px0 = [P11 P12 P13;
%        P12 P22 P23;
%        P13 P23 P33];
Px0(4,4) = 2*P11/(0.0085^2);
Px0(5,5) = 2*P11/(0.0085^2);
Px0(6,6) = 2*P11/(0.0085^2);
Px0(7,7) = 10*(pi/180)^2;%0.5^2;
Px0(1,4) = P11;
Px0(4,1) = Px0(1,4);
Px0(2,5) = P11;
Px0(5,2) = Px0(2,5);
Px0(3,6) = P11;
Px0(6,3) = Px0(3,6);


% array to storge E{x_k} and P_k after update 
Ex_k = zeros(length(Ex0), length(xm)+1);
P_k = zeros(length(Ex0), length(Ex0), length(xm)+1);
S = 0;  % compute c0 part2
i_sum = 0;

% Use Extend Kalman Filter (EKF)
Ex_k(:,45) = Ex0;
P_k(:,:,45) = Px0;

% start from 45
for k = 45 : length(zm)
    % Prediction
    % discrete method
    
%     pit = Ex_k(4, k);
%     rol = Ex_k(5, k);
%     yaw = Ex_k(6, k);

    x4 = Ex_k(4, k);     
    x5 = Ex_k(5, k);
    x6 = Ex_k(6, k);
    rol = roll(k);
    pit = pitch(k);
    x7 = yawm(k);
    dt = dtrec(k);
    
    s = (255/6000) * thrust(k);
    c0 = 5*1e+4; % Part2: tuning c0 to get min. residual
    c = (0.000409*s^2+0.1405*s-0.099)/c0;
    c = 1;
    f = [x4; x5; x6;
         c*g*( sind(x7)*sind(rol) + cosd(x7)*cosd(rol)*sind(pit) );
         c*g*( sind(x7)*cosd(rol)*sind(pit) -cosd(x7)*sind(rol) );
         c*(g*cosd(rol)*cosd(pit) - g);
         0]*dt; 
% [c*g*( sin(yaw)*sin(rol) + cos(yaw)*cos(rol)*sin(pit) );
%  c*g*( -1*cos(yaw)*sin(rol) + sin(yaw)*cos(rol)*sin(pit) );
%  c*g*cos(rol)*cos(pit) - g;
%  0;
%  0;
%  0];
    
    phi = [1 0 0 dt 0 0 0;
           0 1 0 0 dt 0 0;
           0 0 1 0 0 dt 0;
           0 0 0 1 0 0 dt*c*g*( cosd(x7)*sind(rol)-sind(x7)*cosd(rol)*sind(pit) );
           0 0 0 0 1 0 dt*c*g*( cosd(x7)*cosd(rol)*sind(pit)+sind(x7)*sind(pit) );
           0 0 0 0 0 1 0;
           0 0 0 0 0 0 1];
    
% [0, 0, 0, c*g*cos(pit)*cos(rol)*cos(yaw),  c*g*(cos(rol)*sin(yaw) - cos(yaw)*sin(pit)*sin(rol)), c*g*(cos(yaw)*sin(rol) - cos(rol)*sin(pit)*sin(yaw)) ;
% 0, 0, 0, c*g*cos(pit)*cos(rol)*sin(yaw),  c*g*(-cos(rol)*cos(yaw) + sin(pit)*sin(rol)*sin(yaw)), c*g*(sin(rol)*sin(yaw) + cos(rol)*cos(yaw)*sin(pit)) ;
% 0, 0, 0,     -c*g*cos(rol)*sin(pit),                     -c*g*cos(pit)*sin(rol),                                          0 ;
% 0, 0, 0, 1, 0, 0;
% 0, 0, 0, 0, 1, 0;
% 0, 0, 0, 0, 0, 1];
    ga = [0 0 0 0;
          0 0 0 0;
          0 0 0 0;
          .1*sqrt(dt) 0 0 0;
          0 .5*sqrt(dt) 0 0;
          0 0 .1*sqrt(dt) 0;
          0 0 0 0.01*sqrt(dt)];
    
    Ex_k(:, k+1) = Ex_k(:, k)+f;
    P_k(:,:,k+1) = phi*P_k(:,:,k)*phi'+ ga*Q*ga';
    
    % Observation
    H = [ 1, 0, 0, 0, 0, 0, 0;
          0, 1, 0, 0, 0, 0, 0;
          0, 0, 1, 0, 0, 0, 0;
          0, 0, 0, 0, 0, 0, 1];
    x_k1 = [Ex_k(1,k+1) Ex_k(2,k+1) Ex_k(3,k+1) 0 0 0 Ex_k(7,k+1)]';
    z_k1 = [xm(k) ym(k) zm(k) yawm(k)]';
    
    % part2 compute c0
    S = S + sum( (z_k1(1:3)-x_k1(1:3)).^2 );
    i_sum = i_sum+1;
    
    %update K and P(+)
    K = P_k(:,:,k+1) * H'/(H*P_k(:,:,k+1)*H'+R);
    P_k(:,:,k+1) = (eye(length(Ex0)) - K*H)*P_k(:,:,k+1);
    Ex_k(:,k+1) = Ex_k(:,k+1) + K*(z_k1 - H*x_k1);
    
end

% get average of S (part2)
S_ave = S/i_sum;

% plot results of EKF
figure;
%subplot(2,1,1);
plot3(Ex_k(1,:), Ex_k(2,:), Ex_k(3,:),'--','DisplayName','Predicted Position');
hold on;
plot3(xm, ym, zm, 'r', 'DisplayName','Measurement Position');
axis equal;
xlabel('X (M)','FontSize',12,'FontName','Times');
ylabel('Y (M)','FontSize',12,'FontName','Times');
zlabel('Z (M)','FontSize',12,'FontName','Times');
h_f1 = gca;         % get handle value (pointer) of above canvs
set(h_f1, 'fontsize', 12, 'FontName','Times');
%legend('show');
title('Robot Position','FontSize',16)



figure;
subplot(4,1,1);
plot(Ex_k(1,:),'--','DisplayName','Predicted X');
hold on
plot(xm, 'r', 'DisplayName','Measurement X');
xlabel('Time (sec)','FontSize',12,'FontName','Times');
ylabel('X (M)','FontSize',12,'FontName','Times');
h_f2 = gca;         % get handle value (pointer) of above canvs
set(h_f2, 'fontsize', 12, 'FontName','Times');
%legend('show')
title('X position','FontSize',16)

subplot(4,1,2);
plot(Ex_k(2,:),'--','DisplayName','Predicted Y');
hold on
plot(ym, 'r', 'DisplayName','Measurement Y');
xlabel('Time (sec)','FontSize',12,'FontName','Times');
ylabel('Y (M)','FontSize',12,'FontName','Times');
h_f3 = gca;         % get handle value (pointer) of above canvs
set(h_f3, 'fontsize', 12, 'FontName','Times');
%legend('show')
title('Y position','FontSize',16)

subplot(4,1,3);
plot(Ex_k(3,:),'--','DisplayName','Predicted Z');
hold on
plot(zm, 'r', 'DisplayName','Measurement Z');
xlabel('Time (sec)','FontSize',12,'FontName','Times');
ylabel('Z (M)','FontSize',12,'FontName','Times');
h_f4 = gca;         % get handle value (pointer) of above canvs
set(h_f4, 'fontsize', 12, 'FontName','Times');
%legend('show')
title('Z position','FontSize',16)

subplot(4,1,4);
plot(Ex_k(7,:),'--','DisplayName','Predicted Yaw');
hold on
plot(yawm, 'r', 'DisplayName','Measurement Yaw');
xlabel('Time (sec)','FontSize',12,'FontName','Times');
ylabel('Angle (degree)','FontSize',12,'FontName','Times');
h_f5 = gca;         % get handle value (pointer) of above canvs
set(h_f5, 'fontsize', 12, 'FontName','Times');
%legend('show')
title('Yaw angle (\psi)','FontSize',16)


% plot variance of expected yaw
var_x = reshape(P_k(1,1,:), 1, []);
var_y = reshape(P_k(2,2,:), 1, []);
var_z = reshape(P_k(3,3,:), 1, []);
var_yaw = reshape(P_k(7,7,:), 1, []);
figure
plot(var_yaw, 'k-.','DisplayName','VAR \psi');
hold on
plot(var_x, 'DisplayName','VAR x');
plot(var_y, 'r', 'DisplayName','VAR y');
plot(var_z, 'g','DisplayName','VAR z');


   
