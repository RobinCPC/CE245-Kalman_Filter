%% Project 2 Part 3
% -------------------------
% CMPE 245 - Spring, 2016
% Chien-Pin Chen
% 05/30/2016
% MATLAB script for Project 2 Part 3
% --------------------------
clear all;

load dataSet2.mat;

% set environment parameter
Ts = 0.1;
g = 9.81;
s_pnt = 1;

R = [.1 0 0;
     0 .1 0;
     0 0 1];
 
 Q = eye(3);
 
 % setting initial x_0, P_0
x0 = xm(s_pnt);
y0 = ym(s_pnt);
dx0 = (xm(s_pnt+1)-xm(s_pnt))/Ts;
dy0 = (ym(s_pnt+1)-ym(s_pnt))/Ts;
yaw0 = yaw(s_pnt);
Ex0 = [x0 y0 dx0 dy0 yaw0]';

P11 = 1.225*.01;    % var{x}
Px0 = eye(length(Ex0)) * P11;

Px0(1,3) = 2*P11/(Ts^2);
Px0(3,1) = Px0(1,3);
Px0(2,4) = 2*P11/(Ts^2);
Px0(4,2) = Px0(2,4);

% array to storge E{x_k} and P_k after update 
Ex_k = zeros(length(Ex0), length(xm)+1);
P_k = zeros(length(Ex0), length(Ex0), length(xm)+1);

% catch corrupt data and compute residual
S = 0;  % compute the average of residual
i_sum = 0;
corrupt_yaw = [];    % 1st parameter is time point, 2nd is the value of corrupt measurement yaw 


% Use Extend Kalman Filter (EKF)
Ex_k(:,s_pnt) = Ex0;
P_k(:,:,s_pnt) = Px0;

for k = s_pnt : length(xm)
    % Prediction
    % discrete method
    
    x3 = Ex_k(3, k);     
    x4 = Ex_k(4, k);
    rol = roll(k);
    pit = pitch(k);
    x5 = yaw(k);
    
    s = (255/6000) * thrust(k);
    c0 = 10000; % Part2: tuning c0 to get min. residual
    c = (0.000409*s^2+0.1405*s-0.099)/c0;
    c = 1;
    
    f = [x3; x4;
         c*g*( sind(x5)*sind(rol) + cosd(x5)*cosd(rol)*sind(pit) );
         c*g*( sind(x5)*cosd(rol)*sind(pit) -cosd(x5)*sind(rol) );
         0]*Ts;
    
     phi = [1 0 Ts 0 0;
            0 1 0 Ts 0;
            0 0 1 0 Ts*c*g*( cosd(x5)*sind(rol)-sind(x5)*cosd(rol)*sind(pit) );
            0 0 0 1 Ts*c*g*( cosd(x5)*cosd(rol)*sind(pit)+sind(x5)*sind(pit) );
            0 0 0 0 1];
        
    ga = [0 0 0;
          0 0 0;
          1*sqrt(Ts) 0 0;
          0 1*sqrt(Ts) 0;
          0 0 .5*sqrt(Ts)];
        
    Ex_k(:, k+1) = Ex_k(:, k)+f;
    P_k(:,:,k+1) = phi*P_k(:,:,k)*phi'+ ga*Q*ga'; 

    % Observation
    H = [ 1, 0, 0, 0, 0;
          0, 1, 0, 0, 0;
          0, 0, 0, 0, 1];
    
    x_k1 = [Ex_k(1,k+1) Ex_k(2,k+1) 0 0 Ex_k(5,k+1)]';
    z_k1 = [xm(k) ym(k) yaw(k)]';
    
    if sum( (z_k1(1:3)-x_k1([1:2,5])).^2 ) > 498.16
        % ignore the corrupted measurement of yaw
        x_k1 = [Ex_k(1,k+1) Ex_k(2,k+1) 0 0 0]';
        z_k1 = [xm(k) ym(k) 0]';
        corrupt_yaw = [corrupt_yaw; k, yaw(k)];
    end
    
    % compute the total of residual
    S = S + sum( (z_k1(1:3)-x_k1([1:2,5])).^2 );
    i_sum = i_sum+1;
  
    %update K and P(+)
    K = P_k(:,:,k+1) * H'/(H*P_k(:,:,k+1)*H'+R);
    P_k(:,:,k+1) = (eye(length(Ex0)) - K*H)*P_k(:,:,k+1);
    Ex_k(:,k+1) = Ex_k(:,k+1) + K*(z_k1 - H*x_k1);    
    
end

% get average of residual
S_ave = S/i_sum;

figure;
plot(Ex_k(5,:), '--','DisplayName', 'Predicted Yaw (\psi)');
hold on;
plot(yaw, 'r','DisplayName', 'Measurement Yaw (\psi)');
if isempty(corrupt_yaw) == 0
    plot(corrupt_yaw(:,1), corrupt_yaw(:,2), 'bo', 'MarkerFaceColor', 'green', 'DisplayName', 'Corrupted Yaw (\psi)');
end
xlabel('Time (sec)','FontSize',12,'FontName','Times');
ylabel('Angle (degree)','FontSize',12,'FontName','Times');
h_f1 = gca;         % get handle value (pointer) of above canvs
set(h_f1, 'fontsize', 12, 'FontName','Times');
legend('show', 'Location', 'NorthWest');
title('Yaw angle (\psi)','FontSize',16)

figure
subplot(2,1,1);
plot(Ex_k(1,:),'--','DisplayName','Predicted X');
hold on
plot(xm, 'r', 'DisplayName','Measurement X');
xlabel('Time (sec)','FontSize',12,'FontName','Times');
ylabel('X (M)','FontSize',12,'FontName','Times');
h_f2 = gca;         % get handle value (pointer) of above canvs
set(h_f2, 'fontsize', 12, 'FontName','Times');
legend show;
title('X position','FontSize',16);

subplot(2,1,2);
plot(Ex_k(2,:),'--','DisplayName','Predicted Y');
hold on
plot(ym, 'r','DisplayName','Measurement Y');
xlabel('Time (sec)','FontSize',12,'FontName','Times');
ylabel('Y (M)','FontSize',12,'FontName','Times');
h_f3 = gca;         % get handle value (pointer) of above canvs
set(h_f3, 'fontsize', 12, 'FontName','Times');
legend show;
title('Y position','FontSize',16);


figure
plot(xm, ym, 'DisplayName', 'measure x y');
hold on
plot(Ex_k(1,:), Ex_k(2,:), 'r','DisplayName','predict x y')
legend show