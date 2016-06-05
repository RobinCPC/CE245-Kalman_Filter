%% Project 2 Part 3
% -------------------------
% CMPE 245 - Spring, 2016
% Chien-Pin Chen
% 05/30/2016
% MATLAB script for Project 2 Part 3
% --------------------------
clear all;

load dataSet2.mat;

plot(yaw, 'DisplayName', 'measure yaw');

% set environment parameter
Ts = 0.1;
g = 9.81;
s_pnt = 1;

R = [1 0 0;
     0 1 0;
     0 0 1];
 
 Q = eye(3);
 
 % setting initial x_0, P_0
x0 = xm(s_pnt);
y0 = ym(s_pnt);
dx0 = 0;%(xm(s_pnt+1)-xm(s_pnt))/Ts;
dy0 = 0;%(ym(s_pnt+1)-ym(s_pnt))/Ts;
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

% Use Extend Kalman Filter (EKF)
Ex_k(:,s_pnt) = Ex0;
P_k(:,:,s_pnt) = Px0;

% start from 45
for k = s_pnt : length(xm)
    % Prediction
    % discrete method
    if k == 61
        continue;
    end
    
    x3 = Ex_k(3, k);     
    x4 = Ex_k(4, k);
    rol = roll(k);
    pit = pitch(k);
    yaw_k = yaw(k);
    
    s = (255/6000) * thrust(k);
    c0 = 10000; % Part2: tuning c0 to get min. residual
    c = (0.000409*s^2+0.1405*s-0.099)/c0;
    c = 1;
    
    f = [x3; x4;
         c*g*( sind(yaw_k)*sind(rol) + cosd(yaw_k)*cosd(rol)*sind(pit) );
         c*g*( sind(yaw_k)*cosd(rol)*sind(pit) -cosd(yaw_k)*sind(rol) );
         0]*Ts;
    
     phi = [1 0 Ts 0 0;
            0 1 0 Ts 0;
            0 0 1 0 0;
            0 0 0 1 0;
            0 0 0 0 1];
        
    ga = [0 0 0;
          0 0 0;
          1*sqrt(Ts) 0 0;
          0 1*sqrt(Ts) 0;
          0 0 1*sqrt(Ts)];
        
    Ex_k(:, k+1) = Ex_k(:, k)+f;
    P_k(:,:,k+1) = phi*P_k(:,:,k)*phi'+ ga*Q*ga'; 

    % Observation
    H = [ 1, 0, 0, 0, 0;
          0, 1, 0, 0, 0;
          0, 0, 0, 0, 1];
    
    x_k1 = [Ex_k(1,k+1) Ex_k(2,k+1) 0 0 Ex_k(5,k+1)]';
    z_k1 = [xm(k) ym(k) yaw(k)]';
  
    %update K and P(+)
    K = P_k(:,:,k+1) * H'/(H*P_k(:,:,k+1)*H'+R);
    P_k(:,:,k+1) = (eye(length(Ex0)) - K*H)*P_k(:,:,k+1);
    Ex_k(:,k+1) = Ex_k(:,k+1) + K*(z_k1 - H*x_k1);    
    
end

hold on;
plot(Ex_k(5,:), 'r','DisplayName', 'predict yaw');
legend show

figure
plot(xm, 'DisplayName','meausre x');
hold on
plot(Ex_k(1,:),'r','DisplayName','predict x');
legend show

figure
plot(ym, 'DisplayName','meausre y');
hold on
plot(Ex_k(2,:),'r','DisplayName','predict y');
legend show