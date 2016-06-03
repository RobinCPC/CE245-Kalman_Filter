%% Project 2
% -------------------------
% CMPE 245 - Spring, 2016
% Chien-Pin Chen
% 05/30/2016
% MATLAB script for Project 2
% --------------------------
clear all;

load data3D.mat;

plot3(xm, ym, zm,'DisplayName','xm,ym,zm');
title('xm,ym,zm');

%T = 0.1;
%c = 1e-16;
g = 9.81;

R = [1 0 0 0;
     0 1 0 0;
     0 0 0.1 0;
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

P11 = 1.225*.01;    % var{x}
P22 = P11;      % var{y}
P33 = P11;      % var{z}
P12 = 0;        % temp set to zero
P13 = 0;        % temp set to zero
P23 = 0;        % temp set to zero
Px0 = eye(length(Ex0)) * P11;
% Px0 = [P11 P12 P13;
%        P12 P22 P23;
%        P13 P23 P33];

Px0(1,4) = 2*P11/(0.0085^2);
Px0(4,1) = Px0(1,4);
Px0(2,5) = 2*P11/(0.0085^2);
Px0(5,2) = Px0(2,5);
Px0(3,6) = 2*P11/(0.0085^2);
Px0(6,3) = Px0(3,6);

% array to storge E{x_k} and P_k after update 
Ex_k = zeros(length(Ex0), length(xm)+1);
P_k = zeros(length(Ex0), length(Ex0), length(xm)+1);

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
    yaw = yawm(k);
    dt = dtrec(k);
    
    s = (255/6000) * thrust(k);
    c0 = 1;%100000000*5;
    c = (0.000409*s^2+0.1405*s-0.099)/c0;
    c=1;
    f = [x4; x5; x6;
         c*g*( sind(yaw)*sind(rol) + cosd(yaw)*cosd(rol)*sind(pit) );
         c*g*( sind(yaw)*cosd(rol)*sind(pit) -cosd(yaw)*sind(rol) );
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
           0 0 0 1 0 0 0;
           0 0 0 0 1 0 0;
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
          10*sqrt(dt) 0 0 0;
          0 10*sqrt(dt) 0 0;
          0 0 1*sqrt(dt) 0;
          0 0 0 1*sqrt(dt)];
    
    Ex_k(:, k+1) = Ex_k(:, k)+f;
    P_k(:,:,k+1) = phi*P_k(:,:,k)*phi'+ ga*Q*ga';
    
    % Observation
    H = [ 1, 0, 0, 0, 0, 0, 0;
          0, 1, 0, 0, 0, 0, 0;
          0, 0, 1, 0, 0, 0, 0;
          0, 0, 0, 0, 0, 0, 1];
    x_k1 = [Ex_k(1,k+1) Ex_k(2,k+1) Ex_k(3,k+1) 0 0 0 Ex_k(7,k+1)]';
    z_k1 = [xm(k) ym(k) zm(k) yawm(k)]';
    
    %update K and P(+)
    K = P_k(:,:,k+1) * H'/(H*P_k(:,:,k+1)*H'+R);
    P_k(:,:,k+1) = (eye(length(Ex0)) - K*H)*P_k(:,:,k+1);
    Ex_k(:,k+1) = Ex_k(:,k+1) + K*(z_k1 - H*x_k1);
    
end

%figure;
hold on;
plot3(Ex_k(1,:), Ex_k(2,:), Ex_k(3,:),'r','DisplayName','Ex,Ey,Ez');
title('Ex,Ey,Ez');

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


figure
plot(zm, 'DisplayName','meausre z');
hold on
plot(Ex_k(3,:),'r','DisplayName','predict z');
legend show

figure
plot(yawm, 'DisplayName','meausre yaw');
hold on
plot(Ex_k(7,:),'r','DisplayName','predict yaw');
legend show

% f = [c*g*( sin(yaw)*sin(rol) + cos(yaw)*cos(rol)*sin(pit) );
%      c*g*( -1*cos(yaw)*sin(rol) + sin(yaw)*cos(rol)*sin(pit) );
%      c*g*cos(rol)*cos(pit) - g;
%      0;
%      0;
%      0];
% 
% phi = [0, 0, 0, c*g*cos(pit)*cos(rol)*cos(yaw),  c*g*(cos(rol)*sin(yaw) - cos(yaw)*sin(pit)*sin(rol)), c*g*(cos(yaw)*sin(rol) - cos(rol)*sin(pit)*sin(yaw)) ;
%        0, 0, 0, c*g*cos(pit)*cos(rol)*sin(yaw),  c*g*(-cos(rol)*cos(yaw) + sin(pit)*sin(rol)*sin(yaw)), c*g*(sin(rol)*sin(yaw) + cos(rol)*cos(yaw)*sin(pit)) ;
%        0, 0, 0,     -c*g*cos(rol)*sin(pit),                     -c*g*cos(pit)*sin(rol),                                          0 ;
%        0, 0, 0, 1, 0, 0;
%        0, 0, 0, 0, 1, 0;
%        0, 0, 0, 0, 0, 1];
%  
%  
% H = [ 1, 0, 0, 0, 0, 0;
%       0, 1, 0, 0, 0, 0;
%       0, 0, 1, 0, 0, 0;
%       0, 0, 0, 0, 0, 1];
   
