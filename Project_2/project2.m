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

c = 1;
g = 9.81;

R = eye(4)*12.25; 
% [12.25 0;
%      0 12.25]; % do we get this ?? or guess by ourself

Q = 12.25/2;

% setting initial x_0, P_0
x0 = xm(1);
y0 = ym(1);
z0 = zm(1);
pit0 = pitch(1);
rol0 = roll(1);
yaw0 = yawm(1);
Ex0 = [ x0 y0 z0 pit0 rol0 yaw0]';

P11 = 12.25;    % var{x}
P22 = P11;      % var{y}
P33 = P11;      % var{z}
P12 = 0;        % temp set to zero
P13 = 0;        % temp set to zero
P23 = 0;        % temp set to zero
Px0 = eye(length(Ex0)) * P11;
% Px0 = [P11 P12 P13;
%        P12 P22 P23;
%        P13 P23 P33];
   

% array to storge E{x_k} and P_k after update 
Ex_k = zeros(length(Ex0), length(xm)+1);
P_k = zeros(length(Ex0), length(Ex0), length(xm)+1);

% Use Extend Kalman Filter (EKF)
Ex_k(:,1) = Ex0;
P_k(:,:,1) = Px0;


for k = 1 : length(zm)
    % Prediction
    % discrete method
    pit = Ex_k(4, k);
    rol = Ex_k(5, k);
    yaw = Ex_k(6, k);
    
    f = [c*g*( sin(yaw)*sin(rol) + cos(yaw)*cos(rol)*sin(pit) );
         c*g*( -1*cos(yaw)*sin(rol) + sin(yaw)*cos(rol)*sin(pit) );
         c*g*cos(rol)*cos(pit) - g;
         0;
         0;
         0];
    
    phi = [0, 0, 0, c*g*cos(pit)*cos(rol)*cos(yaw),  c*g*(cos(rol)*sin(yaw) - cos(yaw)*sin(pit)*sin(rol)), c*g*(cos(yaw)*sin(rol) - cos(rol)*sin(pit)*sin(yaw)) ;
           0, 0, 0, c*g*cos(pit)*cos(rol)*sin(yaw),  c*g*(-cos(rol)*cos(yaw) + sin(pit)*sin(rol)*sin(yaw)), c*g*(sin(rol)*sin(yaw) + cos(rol)*cos(yaw)*sin(pit)) ;
           0, 0, 0,     -c*g*cos(rol)*sin(pit),                     -c*g*cos(pit)*sin(rol),                                          0 ;
           0, 0, 0, 1, 0, 0;
           0, 0, 0, 0, 1, 0;
           0, 0, 0, 0, 0, 1];
    
    Ex_k(:, k+1) = Ex_k(:, k)+f;
    P_k(:,:,k+1) = phi*P_k(:,:,k)*phi'+ eye(length(Ex0))*Q;
    
    % Observation
    H = [ 1, 0, 0, 0, 0, 0;
      0, 1, 0, 0, 0, 0;
      0, 0, 1, 0, 0, 0;
      0, 0, 0, 0, 0, 1];
    x_k1 = [Ex_k(1,k+1) Ex_k(2,k+1) Ex_k(3,k+1) 0 0 Ex_k(6,k+1)]';
    z_k1 = [xm(k) ym(k) zm(k) yawm(k)]';
    
    %update K and P(+)
    K = P_k(:,:,k+1) * H'* inv(H*P_k(:,:,k+1)*H'+R);
    P_k(:,:,k+1) = (eye(length(Ex0)) - K*H)*P_k(:,:,k+1);
    Ex_k(:,k+1) = Ex_k(:,k+1) + K*(z_k1 - H*x_k1);
    
end

%figure;
hold on;
plot3(Ex_k(1,:), Ex_k(2,:), Ex_k(3,:),'r','DisplayName','Ex,Ey,Ez');
title('Ex,Ey,Ez');
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
   
