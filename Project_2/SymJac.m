%% Project2
% -------------------------
% CMPE 245 - Spring, 2016
% Chien-Pin Chen
% 05/30/2016
% MATLAB script for Project2
% --------------------------
% script for do jacobian with symbolic
clear all;

syms dx dy dz rol pit yaw
syms c g

X = [dx dy dz pit rol yaw];
f = [c*g*( sin(yaw)*sin(rol) + cos(yaw)*cos(rol)*sin(pit) );
     c*g*( -1*cos(yaw)*sin(rol) + sin(yaw)*cos(rol)*sin(pit) );
     c*g*cos(rol)*cos(pit) - g];
A = jacobian(f,X);

h = [dx;
     dy;
     dz;
     yaw];
 H = jacobian(h,X);