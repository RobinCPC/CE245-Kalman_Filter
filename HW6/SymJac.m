%% HW6
% -------------------------
% CMPE 245 - Spring, 2016
% Chien-Pin Chen
% 05/17/2016
% MATLAB script for HW6
% --------------------------
% script for do jacobian with symbolic
clear all;

syms x dx beta
syms mu0 c r1 r2
%dx = diff(x);
X = [x, dx, beta];
f = [dx; mu0*exp(-x/c)*dx^2/(2*beta); 0];
A = jacobian(f, X);
h = sqrt(r1^2 + (x - r2)^2);
H = jacobian(h, X);
