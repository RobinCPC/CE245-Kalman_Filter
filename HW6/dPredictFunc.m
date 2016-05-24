%% HW6
% -------------------------
% CMPE 245 - Spring, 2016
% Chien-Pin Chen
% 05/17/2016
% MATLAB subfunction for HW6
% --------------------------
function dECovX = dPredictFunc(t, Z)
% setting initial value
g = 9.8;
intensity = 1000;
mu0 = 1220;
c = 10263;
u = [0; -g; 0];
b = [0; 0; 1];

r1 = 1000;
r2 = 500;
R = 5;


% compute A(t) (use symbolic to derive first)
%syms x dx beta
%X = [x, dx, beta];
x = Z(1,1);
dx = Z(2,1);
beta = Z(3,1); 

f = [dx; mu0*exp(-x/c)*dx^2/(2*beta); 0];
Ex_t = f + u;
A = [0 ,1 ,0;
     -(610*dx^2)/(10263*beta*exp(x/10263)), (1220*dx)/(beta*exp(x/10263)), -(610*dx^2)/(beta^2*exp(x/10263));
     0, 0,0];

P = zeros(3,3);
P(1,1) = Z(4,1);
P(2,2) = Z(5,1);
P(3,3) = Z(6,1);
P(1,2) = Z(7,1);
P(2,1) = P(1,2);
P(1,3) = Z(8,1);
P(3,1) = P(1,3);
P(2,3) = Z(9,1);
P(3,2) = P(2,3);

dP = A*P+P*A'+b*b'*intensity;

dECovX = [Ex_t(1); Ex_t(2); Ex_t(3); dP(1,1); dP(2,2); dP(3,3); dP(1,2); dP(1,3); dP(2,3)];