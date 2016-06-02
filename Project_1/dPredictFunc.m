%% HW6
% -------------------------
% CMPE 245 - Spring, 2016
% Chien-Pin Chen
% 05/17/2016
% MATLAB subfunction for HW6
% --------------------------
function dECovX = dPredictFunc(t, Z)
% setting initial value
b = [0; 0; 1; 1];
intensity = 12.25^2;


% compute A(t) (use symbolic to derive first)
%syms x dx beta
%X = [x, dx, beta];
x1 = Z(1,1); x2 = Z(2,1); x3 = Z(3,1); x4 = Z(4,1);

f = [x3*cos(x4); x3*sin(x4); 0; 0];
Ex_t = f;
A = [0, 0, cos(x4), -x3*sin(x4);
     0, 0, sin(x4), x3*cos(x4);
     0, 0, 0, 0;
     0, 0, 0, 0];

P = zeros(4,4);
P(1,1) = Z(5,1);
P(2,2) = Z(6,1);
P(3,3) = Z(7,1);
P(4,4) = Z(8,1);

P(1,2) = Z(9,1);
P(2,1) = P(1,2);
P(1,3) = Z(10,1);
P(3,1) = P(1,3);
P(1,4) = Z(11,1);
P(4,1) = P(1,4);

P(2,3) = Z(12,1);
P(3,2) = P(2,3);
P(2,4) = Z(13,1);
P(4,2) = P(2,4);
P(3,4) = Z(14,1);
P(4,3) = P(3,4);


dP = A*P+P*A'+b*b'*intensity;

dECovX = [Ex_t(1); Ex_t(2); Ex_t(3); Ex_t(4); 
          dP(1,1); dP(2,2); dP(3,3); dP(4,4); 
          dP(1,2); dP(1,3); dP(1,4); dP(2,3);
          dP(2,4); dP(3,4)];
 
