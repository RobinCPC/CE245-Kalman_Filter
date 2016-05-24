%% HW6
% -------------------------
% CMPE 245 - Spring, 2016
% Chien-Pin Chen
% 05/17/2016
% MATLAB script for HW6
% --------------------------
clear all;

% setting initial value
g = 9.8;
intensity = 1000;
mu0 = 1220;
c = 10263;
T = 1;
u = [0; -g; 0];
b = [0; 0; 1];

r1 = 1000;
r2 = 500;
R = 5;
z_k  = [9055, 8560, 7963, 7467, 7000, 6378, 5885, 5400, 4928, 4503];

% syms x dx beta
% %dx = diff(x);
% X = [x, dx, beta];
% f = [dx; mu0*exp(-x/c)*dx^2/(2*beta); 0];
% A = jacobian(f, X);
% h = sqrt(r1^2 + (x - r2)^2);
% H = jacobian(h, X);


Ex0 = 10000;        % E{x(0)}
Edx0 = -500;        % E{dx(0)}
Ebeta0 = 6*10^7;    % E{beta(0)}

P11 = 50;       % var{x_0}
P22 = 200;      % var{dx_0}
P33 = 2*10^12;  % var{beta0}
P12 = 0;        % temp set to zero, but x, dx should correlated
P13 = 0;        % zero, beta and x are independant
P23 = 0;        % zero, beta and dx are independant 

% array to storge E{x_k} and P_k after update 
Ex_k = zeros(3, length(z_k)+1);
P_k = zeros(3, 3, length(z_k)+1);

% Use Extend Kalman Filter (EKF)
Ex_k(:,1) = [Ex0; Edx0; Ebeta0];
P_k(:,:,1) = [P11 P12 P13;
              P12 P22 P23;
              P13 P23 P33];
X_in = [Ex0; Edx0; Ebeta0; P11; P22; P33; P12; P13; P23];

% extra array to collect all result from ode45
tot_T=[0];
tot_X=X_in';
for k = 1: length(z_k)
    % prediction
    [t, solXs] = ode45(@dPredictFunc, [(k-1)*T k*T], X_in);
    solX = solXs(end,:)';
    Ex_k(:,k+1) = [solX(1); solX(2); solX(3)];
    P_k(:,:,k+1) = [solX(4) solX(7) solX(8);
                    solX(7) solX(5) solX(9);
                    solX(8) solX(9) solX(6)];
    % observation
    x_k1 = solX(1);
    H = [ (x_k1 -r2)/( r1^2 + (x_k1 - r2)^2 )^(1/2), 0, 0];
    %H = [(2*x_k1 - 1000)/(2*((x_k1 - 500)^2 + 1000000)^(1/2)) 0 0];
    
    % update K and P(+)
    K = P_k(:,:,k+1) * H'* inv(H*P_k(:,:,k+1)*H'+R);
    P_k(:,:,k+1) = (eye(3) - K*H)*P_k(:,:,k+1);
    %Ex_k(:,k+1) = Ex_k(:,k+1) + K*([z_k(k) 0 0]' - [sqrt(r1^2+(x_k1 - r2)^2) 0 0]');
    Ex_k(:,k+1) = Ex_k(:,k+1) + K*(z_k(k) - sqrt(r1^2+(x_k1 - r2)^2));
    
    X_in = [Ex_k(1,k+1); Ex_k(2,k+1); Ex_k(3,k+1); P_k(1,1,k+1); P_k(2,2,k+1); 
            P_k(3,3,k+1); P_k(1,2,k+1); P_k(1,3,k+1); P_k(2,3,k+1)];
        
   % collecting all ode results
   tot_T = cat(1, tot_T, t(2:end));
   tot_X = cat(1, tot_X, solXs(2:end, :));
end


% compute x_m with z_k use equation (6)
x_m = zeros(1, length(z_k));
for i = 1 : length(z_k)
    x_m(i) = sqrt(z_k(i)^2 - r1^2) + r2;
end

% ploting measurement x_m and the EKF estimated x(t)
t = (0:1:10);

% initalize figure and its handles for ploting line later
f1 = figure;        % delcare an empty canvas
h_f1 = gca;         % get handle value (pointer) of above canvs
set(h_f1, 'fontsize', 12, 'FontName','Times New Roman');
set(h_f1,'NextPlot', 'add');
title('The measurement x_m v.s. the estimate x(t)');
xlabel('Time (sec)');
ylabel ('Altitude (m)');

plot(h_f1, t, [Ex0 x_m], 'DisplayName','Measurement x_m (k)');
plot(h_f1, t, Ex_k(1,:), 'r', 'DisplayName','Estimation of x(t)');
plot(h_f1, tot_T, tot_X(:,1), 'g', 'DisplayName','ODE results of x(t)');
legend('show');

% initalize figure and its handles for ploting line later
f2 = figure;        % delcare an empty canvas
h_f2 = gca;         % get handle value (pointer) of above canvs
set(h_f2, 'fontsize', 12, 'FontName','Times New Roman');
set(h_f2,'NextPlot', 'add');
title('The estimated velocity dx(t)');
xlabel('Time (sec)');
ylabel ('Velocity (m/s)');
plot(h_f2, t, Ex_k(2,:), 'DisplayName','Estimation of  dx(t)');
plot(h_f2, tot_T, tot_X(:,2), 'g', 'DisplayName','ODE results of dx(t)');
legend('show');

% initalize figure and its handles for ploting line later
f3 = figure;        % delcare an empty canvas
h_f3 = gca;         % get handle value (pointer) of above canvs
set(h_f3, 'fontsize', 12, 'FontName','Times New Roman');
set(h_f3,'NextPlot', 'add');
title('The estimated ballistic coefficient \beta (t)');
xlabel('Time (sec)');
ylabel ('Ballistic Coefficient (g/ms^2)');
plot(h_f3, t, Ex_k(3,:), 'DisplayName','Estimation of \beta (t)');
plot(h_f3, tot_T, tot_X(:,3), 'g', 'DisplayName','ODE results of \beta (t)');
legend('show');
