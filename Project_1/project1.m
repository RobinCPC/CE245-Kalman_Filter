%% Project 1
% -------------------------
% CMPE 245 - Spring, 2016
% Chien-Pin Chen
% 05/26/2016
% MATLAB script for Project 1
% --------------------------
clear all;

% setting prediction function


% read the measurement data
z_k = load('XYData_cm.csv'); %xlsread('XYData_cm.csv'); % another functio to use: dataset

plot( z_k(:,1), z_k(:,2));

%%  EKF
T = 1/3;   % 30 frames per sec, data are acquired per 10 frames

b = [0 0 1 1]';

R = [12.25 0;
     0 12.25]; % do we get this ?? or guess by ourself

x0 = z_k(1,1);
y0 = z_k(1,2);
v0 = 6;     % (z_k(2,1)-z_k(1,1)) / 0.3 
theta0 = atan2( (z_k(2,2)-z_k(1,2)), (z_k(2,1)-z_k(1,1))  );

Ex0 = [ x0 y0 v0 theta0 ]';   %

P11 = 12.25*1;    % var{x}
P22 = P11;    % var{y}
P33 = P11;    % var{v}
P44 = (180/pi)*2*P11;    % var{theta}
P12 = 0;
P13 = 0;
P14 = 0;
P23 = 0;
P24 = 0;
P34 = 0;

Px0 =[P11 P12 P13 P14;
      P12 P22 P23 P24;
      P13 P23 P33 P34;
      P14 P24 P34 P44];

% array to storge E{x_k} and P_k after update 
Ex_k = zeros(4, length(z_k)+1);
P_k = zeros(4, 4, length(z_k)+1);

% Use Extend Kalman Filter (EKF)
Ex_k(:,1) = Ex0; %[Ex0; Edx0; Ebeta0];
P_k(:,:,1) = Px0; %[P11 P12 P13;
                   % P12 P22 P23;
                   % P13 P23 P33];

X_in = [x0; y0; v0; theta0; 
        P11; P22; P33; P44; 
        P12; P13; P14; P23;
        P24; P34];
    
% extra array to collect all result from ode45
tot_T=[0];
tot_X=X_in';

for k = 1 : length(z_k)
    % prediction
    [t, solXs] = ode45(@dPredictFunc, [(k-1)*T k*T], X_in);
     solX = solXs(end,:)';
     Ex_k(:,k+1) = [solX(1);  solX(2); solX(3); wrapToPi(solX(4))];
     
     P_k(:,:,k+1) = [solX(5)  solX(9)  solX(10) solX(11); 
                     solX(9)  solX(6)  solX(12) solX(13);
                     solX(10) solX(12) solX(7)  solX(14);
                     solX(11) solX(13) solX(14) solX(8)];
    % observation
    H = [1 0 0 0;
         0 1 0 0];
    
    %update K and P(+)
    K = P_k(:,:,k+1) * H'* inv(H*P_k(:,:,k+1)*H'+R);
    P_k(:,:,k+1) = (eye(4) - K*H)*P_k(:,:,k+1);
    Ex_k(:,k+1) = Ex_k(:,k+1) + K*(z_k(k,:)' - [Ex_k(1,k+1); Ex_k(2,k+1)]);
    
    X_in = [Ex_k(1,k+1); Ex_k(2,k+1); Ex_k(3,k+1); wrapToPi(Ex_k(4,k+1)); 
            P_k(1,1,k+1); P_k(2,2,k+1); P_k(3,3,k+1); P_k(4,4,k+1);
            P_k(1,2,k+1); P_k(1,3,k+1); P_k(1,4,k+1); P_k(2,3,k+1);
            P_k(2,4,k+1); P_k(3,4,k+1)];
    
    % collecting all ode results
    tot_T = cat(1, tot_T, t(2:end));
    tot_X = cat(1, tot_X, solXs(2:end, :));
    if k == 45
        break;
    end
end

% ploting measurement x_m and the EKF estimated x(t)
t = (0: T: T* length(z_k));

hold on
plot(tot_X(:,1),tot_X(:,2), 'r')
figure
plot(tot_T,tot_X(:,3), 'g')
figure
plot(tot_T,tot_X(:,4), 'k')
