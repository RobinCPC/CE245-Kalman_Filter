%% HW3-2
% -------------------------
% CMPE 245 - Spring, 2016
% Chien-Pin Chen
% 04/22/2016
% MATLAB script for HW3
% --------------------------
clear all;

% setting initial value
x_k = zeros(1, 100, 100);
y_k = zeros(1, 100, 100);
r_k = zeros(1, 100, 100);
r_avg = zeros(1, 100);

dt = 0.1;   % 0.1s
T = [1:1:100]*dt;


% initalize figure and its handles for ploting line later
f1 = figure;        % delcare an empty canvas fb: figure for black line 
h_f1 = gca;         % get handle value (pointer) of above canvs
set(gca, 'fontsize', 12, 'FontName','Times New Roman');
set(h_f1,'NextPlot', 'add');
xlabel ('dx');
ylabel ('dy');
title('2D random walk');

f2 = figure;
h_f2 = gca;
set(gca, 'fontsize', 12, 'FontName','Times New Roman');
set(h_f2,'NextPlot', 'add');
xlabel ('Time (s)');
ylabel ('r(k\DeltaT)');
title('Time v.s. r(k\DeltaT)');


% compute stochastic differential equation in 10 sec (k 1:100) for 100 trajectories (traj 1:100)
for traj=1:100
    
    x_k(1,1,traj) = 0 + random('Normal', 0, dt*1);
    y_k(1,1,traj) = 0 + random('Normal', 0, dt*1);
    r_k(1,1,traj) = sqrt(x_k(1)^2+y_k(1)^2);
    
    for k = 2:100
        x_k(1,k,traj) = x_k(1, k-1, traj) + random('Normal', 0, dt*k);
        y_k(1,k,traj) = y_k(1, k-1, traj) + random('Normal', 0, dt*k);

        r_k(1,k,traj) = sqrt(x_k(1,k,traj)^2+y_k(1,k,traj)^2);
    end

    % plot each trajectory of x-y coordination and r(kdt)
    plot(h_f1, x_k(1,:,traj), y_k(1,:,traj));
    plot(h_f2, T, r_k(1,:,traj),'g');

end

% compute average of r(kdt) in each time k
for k = 1:100
    r_avg(k) = sum(r_k(1,k,:))/100;
end

plot(h_f2, T, r_avg,'b');

%plotyy(h_f2, T,r_k, T, r_avg);
%axis(h_f1, 'tight');
%axis(h_f2, 'tight');
