%% HW3-1
% -------------------------
% CMPE 245 - Spring, 2016
% Chien-Pin Chen
% 04/22/2016
% MATLAB script for HW3
% --------------------------
clear all;

% setting initial value
ph = [1.5 1; -0.7 0];
ga = [1.0; 0.5];

% vector to stoarge x(k) for k = 1..500
x_k = zeros(2, 500);
x_0 = [0;0];        % x(0)
w_k = random('Normal', 0, 1, 2,1);



% initalize figure and its handles for ploting line later
f1 = figure;        % delcare an empty canvas fb: figure for black line 
h_f1 = gca;         % get handle value (pointer) of above canvs
set(gca, 'fontsize', 12);
set(h_f1,'NextPlot', 'add');
xlabel ('k');
ylabel ('x_1 (k)');
title('Stochastic Differential Equation: x_1');

f2 = figure;
h_f2 = gca;
set(gca, 'fontsize', 12);
set(h_f2,'NextPlot', 'add');
xlabel ('k');
ylabel ('x_2 (k)');
title('Stochastic Differential Equation: x_2');

cmap = colormap(lines);    % plot 5 seqent in different colors.

% compute stochastic differential equation 500 times for 5 sequences
for seq = 1:5
    
    x_k(:,1) = ph*x_0 + ga* random('Normal', 0, 1);
    
    for k = 2:500
        x_k(:,k) = ph*x_k(:,k-1) + ga* random('Normal', 0, 1);
    end
    t = 1:1:500;
    
    pl1 = plot (h_f1, t, x_k(1,:), 'color', cmap(seq,:));
    %set(pl1,'DisplayName',['seq' num2str(seq)]);
    pl2 = plot (h_f2, t, x_k(2,:), 'color', cmap(seq,:));
    %set(pl2,'DisplayName',['seq' num2str(seq)]);
end
%legend(h_f1,'show');
%legend(h_f2,'show');

% HW3-1c
% vector to stoarge x(k) for k = 1..500 and 50 tracjetories
x_k = zeros(2, 500, 50);

% compute equation 1 for 50 times
for seq= 1:50
    x_k(:,1, seq) = ph*x_0 + ga* random('Normal', 0, 1, 1,1);
    
    for k = 2:500
        x_k(:,k, seq) = ph*x_k(:,k-1) + ga* random('Normal', 0, 1, 1,1);
    end
end

% standard deviation of x1(k) and x2(k)
x1= reshape(x_k(1,:,:), 1, []);
std_x1 = std(x1);
x2= reshape(x_k(2,:,:), 1, []);
std_x2 = std(x2);

% compute covariance of x1(k) and x2(k)
P_k_num = cov([x1' x2']);

% HW3-1d
% set initial value of P and Q
P_0 = [0 0; 0 0];
rQr = ga*1*ga';

% use eqaution 6 to compute P_k
P_k = zeros(2,2);
for k = 1:500
    P_k = ph*P_k*ph'+ rQr;
end

std_x1k = sqrt(P_k(1,1));
std_x2k = sqrt(P_k(2,2));



