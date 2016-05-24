%% CMPE 245 HW2 2-b
% -------------------------
% CMPE 245 - Spring, 2016
% Chien-Pin Chen
% 04/13/2016
% MATLAB script for HW2 b
% reference of matlab code:
% https://www.youtube.com/watch?v=LZ0qjZezGkQ&nohtml5=False
% --------------------------
clc;
%clf;

% set time vector
t = linspace(-2, 2, 8000);

% construct the uniform-shaped p_n(n)
x = 0*t;
ind = t >= -1/4 & t <=1/4;
x(ind) = 1/(1/2);

% plot uniform-shaped plot
figure;
plot(t, x, 'linewidth', 2);
ylim([-0.25 2.25]);
set(gca, 'fontsize', 12);
grid on;
xlabel('n');
ylabel('p_n(n)');
title('Uniform Function p_n(n)');


% compute the convolution 
z = median(diff(t))*conv(x,x,'same');

% plot the matlab convolution
figure;
plot(t, z,'linewidth', 2);
ylim([-0.25 2.25]);
set(gca, 'fontsize', 12);
grid on;
xlabel('m_k');
ylabel('p(m_k)');
title('Distribution Function of p(m_k)');

% get p_m probability distribution
p_m = z(3000:5000);
cdf_pm = (0);
for i=2:length(p_m)
    cdf_pm(i) = cdf_pm(i-1)+p_m(i);
end
cdf_pm = cdf_pm./length(cdf_pm);

% plot cumulated distribution function (cdf) of p_m
figure;
plot(t(3000:5000), cdf_pm,'linewidth', 2);
xlim([-2.0 2.0]);
ylim([-0.25 1.25]);
set(gca, 'fontsize', 12);
grid on;
xlabel('m_k');
ylabel('cdf of p(m_k)');
title('Cumulated Distribution Function of p(m_k)');

% simulate the stochastic process of problem 1.
m_k = linspace(-0.5, 0.5, length(p_m));
x=0;
for i=1:100
    x = x + m_k(randi(length(m_k)));
end

x_k = x