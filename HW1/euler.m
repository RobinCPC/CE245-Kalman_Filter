% euler method
% 例：設f(t,y)=-10y，初值t=0, y0=2。且y'=-10y
% 之真實解為y(t)=2e-10t。設Δt＝0.02，以此利用MATLAB求解。

function y = euler(y0, time, delta)
% Using Euler method to solve differential eqs.
%   y0: initial value
%   time: time limit
%   delta: step size
% Example: y = euler(2, 0.5, 0.02)

r= -10; k=0;
t = 0:delta:time; y = zeros(size(t));
y(1) = y0;

for i=2:length(t)
    y(i) = y(i-1) + r*y(i-1)*delta;
end

y_true = 2*exp(-10*t);
plot(t, y, 'o', t, y_true), xlabel('t'), ylabel('y');
