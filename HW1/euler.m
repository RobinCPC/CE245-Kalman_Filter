% euler method
% �ҡG�]f(t,y)=-10y�A���t=0, y0=2�C�By'=-10y
% ���u��Ѭ�y(t)=2e-10t�C�]�Gt��0.02�A�H���Q��MATLAB�D�ѡC

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
