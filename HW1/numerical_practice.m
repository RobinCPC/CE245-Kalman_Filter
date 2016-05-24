% numerical method practice
clear all % clc
%% test diff
h = 0.01;
x = 0:h:pi;
d = diff(sin(x.^2))/h;
plot(x, sin(x.^2), x, [0 d], 'r:')

x = (1:10).^2
x1 = diff(x), x2=diff(x,2)%, x3=diff(y)

%% test diff 2

t = [0:pi/50:pi]; nn = length(t); td=cos(t);
y =sin(t)+0.05*(randn(1,nn)-0.5);
dt1 = diff(y)./diff(t); % backward(+)
dt2 = (y(3:nn)-y(1:nn-2))./(t(3:nn)-t(1:nn-2));

plot(t,sin(t)); hold on; plot(t,y,'b.');
plot(t,cos(t),'r:');
plot(t(1:nn-1), dt1, 'k+');
plot(t(1:nn-2), dt2, 'bo');

[r1,p1]=corrcoef(t(1:nn-1),dt1)
[r2,p2]=corrcoef(t(1:nn-2),dt2)

%%  solving method of 1st-order differential equation

% euler method
% 例：設f(t,y)=-10y，初值t=0, y0=2。且y'=-10y
% 之真實解為y(t)=2e-10t。設Δt＝0.02，以此利用MATLAB求解。

y = euler(2, 0.5, 0.02);

% use ode45 to solve the same question. First, need to set function
figure;
fun = @(t,y)(-10*y);
[t,y] = ode45(fun, [0 0.5], 2);
plot(t,y,'ro', t, 2*exp(-10*t)), xlabel('t'),ylabel('y');


%%  non-linear equation solver
% example of sovling 2nd order
% example detial please see:
% http://bime-matlab.blogspot.com/search/label/Chap12
% Ch12.8 

clear all
% solve: 5y"+7y'+4y = f(t)

% rewrite y"=(1/5)f(t)-(4/5)y-(7/5)y'
% rewrite to fit 1st order form
% x1 = y
% x2 = y'

% diff again
% x1'= y'
% x2'= y"

% final two equantions:
% x1'= x2
% x2'= y" =(1/5)f(t)-(4/5)x1-(7/5)x2

fun = @(t,x) [x(2); (1/5)*sin(t)-(4/5)*x(1)-(7/5)*x(2)];
[t,x] = ode23(fun, [0 6], [3 9]);
plot(t, x(:,1),'r-', t, x(:,2),'bo', t,x(:,2));

%%  HW1 a test

% mLx'' = -mgsin(x)
% x'' = -(g/L)sin(x)

% x1' = x2
% x2' = x1'' = -(g/L)*sin(x1)   , where (g/L)=0.1
% I.C.: x1 = pi/6;   x2 = 0

fun = @(t,x) [x(2); -0.1*sin(x(1))];
[t,x] = ode23(fun, [0 100], [pi/6 0]);
plot(t, x(:,1),'r')
figure
plot(t, x(:,2))
figure
plot(t, x(:,1),'r-', t, x(:,2),'b', t,x(:,2));

%% HW1 b

% x1' = x2
% x2' = x1'' = -(g/L)*sin(x1) + f(t)  , where (g/L)=0.1
clear all;

T = 50;    % Total time interval
la = 0.02;  % parameter of exponential distribution

fun = @(t,x) [x(2); -0.1*sin(x(1))];
% get # of event happenned by p(t_k - t_k-1)

t_p = (0);
k = 1;
while (t_p(end) < T)
    k = k + 1;
    t_p(end + 1) = t_p(end) + random('exp', la);
end
n_event = k;
t_p(end) = T;

x0 = [0 0];     % Initial Condition x1=0 x2=0
t=[];
x=[];

for k=1:(n_event-1)
    [ts, xs] = ode23(fun, [t_p(k) t_p(k+1)], x0);
    t = [t:ts];
    x = [x:xs];
    ts=0,xs=0;
    f = random('norm', 0, 0.2);
    x0 = x(end) + [0 f];
end

