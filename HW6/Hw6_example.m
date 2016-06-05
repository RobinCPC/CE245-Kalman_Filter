clear all;
close all;

r1 = 1000;
r2 = 500;
var_r = 5;

x0 = 10000;
varx0 = 50;
dotx0 = -500;
vardotx0 = 200;
beta0 = 6e7;
varbeta0 = 2e12;

z = [ 9055, 8560, 7963, 7467, 7000, 6378, 5885, 5400, 4928, 4503];

P = [varx0 0 0; 0 vardotx0 0 ;0 0 varbeta0];

XRec = zeros(3,length(z));
Krec = zeros(3,1,10);
P11Rec=P(1,1);
P22Rec=P(2,2);
P33Rec=P(3,3);
R = 5;
X = [x0;dotx0;beta0];
xRec(:,1)=X;
tRec(1)=0;
for k=1:length(z),
  Sol0=[X(1,1);X(2,1);X(3,1);...
        P(1,1);P(1,2);P(1,3);P(2,2);P(2,3);P(3,3)];
  [tSol Sol]=ode45(@odeFun,[0 1],Sol0);
  Sol=Sol';
  tRec=[tRec; tSol(2:end)+(k-1)*ones((length(tSol)-1),1)];
  xRec=[xRec, Sol(1:3,2:end)];
  size(Sol)
  Xm=[Sol(1,end);Sol(2,end);Sol(3,end)];
  Pm=[Sol(4,end) Sol(5,end) Sol(6,end);
      Sol(5,end) Sol(7,end) Sol(8,end);
      Sol(6,end) Sol(8,end) Sol(9,end)];
  P11Rec=[P11Rec Sol(4,2:end)];
  P22Rec=[P22Rec Sol(7,2:end)];
  P33Rec=[P33Rec Sol(9,2:end)]; 
  h = sqrt(r1^2 + (Xm(1,1) - r2)^2);
  H = [(Xm(1,1) - r2)/(sqrt(r1^2 + (Xm(1,1) - r2)^2)) 0 0];
  K= Pm * H' * inv(H * Pm * H' + R);
  X= Xm + K * (z(k) - h)
  P = (eye(3) - K*H) *  Pm ;
  clear tSol;
  clear Sol;
end
tRec=[tRec;tRec(end)];
xRec=[xRec, X];
P11Rec=[P11Rec P(1,1)];
P22Rec=[P22Rec P(2,2)];
P33Rec=[P33Rec P(3,3)]; 

figure(1)
plot(tRec,xRec(1,:)), hold on
for k=1:10,
  xm(k)=sqrt(z(k)^2-r1^2)+r2
end    
plot([1:10],xm,'ro'), hold on
ylabel('xm'); 
figure(2)
plot(tRec,xRec(2,:)), hold on
ylabel('Velocity');
figure(3)
plot(tRec,xRec(3,:)), hold on
ylabel('Ballistic coefficient  \beta');
figure(4)
plot(tRec,P11Rec), hold on
ylabel('P11');
figure(5)
plot(tRec,P22Rec), hold on
ylabel('P22');
figure(6)
plot(tRec,P33Rec), hold on
ylabel('P33');


