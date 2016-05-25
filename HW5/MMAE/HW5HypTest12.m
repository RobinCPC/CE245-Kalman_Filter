clear all
close all
load('robotMes.mat');
%part b
dt=1;
%Phi
Phi1=[1,1,0;
   0,0,0;
   0,0,1];
Phi2=[1,1,0;
   0,1,0;
   0,0,1];
%Deterministic input
U1=[0;
   0.8;
   0];
U2=[0;
   0;
   0];
%Matrix multiplying the noise
Gamma1=[0;
   0.1;
   0];
Gamma2=[0;
   0.1;
   0];
%Observation model
H=[1,0,0;
   -1,0,1];
%Measurement noise covariance
R=[10 0; 0 10];
%Initial guess, estimated from the first two data points of y1
xhat1=[y(1,2),(y(1,2)-y(1,1))/dt,y(1,1)+y(2,1)]';
xhat2=[y(1,2),(y(1,2)-y(1,1))/dt,y(1,1)+y(2,1)]';
%xhat1=[y(1,2),0.8,y(1,1)+y(2,1)]';
%xhat2=[y(1,2),0.8,y(1,1)+y(2,1)]';
P1=[R(1,1) 0 0;
  0 2*R(1,1)/dt^2 0;
  0 0 R(1,1)+R(2,2)];
%P1=[R(1,1) 0 0;
%  0 0.01 0;
%  0 0 R(1,1)+R(2,2)];
P2=P1;
%time vector
t_Rec=[2:length(y)];
xhat1_Rec(:,1)=xhat1; % <- 1 corresponds to t=2
xhat2_Rec(:,1)=xhat2; % <- 1 corresponds to t=2
Pr1=0.5;
Pr2=1-Pr1;
Pr1_Rec(:,1)=Pr1;
Pr2_Rec(:,1)=Pr2;
for t=3:length(y),
  xhat1P=Phi1*xhat1+U1;
  xhat2P=Phi2*xhat2+U2;
  P1P=Phi1*P1*Phi1'+Gamma1*Gamma1';
  P2P=Phi2*P2*Phi2'+Gamma2*Gamma2';
  K1=P1P*H'*inv(H*P1P*H'+R);
  K2=P2P*H'*inv(H*P2P*H'+R);
  ym=[y(1,t);
      y(2,t)];
  res1=(ym-H*xhat1P);
  res2=(ym-H*xhat2P);
  A1=H*P1P*H'+R;
  A2=H*P2P*H'+R;
  beta1=1/((2*pi)*sqrt(det(A1)));
  beta2=1/((2*pi)*sqrt(det(A2)));
  prob1=beta1*exp(-0.5*res1'*inv(A1)*res1);
  prob2=beta2*exp(-0.5*res2'*inv(A2)*res2);
  Pr1U=prob1/(prob1*Pr1+prob2*Pr2)*Pr1;
  Pr2U=prob2/(prob1*Pr1+prob2*Pr2)*Pr2;
  Pr1_Rec(:,t-1)=Pr1U;
  Pr2_Rec(:,t-1)=Pr2U;
  Pr1=Pr1U;
  Pr2=Pr2U;
  xhat1=xhat1P+K1*res1;
  xhat2=xhat2P+K2*res2;
  P1=(eye(3,3)-K1*H)*P1P; 
  P2=(eye(3,3)-K2*H)*P2P;
end
figure(1)
plot(t_Rec, Pr1_Rec','b.-')
hold on
plot(t_Rec, Pr2_Rec','r.-');
legend('Model 1 probability','Model 2 probability');
