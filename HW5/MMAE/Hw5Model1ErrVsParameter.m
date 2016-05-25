clear all
close all
load('robotMes.mat');
figure(1)
plot(y(1,:),'b.-')
hold on
plot(y(2,:),'g.-')
xlabel('Time')
legend('y_1=Distance to corner','y_2=Distance to wall')
axis([0 100 40 240]);
%part b
detInp=[0:0.01:4];
for iter=1:length(detInp),
  %Run KF for each parameter value
  dt=1;
  %Phi
  Phi=[1,1,0;
     0,0,0;
     0,0,1];
  %Deterministic input
  U=[0;
     detInp(iter);
     0]; 
  %Matrix multiplying the noise
  Gamma=[0;
     0.1;
     0];
  %Observation model
  H=[1,0,0;
     -1,0,1];
  %Measurement noise covariance
  R=[10 0; 0 10];
  %Initial guess, estimated from the first two data points of y1
  %xhat=[y(1,1),(y(1,2)-y(1,1))/dt,y(1,1)+y(2,1)]';
  %P=[R(1,1) 0 0;
  %  0 2*R(1,1)/dt^2 0;
  %  0 0 R(1,1)+R(2,2)];
  xhat=[y(1,1), detInp(iter),y(1,1)+y(2,1)]';
  P=[R(1,1) 0 0;
    0 0.01 0;
    0 0 R(1,1)+R(2,2)];
  %time vector
  t_Rec=[2:length(y)];
  clear xhat_Rec; 
  clear p11_Rec;
  clear p_33_Rec;
  xhat_Rec(:,1)=xhat; % <- 1 corresponds to t=2
  p11_Rec(1)=P(1,1);     % <- 1 corresponds to t=2
  p33_Rec(1)=P(3,3);     % <- 1 corresponds to t=2
  sumErr=0;
  sumErrW=0;
  sumLogErr=0;
  for t=3:length(y),
    xhatP=Phi*xhat+U;
    PP=Phi*P*Phi'+Gamma*Gamma';
    K=PP*H'*inv(H*PP*H'+R);
    ym=[y(1,t);
        y(2,t)];
    res=(ym-H*xhatP);
    xhat=xhatP+K*res;
    %Sum of weighted residuals
    sumErr=sumErr+res'*res;
    sumErrW=sumErrW+res'*inv(H*PP*H'+R)*res;
    sumLogErr=sumLogErr+res'*inv(H*PP*H'+R)*res+det(H*PP*H'+R);
    P=(eye(3,3)-K*H)*PP; 
    xhat_Rec(:,t-1)=xhat; % <- corresponds to t after update
    p11_Rec(t-1)=P(1,1);  % <- corresponds to t after update
    p33_Rec(t-1)=P(3,3);  % <- corresponds to t after update
  end 
  sumErr_Rec(iter)=sumErr;
  sumErrW_Rec(iter)=sumErrW;
  sumLogErr_Rec(iter)=sumLogErr;
end
figure(2)
plot(detInp,sumErr_Rec);
figure(3)
plot(detInp,sumErrW_Rec);
figure(4)
plot(detInp,sumLogErr_Rec);
