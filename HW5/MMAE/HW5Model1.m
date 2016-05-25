clear all
close all
load('robotMes.mat');
%part b
dt=1;
%Phi
Phi=[1,1,0;
   0,0,0;
   0,0,1];
%Deterministic input
U=[0;
   0.8;
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
xhat=[y(1,1),0.8,y(1,1)+y(2,1)]';
%P=[R(1,1) 0 0;
%  0 2*R(1,1)/dt^2 0;
%  0 0 R(1,1)+R(2,2)];
P=[R(1,1) 0 0;
  0 0.01 0;
  0 0 R(1,1)+R(2,2)];

%time vector
t_Rec=[2:length(y)];
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
  res1_Rec(t-2)=res(1);
  res2_Rec(t-2)=res(2);
  xhat=xhatP+K*res;
  sumErr=sumErr+res'*res; % Sum of res^1
  sumErrW=sumErrW+res'*inv(H*PP*H'+R)*res; %Weighted sume of res^2
  sumLogErr=sumLogErr+res'*inv(H*PP*H'+R)*res+det(H*PP*H'+R); %log
  P=(eye(3,3)-K*H)*PP; 
  xhat_Rec(:,t-1)=xhat; % <- corresponds to t after update
  p11_Rec(t-1)=P(1,1);     %<- corresponds to t after update
  p33_Rec(t-1)=P(3,3);     %<- corresponds to t after update
end
MeanSqError=sumErr/length(res1_Rec)
MeanSqWError=sumErrW/length(res1_Rec)
MeanLogError=sumLogErr/length(res1_Rec)
%mean(res1_Rec)
%mean(res2_Rec)
figure(1)
xcorr_res1=xcorr(res1_Rec,res1_Rec);
subplot(2,3,1),plot(xcorr_res1),title('r11'),ylabel('Model 1');
xcorr_res2=xcorr(res2_Rec,res2_Rec);
subplot(2,3,2),plot(xcorr_res2),title('r22')
xcorr_res12=xcorr(res1_Rec,res2_Rec);
subplot(2,3,3),plot(xcorr_res12),title('r12')