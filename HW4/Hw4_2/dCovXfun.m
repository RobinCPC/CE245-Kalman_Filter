%% CMPE 245 HW4-2
% -------------------------
% CMPE 245 - Spring, 2016
% Written by Prof. Dejan Milutinovic
% Edited by Chien-Pin Chen
% 05/08/2016
% MATLAB subfunction for HW4-2
% --------------------------
function dCovX=dCovXfun(t,X)
  P(1,1)=X(1,1);    % the differential of the variance of s(t)
  P(1,2)=X(2,1);    % the differential of the covariance of s(t) and v(t)
  P(2,1)=P(1,2);
  P(2,2)=X(3,1);    % the differential of the variance of v(t)
  A=[0 1; 0 -2];
  Gamma=[0; 20];
  dP=A*P+P*A'+Gamma*Gamma'; % the differential equation describeing dynamics of the variance of s(t) and v(t)
  dCovX=[dP(1,1);dP(1,2);dP(2,2)];  % return the result of the covriance matrix of s(t) and v(t)
end